# outbreakSim

Worked examples of simulating a SARS-CoV-2 outbreak in a real city and watching it through clinical and wastewater surveillance. This is a [Modjadji](https://modjadji.info) project built on [covasim4wastewater](https://github.com/andersen-lab/covasim4wastewater), the Andersen lab's wastewater extension of [Covasim](https://docs.covasim.org).

## Examples

Each example is a script in `src/simulation_examples/<place>/` that writes its figures next to itself:

```bash
python src/simulation_examples/sandiego/sandiego_simulation.py
python src/simulation_examples/lusakacity/lusakacity_simulation.py
```

Both follow the same steps: simulate an outbreak on a map, track viral shedding by sewershed, introduce variants with real genomes, simulate amplicon sequencing reads, compare clinical and wastewater surveillance, and draw the transmission phylogeny.

* **San Diego** (`sandiego/`): San Diego County, with clinical cases by ZIP code and wastewater by treatment-plant sewershed. Written as a step-by-step tutorial. `make_sewershed_map.py` and `clip_worldpop.py` build the sewershed and population-raster inputs.
* **Lusaka** (`lusakacity/`): Lusaka, Zambia, with the same workflow using administrative boundaries, sewershed catchments and a population raster.

To model another place, supply a boundary file, a sewershed file and a density raster to `cv.Sim` (`boundary_shapefile`, `sewershed_file`, `density_raster`).

## Installation

Install covasim4wastewater following the instructions in its repository. The scripts need that build of `covasim`, plus `geopandas`, `rasterio`, `matplotlib`, `pillow`, `baltic` and `bygul`.

## What covasim4wastewater adds

Covasim is an agent-based COVID-19 model. covasim4wastewater adds:

* spatial populations placed using a density raster and region boundaries
* sewershed assignment, so shedding can be summed per treatment plant
* per-agent viral shedding over time
* sequence evolution, with variants seeded from real genomes and transmissibility derived from their mutations
* clinical and wastewater samplers, and amplicon read simulation with Bygul
* an infection log for building transmission phylogenies


# Covasim Model Architecture & Calculations

## 1. Population Setup & Generation

* **Key Parameters:** `pop_size = 20e3`, `pop_infected = 20`, `pop_type = random`, `location = None (seattle)`
* **Core Functions:** `make_people`, `assign_sewersheds`, `make_synthpop`, `make_randpop`
* **Process:** Generates agents based on the selected population type after verifying user-provided data.
* **Synthpop (`pop_type = 'synthpop'`):**
* Powered by the [synthpops GitHub repository](https://github.com/synthpops/synthpops).
* *Note:* Limited location support and currently unmaintained.
* *Output:* Outputs a structured synthetic population with demographics (age, sex, household size) mapped onto realistic interaction networks across households, schools, workplaces, long-term care facilities, and community spaces using real census data.



---

## 2. Contact Layers & Person Attributes

* **Location Effects:** Uses age and household size data for the specified location (defaults to Seattle if unspecified).
* **Contact Layers:**
* Household = `h`
* Schools = `s`
* Workplace = `w`
* Random contact = `c`
* *(Stored globally in `layers_defaults`)*


* **Person Attributes:** Each individual contains: `uid`, `age`, `sex`, and `contacts`.
* **Prognoses & Spatial Data:** Default prognoses assigned by age bins alongside regions, coordinates, and sewershed assignments.

---

## 3. Contact Generation Strategies

* **Core Functions:** `make_random_contacts`, `make_hybrid_contacts`, `make_rank_based_contacts`

### A. Random Contacts (`make_random_contacts`)

* Assigns contacts around a mean value $N$ using a **Poisson distribution**.
* Switches to a **negative binomial distribution** if a dispersion factor is active (default is `None`) to simulate superspreading events.

### B. Hybrid Contacts (`make_hybrid_contacts`)

* Sets exact contact rules per layer via `reset_layer_pars` (`h`, `s`, `w`, `c`), with school age defined as `[6, 22]` and work age as `[22, 65]`.
* **Structure:**
1. Microstructured household clustering (cluster size drawn from a Poisson distribution).
2. Random connections made across the broader population.
3. Random contacts assigned specifically within school and work age brackets (household contacts are structured; school/work networks are random).



### C. Rank-Based Contacts (`make_rank_based_contacts`)

* Requires coordinates; converts longitude/latitude to kilometers via `lonlat_to_km` or generates them via `assign_density_based_coords` using a user-provided density raster file.
* Calculates rank-based contacts with a default mean of `20` (from `pars['contacts']`).
* Dispersion parameter defaults to `None`.

### D. Dynamic Layers Note

* Layers can be set to dynamic, but the default for both random and hybrid models is `0` (static contacts that do not change over simulation time).

---

## 4. Rescaling

* **Parameters:** `pop_scale = 1`, `scaled_pop = None`, `rescale = TRUE`, `rescale_threshold = 0.05`, `rescale_factor = 1.02`
* **Core Function:** `rescale(self)`
* **Dependencies:** `self.people.false('naive')`, `self['pop_size']`
* **Mechanism:**
* `pop_scale` defines how many real individuals each simulated agent represents, determining the total `scaled_pop` size.
* When `rescale` is enabled, the model monitors non-susceptible agents. Once infections/vaccinations hit the `rescale_threshold`, it applies the `rescale_factor` to prune agents and increase scale, optimizing memory and execution speed.



---

## 5. Disease Transmission Mechanics

* **Parameters:** `beta = 0.016`, `n_import = 0` (random Poisson-drawn cases may still occur), `beta_dist = 'neg_binomial'`, `par1 = 1.0`, `par2 = 0.45`, `step = 0.01`, `viral_dist = frac_time = 0.3`, `load_ratio = 2`, `high_cap = 4`, `asymp_factor = 1`
* **Core Concept:** Transmission probability from an infectious person to a susceptible person depends on network layer connectivity, viral load, susceptible age, and other modifiers.

### A. Viral Load Calculation (`compute_viral_load`)

* Dynamically calculates relative transmissibility at timestep $t$ based on infection duration (start date to recovery/death endpoint).
* Splits timelines into a high-load early phase and a low-load late phase, bounded by a maximum cap (`high_cap`) to prevent unrealistic peak lengths.
* Automatically sets viral load to `0` for uninfected, recovered, dead, or unaffected individuals.

### B. Transmission Parameters & Modifiers

* **Initial Draw:** `rel_trans` is drawn from `beta_distribution` and multiplied by individual age-bin odds ratios (default is `1`).
* **Sequence Tracker Override:** If active, `rel_trans` is multiplied by the calculated viral fitness value.
* **Base Beta (`beta`):** Fixed default of `0.016` (beta per symptomatic contact).
* **Variant Adjustments:** Recalculated per variant using:

$$\text{beta} \times \text{rel\_beta} \times \text{variant\_pars}[\text{variant\_label}][\text{rel\_beta}]$$

* **Layer-Specific Betas:** Individual beta values are assigned per network layer (household, school, workplace, community).

---

## 6. Final Transmissibility & Susceptibility (`compute_trans_sus`)

Calculates final relative transmission and susceptibility values using the following core multipliers:

* **Asymptomatic Factor:**

$$f_{\text{asymp}} = \text{symp} + \sim\text{symp} \times \text{asymp\_factor}$$

* **Isolation Factor:**

$$f_{\text{iso}} = \sim\text{iso} + \text{iso} \times \text{iso\_factor}$$

* **Quarantine Factor:**

$$f_{\text{quar}} = \sim\text{quar} + \text{quar} \times \text{quar\_factor}$$

* **Final Transmissibility:**

$$\text{rel\_trans} = \text{rel\_trans} \times \text{inf} \times f_{\text{quar}} \times f_{\text{asymp}} \times f_{\text{iso}} \times \text{beta\_layer} \times \text{viral\_load}$$

* **Final Susceptibility:**

$$\text{rel\_sus} = \text{rel\_sus} \times \text{sus} \times f_{\text{quar}} \times (1 - \text{immunity\_factors})$$

---

## 7. Population Trackers & Epidemiological Outputs

### A. Population State Metrics

* **Alive (`n\_alive`):** $\text{Total Population} - \text{Total Deaths}$
* **Naive (`n\_naive`):** $\text{Total Population} - \text{Deaths} - \text{Recovered} - \text{Exposed}$
* **Susceptible (`n\_susceptible`):** $\text{Alive} - \text{Exposed} - (\text{count\_recov} \times \text{Cumulative Recoveries})$
* **Pre-infectious (`n\_preinfectious`):** $\text{Exposed} - \text{Currently Infectious}$
* **Removed (`n\_removed`):** $(\text{count\_recov} \times \text{Cumulative Recoveries}) + \text{Cumulative Deaths}$
* **Prevalence:** $\frac{\text{Currently Exposed}}{\text{Alive}}$
* **Incidence:** $\frac{\text{New Infections}}{\text{Susceptible}}$
* **Fraction Vaccinated:** $\frac{\text{Vaccinated}}{\text{Alive}}$
* **By-Variant Metrics:** Variant-specific incidence and prevalence calculated via Einstein summation ($\text{einsum}$).

### B. Epidemic Dynamics

* **Doubling Time:** Calculated over a rolling window using exponential approximation:

$$\text{Days to Double} = \min\left( \text{window} \times \frac{\log(2)}{\log\left(\frac{\text{Infections}_t}{\text{Infections}_{t - \text{window}}}\right)}, \text{max\_doubling\_time} \right)$$

* **Effective Reproduction Number ($R_{\text{eff}}$):**
* **Daily Method:**



$$\text{Mean Infectious Duration} \times \frac{\text{New Local Infections}}{\text{Currently Infectious}}$$

*(Includes kernel smoothing for noise reduction).*

* **Transmission Log Method (`infectious` / `outcome`):** Tracks direct person-to-person infections across the transmission tree, averaged over a rolling window.
* **Generation Time (Serial Interval):**
* **True Interval:** Days elapsed from source exposure to target exposure.
* **Clinical Interval:** Days elapsed from source symptom onset to target symptom onset.