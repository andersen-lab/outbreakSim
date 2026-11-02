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

# Definitions and Calculations

---

## 1. Population Setup & Generation

* **Key Parameters:** `pop_size = 20e3`, `pop_infected = 20`, `pop_type = random`, `location = None (seattle)`
* **Core Functions:** `make_people`, `assign_sewersheds`, `make_synthpop`, `make_randpop`
* **Process:** These functions create people depending on the population type provided. Before generating people, the system checks if the user provided population data.
* **Synthpop:** Population type `synthpop` uses the [synthpops GitHub repository](https://github.com/synthpops/synthpops). This only supports a few locations and is not being maintained anymore, but it outputs a structured synthetic population with realistic demographic details—such as age, sex, and household size—mapped onto realistic interaction networks across households, schools, workplaces, long-term care facilities, and general community spaces. Instead of generating random contacts, it uses real census data to simulate how people actually group and interact in everyday life, providing the underlying contact network that disease models use to track transmission.

---

## 2. Contact Layers & Person Attributes

* **Location Effects:** By default, when a location is specified, age and household size data are also used. If location is not specified, then Seattle data is used.
* **Contact Layers:** There are 4 key contact types in the population:
* Household = `h`
* Schools = `s`
* Workplace = `w`
* Random contact = `c`


* **Note:** `layers_defaults` holds these values for each contact type.
* **Person Attributes:** Once created, each individual contains: `uid`, `age`, `sex`, and `contacts`.
* **Prognoses & Spatial Data:** Prognoses are set by default for different age bins, along with regions and coordinates. Once people are created, they are also assigned a sewershed based on their coordinates.

---

## 3. Contacts

* **Core Functions:** `make_random_contacts`, `make_hybrid_contacts`, `make_rank_based_contacts`
* **`make_random_contacts`:** Every person needs on average $N$ contacts. By default, it uses a Poisson distribution centered at $N$ to decide how many contacts each person gets. If a dispersion factor is used, it switches to a negative binomial distribution to account for superspreading events. Dispersion is `None` by default.
* **`make_hybrid_contacts`:** Makes networks across different social layers. Exact contact information is set by the `reset_layer_pars` function per `pop_type` (`h`, `s`, `w`, `c`), with school and work ages defined as `[6, 22]` and `[22, 65]` respectively.

### Hybrid Contact Structure

* It starts off with making contacts in the household (`make_microstructured_contacts` = clustering people into families with a cluster size drawn from a poisson distribution), then making random contacts between everyone (`make_random_contacts`).
* Then it identifies the individuals that fall into work age and school age. From there, it creates random contacts for each of these groups. So household contacts are structured but not school or work places.

### Rank Contact Structure

* If the `pop_type` rank is used, it first makes sure that people have coordinates, if they do, we will convert lonlat to km using `lonlat_to_km`. If not, it will assign them coordinates using `assign_density_based_coords`.
* Then the rank based contacts are calculated with a mean set as default `pars['contacts']` for "random" which is 20.
* Given the current structure of the code, it does not look like dispersion can be defined as a parameter. The default for dispersion is none, in both `make_rand_contact` and `rank_based_contact`.
* The logic is assigns coordinates based on a density raster file which is provided by the user.

### Dynamic Layers Note

* It is important to note that we can also set layers to be dynamic, meaning that in the specific layer such as "h" the contacts stay the same throughout the simulation.
* Default values for both random and hybrid is set to zero, meaning that in all layers the contacts stay the same and not change throughout the time.

---

## 4. Rescaling

* **Parameters:** `pop_scale = 1`, `scaled_pop = None`, `rescale = TRUE`, `rescale_threshold = 0.05`, `rescale_factor = 1.02`
* **Function:** `rescale(self)`
* **Dependencies:** `self.people.false('naive')` (who can get infected?), `self['pop_size']` (how many agents?)
* **What it does:** In Covasim, **`pop_scale`** determines how many real people each simulated agent represents, giving you the total **`scaled_pop`** size. When **`rescale`** is turned on, the model watches for when too many agents catch the virus or get vaccinated—once that number hits the **`rescale_threshold`**, it uses the **`rescale_factor`** to prune non-susceptible agents and bump up the scale, keeping the simulation fast without running out of memory.

---

## 5. Basic Disease Transmission

* **Parameters:** `beta = 0.016`, `n_import = 0` (there will be a random number of infected cases drawn from a Poisson distribution even if set to zero), `beta_dist = 'neg_binomial'`, `par1 = 1.0`, `par2 = 0.45`, `step = 0.01`, `viral_dist = frac_time = 0.3`, `load_ratio = 2`, `high_cap = 4`, `asymp_factor = 1`
* **What it does:** The probability of transmission from an *infectious* person to a *susceptible* person, also known as infectiousness or transmissibility. The overall transmission probability depends on the network *layer* two people are connected by, as well as the infected person's *viral load*, the susceptible person's age, and other factors.

### Viral Load Calculation (`compute_viral_load`)

* The `compute_viral_load` function dynamically calculates the relative transmissibility of infected individuals at any given timestep `t` by modeling how viral load evolves over the course of an infection. It first determines each person's infection duration by evaluating their start date against their endpoint (either recovery or death), and then divides this timeline into an initial high-load phase and a subsequent low-load phase, subject to a maximum cap (`high_cap`) to prevent unrealistically prolonged peak infectiousness. Finally, the function scales the viral load values based on whether the individual is currently in the early or late phase of infection, automatically setting the load to zero for anyone who has not yet been infected, has already recovered or died, or was never infected at all.

### Transmission Parameter Breakdown

* **`rel_trans` Initial Draw:** `rel_trans` is first drawn from the `beta_distribution` and multiplied by the odds ratio of transmissibility by each individual which currently is 1 for all age cutoff. `beta_dist` default is `dist='neg_binomial', par1=1.0, par2=0.45, step=0.01`.
* **Sequence Tracker:** The `seq_tracker` will overwrite the `rel_trans` to `rel_trans` multiplied by the fitness value calculated. *(Question noted in notes: is this supposed to happen?)*
* **`beta`:** A fixed value set by default to `0.016`, representing beta per symptomatic contact.
* **`beta_distribution`:** Defined as `'neg_binomial'`, `par1=1.0`, `par2=0.45`, `step=0.01`. This is used to calculate the initial relative transmissibility by multiplying the `beta_distribution` by the odds ratio of transmissibility for each individual, which is currently 1 for all age cutoffs (10 bins).

### Important Modifiers and Variant Adjustments

* **Sequence Tracker:** If the sequence tracker is on, `rel_trans` is multiplied by the fitness value.
* **Variant Recalculation:** For each variant, `beta` is recalculated using the formula:
`self['beta'] * self['rel_beta'] * self['variant_pars'][variant_label]['rel_beta']`
Here, `beta` is `0.016` (default), and `rel_beta` is defined per variant (e.g., alpha, delta).
* **Layer-Specific Betas:** For each pair in a contact, betas are defined as the layer's beta (e.g., households, schools, and workplaces have their own specific beta values).

---

## 6. Final Calculations (`compute_trans_sus`)

The final relative transmissibility and susceptibility are calculated through the function:

`compute_trans_sus(prel_trans, prel_sus, inf_variant, sus, beta_layer, viral_load, symp, iso, quar, asymp_factor, iso_factor, quar_factor, sus_imm)`

### Key Calculation Steps Include:

* `f_asymp = symp + ~symp * asymp_factor`
*(Asymptomatic factor; changes e.g. `[0,1]` with a factor of `0.8` to `[0.8, 1.0]`)*
* `f_iso = ~iso + iso * iso_factor`
*(Isolation factor; changes e.g. `[0,1]` with a factor of `0.2` to `[1, 0.2]`)*
* `f_quar = ~quar + quar * quar_factor`
*(Quarantine factor; changes e.g. `[0,1]` with a factor of `0.5` to `[1, 0.5]`)*
* `rel_trans = rel_trans * inf * f_quar * f_asymp * f_iso * beta_layer * viral_load`
*(Recalculate transmissibility)*
* `rel_sus = rel_sus * sus * f_quar * (1-immunity_factors)`
*(Recalculate susceptibility)*