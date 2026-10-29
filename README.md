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

### Population Setup & Generation

* **Key Parameters:** `pop_size = 20e3`, `pop_infected = 20`, `pop_type = random`, `location = None (seattle)`
* **Core Functions:** `make_people`, `assign_sewersheds`, `make_synthpop`, `make_randpop`
* **Process:** These functions create people depending on the population type provided. Before generating people, the system checks if the user provided population data.
* **Synthpop:** Population type `synthpop` uses [synthpops GitHub repository](https://github.com/synthpops/synthpops). This only supports a few locations and is not being maintained anymore, but it outputs a structured synthetic population with realistic demographic details—such as age, sex, and household size—mapped onto realistic interaction networks across households, schools, workplaces, long-term care facilities, and general community spaces. Instead of generating random contacts, it uses real census data to simulate how people actually group and interact in everyday life, providing the underlying contact network that disease models use to track transmission.

---

### Contact Layers & Person Attributes

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

## Contacts

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

## Rescaling

* **Parameters:** `pop_scale = 1`, `scaled_pop = None`, `rescale = TRUE`, `rescale_threshold = 0.05`, `rescale_factor = 1.02`
* **Function:** `rescale(self)`
* **Dependencies:** `self.people.false('naive')` (who can get infected?), `self['pop_size']` (how many agents?)
* **What it does:** In Covasim, **`pop_scale`** determines how many real people each simulated agent represents, giving you the total **`scaled_pop`** size. When **`rescale`** is turned on, the model watches for when too many agents catch the virus or get vaccinated—once that number hits the **`rescale_threshold`**, it uses the **`rescale_factor`** to prune non-susceptible agents and bump up the scale, keeping the simulation fast without running out of memory.

---

## Basic Disease Transmission

* **Parameters:** `beta = 0.016`, `n_import = 0` (there will be a random number of infected cases drawn from a Poisson distribution even if set to zero), `beta_dist = 'neg_binomial'`, `par1 = 1.0`, `par2 = 0.45`, `step = 0.01`, `viral_dist = frac_time = 0.3`, `load_ratio = 2`, `high_cap = 4`, `asymp_factor = 1`
* **What it does:** Transmission is set and controlled using a combination of a baseline transmissibility parameter (`beta`), network layer weights (`beta_layer`), intrahost viral dynamics, and modifying factors.