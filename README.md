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


## Definitions and calculations

* Rescaling: 
- Parameters : `pop_scale`=1, `scaled_pop`= None, `rescale`=TRUE, `rescale_threshold`=0.05, `rescale_factor`=1.02
- Functions: `rescale(self)`
- Depends on: `self.people.false('naive')`= who can get infected?,`self['pop_size']` how many agents?
- What it does: In Covasim, **`pop_scale`** is how many real people each simulated agent represents, giving you the total **`scaled_pop`** size. When **`rescale`** is turned on, the model watches for when too many agents catch the virus or get vaccinated—once that number hits the **`rescale_threshold`**, it uses the **`rescale_factor`** to prune non-susceptible agents and bump up the scale, keeping the simulation fast without running out of memory.

## Basic disease transmission
- Parameters: `beta`=0.016,`n_import`=0 (there will be random number of infected cases drawn from a poisson distribution even if set to zero),`beta_dist`='neg_binomial', par1=1.0, par2=0.45, step=0.01,`viral_dist`=frac_time=0.3, load_ratio=2, high_cap=4, `asymp_factor`=1
- what is does: transmission is set and controlled using a combination of a baseline transmissibility parameter (beta), network layer weights (beta_layer), intrahost viral dynamics, and modifying factor