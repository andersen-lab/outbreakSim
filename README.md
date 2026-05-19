# Covasim Simulation Project

## Overview

This project uses **Covasim** to simulate the spread of COVID-19 within a defined region and estimate viral shedding into wastewater. The region is modeled as a grid (default: 2×2), allowing infections and viral load to be tracked spatially over time.

---

## Installation

Install Covasim via pip:

```bash
pip install -r requirements.txt
```

---

## Project Structure

* `src/run_sim.py` — Main script for running the simulation
* `src/plots` — Viral loads visualized in a heatmaps for example test cases.

---
## How to run the code

```python
python src/run_sim.py --pop_size 100 --pop_type hybrid --n_days 180 --location zambia
```

## Simulation Workflow

The simulation consists of several key steps:

### 1. Define Simulation Parameters

Simulation parameters are configured in the `define_sim_parameters()` function.

Core parameters include:

* `pop_size` — Total population size
* `pop_type` — Population structure (e.g., random, hybrid)
* `location` — Country/region name (used for demographic data)
* `n_days` — Number of simulation days
* `pop_infected` — Initial number of infected individuals
* `n_imports` — Daily number of imported infections

Additional spatial parameters:

* `n_rows` — Number of rows in the region grid (default: 2)
* `n_cols` — Number of columns in the region grid (default: 2)

Infection initialization parameters:

* `n_init_inf` — Number of initial infections in a region
* `r_init_inf` — Row index of initial infections
* `c_init_inf` — Column index of initial infections

Full parameter list:

```python
python src/run_sim.py --help

--pop_size POP_SIZE
--pop_type POP_TYPE
--n_days N_DAYS
--location LOCATION
--pop_infected POP_INFECTED
--n_imports N_IMPORTS
--n_rows N_ROWS
--n_cols N_COLS
--n_init_inf N_INIT_INF
--r_init_inf R_INIT_INF
--c_init_inf C_INIT_INF
```

**Note:**
Covasim automatically includes demographic data (age distribution and household sizes) based on the specified location.

* The function `check_age_household_dist()` validates whether data exists for the given location.
* Errors from this function typically indicate:

  * A typo in the location name, or
  * Missing data for that region

For more details on Covasim parameters, refer to:
[Covasim toturial](https://docs.idmod.org/projects/covasim/en/latest/tutorials/tut_intro.html)

---

### 2. Assign Individuals to Regions

The function `assign_people()` randomly assigns individuals to regions in the grid.

* Uses NumPy’s shuffle for randomness
* Returns an array of length `N` (population size)
* Each value represents a region ID (0 to total_regions − 1)

---

### 3. Initialize Infections

Initial infections are seeded using:

* `n_init_inf`
* `r_init_inf`
* `c_init_inf`

These define:

* How many individuals are infected
* The specific region (by row and column) where infections begin

By default:

* `pop_infected = 0`
* `n_imports = 0`

---

### 4. Run Simulation & Track Infections

The simulation is executed to track infections over time.

The function `calculate_new_infections()`:

1. Determines when each individual becomes infected
2. Maps individuals to regions
3. Aggregates infection counts per:

   * Region
   * Time step

---

### 5. Viral Shedding Estimation

Viral shedding is computed for each region over time using two approaches:

---

#### A. Eclipse Phase Model (Simplified)

A basic model that assumes:

* Uniform viral shedding behavior across individuals
* No variation due to age or other characteristics

This model provides a simplified estimate based on predefined shedding parameters.

---

#### B. Covasim-Based Viral Load

Covasim internally calculates viral load per individual at each time step.

This project leverages that data to compute shedding:

* `get_viral_loads(sim, t_start, t_end)`
* `viral_shedding_covasim(sim, start, end)`

Method:

* Extract per-individual viral loads
* Sum loads across individuals within each region
* Aggregate over time

Reference implementation:
[Covasim original paper](https://github.com/amath-idm/covasim_methods_paper/blob/main/figs/viral_load_plot.py)

---

## Summary

This pipeline:

1. Configures a population and spatial grid
2. Simulates infection spread using Covasim
3. Tracks infections by region and time
4. Estimates wastewater viral shedding using:

   * A simplified model
   * Covasim’s internal viral load calculations
