# Covasim Simulation Project

## Overview

This project uses **Covasim** to simulate the spread of COVID-19 within a defined region and estimate viral shedding into wastewater. The region is modeled as a grid (default: 2×2), allowing infections and viral load to be tracked spatially over time.

---

## Installation

Install Covasim via pip:

```bash
pip install covasim
```

---

## Project Structure

* `src/simulation_examples` example simulation runs in different regions


---
## How to run the code

```python
python src/run_sim.py --pop_size 100 --pop_type hybrid --n_days 180 --location zambia
```
