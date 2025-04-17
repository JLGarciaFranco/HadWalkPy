# Divergent & Non‑Divergent Wind Decomposition

This repository provides a Python implementation of the Helmholtz decomposition for atmospheric wind fields, along with computation of local Hadley and Walker circulations based on Schwendike et al. (2014). You can decompose any 3D wind dataset (time × pressure × latitude × longitude) into its divergent (irrotational) and rotational (non‑divergent) components, then compute associated mass fluxes, vertical velocities, and streamfunctions.

---
## Features

1. **Helmholtz Decomposition**  
   - Separate horizontal wind into divergent (irrotational) and rotational (non‑divergent) parts using `windspharm.xarray`.

2. **Local Hadley & Walker Circulations**  
   - Compute mass streamfunctions for meridional (Hadley) and zonal (Walker) cells via vertical integration.  
   - Calculate vertical mass flux and Omega (vertical velocity) diagnostics for each component.

3. **Flexible Data Handling**  
   - Accepts pressure levels in hPa or Pa (auto‑converted).  
   - Works with monthly climatologies, daily data, or arbitrary time series.  
   - Optional parallel processing for faster decomposition over many time steps.

4. **Single‑Class API**  
   - All functionality exposed through the `Hadley_Walker` class in `decompose_wind.py`.

---
## Installation & Requirements

Clone this repository and install dependencies via `pip`:

```bash
git clone https://github.com/your-org/HadleyWalkerDecomp.git
cd HadleyWalkerDecomp
pip install -r requirements.txt


