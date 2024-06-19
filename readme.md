# Implementation and numerical experiments for Splicing-SIR

This repository contains codes for simulation and real data analysis described in *A Splicing Algorithm for Best Subset Selection in Sliced Inverse Regression*.

## Codes

- `simulation.R`: Code for simulation.
- `analysis_realdata.R`: Code for real-world data analysis.
- `sir_splicing.R`: Implementation of Splicing-SIR.
- `visualization.R`: Code for visualization of simulation results.

## Methods

- Splicing-SIR: Implemented in `sir_splicing.R`.
- SEAS: Source codes may be downloaded from the [Supplemental Material of Subspace Estimation with Automatic Dimension and Variable Selection in Sufficient Dimension Reduction](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2022.2118601?scroll=top). 
  - Please put the scripts `utility.R` and `seas.R` in the `seas/` directory. 

- LASSO-SIR: Implemented in R package `LassoSIR` (0.1.1).

**Note:** Comparison with the [Relaxed Natural SSIR](https://doi.org/10.1214/18-AOS1791) method is omitted in this repository because the codes are not publicly available. They may be requested from the authors. 

## Results

Results of experiments and visualization will appear in the `results/` directory.
