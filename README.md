# spatial_stacking

Code for paper "Bayesian Geostatistics Using Predictive Stacking".

## Roadmap

| Folder Name | Intro |
|:--------------------------|:--------------------------------------------|
| RDA | Data and code for the AOD analysis |
| sim | Code for assumption checking, and example R code for simulation studies |
| sim_carc | Code submitted to cluster for simulation studies in Section 4 |
| sim_julia | Julia code for running time comparison in Section 4.3 |
| utils.R | R functions required in running projects |

## Instructions

The files `Example_code.Rmd` (and its corresponding HTML version `Example_code.html`) provide a concise example to demonstrate the proposed stacking algorithms. The code for all the experiments conducted in the simulation studies, including those using computationally intensive methods for comparison, is available in `sim_cv_experiment_lit.R` and `sim_cv_experiment_prefix.R`, located within the `sim` folder. It's important to note that the simulation studies mentioned in the paper were carried out on a cluster. The examples provided in `sim_cv_experiment_lit.R` and `sim_cv_experiment_prefix.R` were executed with various configurations and seeds, repeated 120 times for the whole simulation studies. The code used for cluster submissions, along with the scripts for generating summary plots, can be found in the `sim_carc` folder. Additionally, the data and code for the real-data analysis, specifically the Aerosol Optical Depth (AOD) analysis, are available in the `RDA` folder.

------------------------------------------------------------------------

### **Roadmap of Figures and Tables**

| **Folder** | **File** | **Figures** | **Tables** |
|--------------|--------------------|------------------------|--------------|
| **sim_carc** | `projects/summary.R` | 1, 2, S3-4, S13–S15 | \- |
|  | `projects/summary2.R` | S16 | \- |
| **sim** | `projects/check_trH22.R` | S2 | \- |
|  | `projects/check_z_var.R` | S1 | \- |
|  | `projects/sim_cv_experiment_prefix` | 3, 5, 6, S5–S12, S17–S21 | \- |
| **sim_julia** | `projects/time_compar.ipynb` | 4 | \- |
| **RDA** | `projects/LA_RDA.R` | 7, 8, S22-23 | 1 |
|  | `projects/H2Otest.R` | \- | 1 |

------------------------------------------------------------------------

## Licensing

-   Code © 2025, licensed under [BSD (3-clause)](https://opensource.org/licenses/BSD-3-Clause).

## Notes

-   You are welcome to contact me for bug reporting.
