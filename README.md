## Input data

Before running the optimization script, two types of input data should be prepared:

1. **Occurrence data**

   Species occurrence records should first be spatially filtered to reduce
   sampling redundancy and spatial autocorrelation. In this study, occurrence
   records were spatially filtered using ENMTools, and the retained coordinates
   were used as input for subsequent MaxEnt optimization.

   The filtered occurrence data are provided in:
   `occurrence_filtered.csv`

2. **Environmental variables**

   Environmental variables should be screened before model optimization.
   In this study, environmental predictors were selected based on pairwise
   correlation analysis together with their contribution to the preliminary
   MaxEnt model. Highly correlated variables were removed to reduce
   multicollinearity, and the retained environmental variables were used for
   subsequent ENMeval analysis.

   The final set of environmental variables used for model optimization is
   provided in:
   `environmental_variable_selection.csv`

   Environmental raster layers are not included in this repository because of
   their file size. They can be obtained from the corresponding public
   environmental databases and prepared using the same spatial extent,
   resolution, and coordinate reference system before running the analysis.

# Asarum sieboldii MaxEnt Optimization

This repository contains the R scripts and associated results used for
ENMeval-based optimization of the MaxEnt model for *Asarum sieboldii*.

## Model optimization

Six feature-class combinations (L, LQ, H, LQH, LQHP, and LQHPT) and
eight regularization multiplier (RM) values ranging from 0.5 to 4.0
were evaluated, resulting in 48 candidate models.

Model selection was primarily based on AICc and ΔAICc. The model with
the lowest AICc (ΔAICc = 0) was selected as the optimal model.

The optimal parameter combination was:

- Feature classes (FC): LQHPT
- Regularization multiplier (RM): 3.5

## Files

- `maxent_optimization.R` – R script for MaxEnt parameter optimization
- `all_candidate_models.csv` – performance metrics of all 48 candidate models
- `selected_optimal_model.csv` – results for the selected optimal model
- `default_vs_optimized_full.csv` – comparison between default and optimized models
- `default_vs_optimized_change_summary.csv` – summary of changes after optimization
- `delta_AICc.pdf` – ΔAICc values across candidate parameter combinations

## Software

Model optimization was performed using the ENMeval package in R.
