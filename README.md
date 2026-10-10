# Trends-and-multidimensional-health-expectancy

## Description
This is the key analysis code for "Diverging trends and drivers in multidimensional health expectancy in China, 2002–2022: a nationwide cohort study".

## Code files
1-01_models. Cox, baseline logistic regression and three-state MSM 
2-02_life_expectancy. ELECT marginal health expectancy, harmonisation and uncertainty summaries
3-03_scenarios. Symmetric decomposition; 100% and 80% attainment

## Requirements
We performed analysis using R (version 4.5.2). No specialised hardware, such as a GPU, is required by the analysis code. Computational requirements depend on the dataset size, the number of model fits and the simulation settings.

## Installation
Download or clone the GitHub repository. 
Install R version 4.5.2.
Install the R packages required by the analysis scripts, including msm, elect, survival and mice, together with any additional dependencies identified in the scripts.
Prepare the harmonised longitudinal input data according to the data layout described below.
Configure the model parameters and input paths before running the analysis scripts.
Installation time varies depending on the R environment, package availability and network speed.

## Demonstration
A minimal demonstration of the health expectancy estimation interface is provided in tests/smoke_elect.R.
The demonstration uses electData, a simulated example dataset distributed with the R package elect.
To run the demonstration: Rscript tests/smoke_elect.R
Run this command from the repository root after installing the required R packages.
Expected output: ELECT interface smoke test passed. Not a CLHLS reproduction test.
The script fits a continuous-time multistate Markov model and verifies the output dimensions, finite numerical estimates and reproducibility of the health expectancy calculation under a fixed random seed.
This demonstration is intended to validate the software interface rather than reproduce the numerical findings of the manuscript.

## Data availability
The data that support the findings of this study are available from Peking University Open Research Data. The data access requires the completion of a detailed application form from the website.
Researchers interested in this study can contact the corresponding author (Yao Yao) for further information.
﻿
## Data layout
One row is a living survey observation or a terminal death observation. There is exactly one baseline per person-period. A person appearing in both periods must keep the same `person_id`; the program constructs a separate trajectory identifier for each period and reports overlap. `period=0` is 2002–2011/2012; `period=1` is 2011/2012–2021/2022. Include only the relevant observation windows. Dates use ISO format. Age is `(date - dob)/365.25`; confirm that this is the original convention. Baseline age determines the 65–74, 75–84 and 85–94 strata. All five health states use 1=favourable, 2=adverse, 3=death. Missing values are
never coded as favourable. Terminal covariates are carried from the preceding living visit only; this does not impute unmeasured exposures after death.

## Technical references
- ELECT official manual: https://cran.r-project.org/web/packages/elect/elect.pdf
- MSM official reference: https://search.r-project.org/CRAN/refmans/msm/html/msm.html
- ELECT authors' methods paper: https://discovery.ucl.ac.uk/10075585/1/CMPB_AvdH_AcceptedVersion.pdf
- ELECT source documentation inspected for output names and random-draw order:
  https://rdrr.io/cran/elect/src/R/elect.r
