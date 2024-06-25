------------------------------------------------------------------------

> Code repository associated with the following manuscript

------------------------------------------------------------------------

## Overview

This repository contains the code for a simulation study investigating the effects of non-normality on estimating causal effects in path and multigroup structural equation models with ordered categorical variables. The study utilizes advanced vine copula methods to simulate multivariate distributions, to investigate the robustness of SEM estimators under different non-normal simulation conditions.

* Simulation of non-normal multivariate distributions using vine copulas
* Implementation of path and multigroup SEM models
* Comparison of ML, MLR, and DWLS estimators
* Analysis of Average Treatment Effects, Average Direct Effects, and Average Indirect Effects
* Evaluation of estimator performance under various degrees of non-normality

### Structure

```
├───study_1
│   ├───main
│   └───margins
├───study_2
│   ├───main
│   └───margins
├───study_path_A
│   ├───main
│   └───margins
└───study_path_B
    ├───main
    └───margins
```
