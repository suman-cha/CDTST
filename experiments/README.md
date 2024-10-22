# Experiments for Cond2ST: Conditional Two-Sample Testing

This folder contains the R code necessary to reproduce the experimental results presented in the paper. The experiments in this directory correspond to the sections and appendices in the paper and are organized into the following categories:

## Structure

### 1. `appendix/`
- This folder contains the code to reproduce the results of **Appendix C.2 to C.5** from the paper.
- Each script corresponds to a specific experiment described in these appendices.

### 2. `real_examples/`
- This folder includes the R scripts used to replicate the real data experiments discussed in **Section 5.3** of the paper.
- The scripts are:
  - **`real_low_dim.R`**: Runs the experiments for low-dimensional real-world data 'diamonds'.
  - **`real_high_dim.R`**: Runs the experiments for high-dimensional real-world data 'superconductivity'.
  
  The necessary data for these experiments can be found in the `data/` folder within this directory.

### 3. `scenarios/`
- This folder contains scripts that replicate the results from **Section 5.2** of the paper. Each script corresponds to a different conditional two-sample test scenario:
  - **Scenario 1**:
    - `scenario1u.R` 
    - `scenario1b.R` 
  - **Scenario 2**:
    - `scenario2u.R` 
    - `scenario2b.R` 
  - **Scenario 3**:
    - `scenario3u.R` 
    - `scenario3b.R` 

These scripts test various assumptions about the underlying distributions of the two populations being compared.

### 4. `CIT_functions.R`
- This script includes the **Conditional Independence Testing** approach, adapted from the [PCM code repository](https://github.com/ARLundborg/pcm_code).
- The method is based on:
  - **Lundborg, A. R., Kim, I., Shah, R. D., and Samworth, R. J. (2022)**. *The Projected Covariance Measure for assumption-lean variable significance testing*. arXiv preprint arXiv:2211.02039.

### 5. `CP_FunFiles.R`
- This script implements a two-sample conditional distribution test using **conformal prediction** and **weighted rank sum** methods, as introduced by:
  - **Hu, X. and Lei, J. (2024)**. *A two-sample conditional distribution test using conformal prediction and weighted rank sum*. Journal of the American Statistical Association, 119(546):1136–1154.

### 6. `all_tests.R`
- This script aggregates all the test methods described in the paper, providing a centralized way to run all of them. It includes tests based on LinearMMD, Classifier-based approach, conformal prediction, conditional independence testing appraoch, and other conditional two-sample testing methods.

### 7. `plot.R`
- This script generates **Figures 1, 2, and 3** from the paper, visualizing the results of the experiments.

### 8. `utils.R`
- Contains utility functions that support the experiment scripts. These include helper functions for data handling, kernel computations, and other necessary utilities for running the experiments.
