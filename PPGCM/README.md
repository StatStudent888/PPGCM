# PPGCM Reproducibility Materials

This file accompanies the reproducibility materials for the paper:
**Prediction-Powered Conditional Independence Testing: A Generalized Covariance Measure Approach**


# Description of the paper

This paper proposes '''PPGCM', a conditional-independence testing procedure for
\[
X \perp\!\!\!\perp Y \mid Z,
\]
that incorporates an auxiliary variable `S` and, when available, unlabeled observations containing `(X, S, Z)` but not `Y`.

The implementation in this repository reports three GCM-based procedures:

1. **PPGCM (opt)**, which estimates the variance-minimizing prediction-powered weight;
2. **PPGCM (`w = 1`)**, which uses a fixed prediction-powered weight of one;
3. **GCM**, the classical cross-fitted generalized covariance measure (GCM) based only on the labeled observations.

The reproducibility code also contains simulation comparisons with **CRT**, **Maxway CRT**, and **SGCM**, oracle experiments studying local power and asymptotic variance under correctly specified and misspecified prediction models, and real-data analyses using a **GEO breast-cancer gene-expression dataset** and **NHANES**.


# Repository structure

The repository is organized around one reusable implementation file together with simulation and real-data reproduction scripts.

## `PPGCM.R`: main method implementation

- `PPGCM/PPGCM.R`: main implementation of PPGCM and the cross-fitted GCM.
  - The main user-facing function is `ppicit_test()`.
  - `label_data` should contain `Y`, `X`, `S`, and `Z` variables.
  - `unlabel_data` is optional and should contain `X`, `S`, and `Z`, but not `Y`.
  - The nuisance/prediction regressions can be fitted by lasso (`glmnet`) or generalized random forests (`grf`).
  - The returned object contains results for PPGCM (opt), PPGCM (`w = 1`), and GCM, including test statistics, estimated asymptotic variances, standard errors, z-scores, p-values, and rejection decisions.

Users who want to apply PPGCM to their own labeled and unlabeled data should start from `PPGCM/PPGCM.R` and the calls to `ppicit_test()` in the simulation or real-data scripts.

## `simulation/`: simulation studies

### `simulation/hardness/`: local-power and variance experiments

- `PPGCM/simulation/hardness/main.R`: main driver for the fixed-`t` local-alternative experiment with `c_n = t / sqrt(n)`. It runs the simulation across labeled sample sizes and produces empirical-power and asymptotic-variance plots.
- `PPGCM/simulation/hardness/ppgcm_hardness_correctly_specified.R`: oracle experiment with a correctly specified prediction model.
- `PPGCM/simulation/hardness/ppgcm_hardness_misspecified.R`: oracle experiment with a misspecified prediction model.
- `PPGCM/simulation/hardness/ppigcm_hardness_plot_functions.R`: plotting functions used by `simulation/hardness/main.R`.

### `simulation/type-I_power/`: type-I error and power simulations

- `PPGCM/simulation/type-I_power/generate_data.R`: generates the labeled and unlabeled samples for the linear and nonlinear simulation designs.
- `PPGCM/simulation/type-I_power/main.R`: runs the main PPGCM/GCM Monte Carlo experiments in Section 4.2.1.
  - Supports both nonlinear and linear designs.
  - The default active block runs the nonlinear design with generalized random forests.
  - A commented block gives the corresponding linear/lasso setting.
- `PPGCM/simulation/type-I_power/CRT_MaxwayCRT/maxway.R`: helper implementation for CRT and Maxway CRT.
- `PPGCM/simulation/type-I_power/CRT_MaxwayCRT/main.R`: runs CRT and Maxway CRT under the same linear/nonlinear simulation designs used for the PPGCM study.
- `PPGCM/simulation/type-I_power/sgcm/sgcm_gcm_ppgcm_comparison.R`: compares oracle GCM, cross-fitted GCM, oracle SGCM, and PPGCM under a null distribution satisfying `X ⟂ Y | Z` for which directly adjusting for `S` can invalidate SGCM.

## `real data/`: real-data analyses

### `real data/GEO_breast_cancer/`

- `PPGCM/real data/GEO_breast_cancer/GEO_real_data.R`: reproduces the GEO breast-cancer gene-expression analysis.
- `GEO_data.csv`: processed GEO breast-cancer gene-expression data used by the analysis script. Because of its large file size, this dataset is provided separately as a GitHub Release asset.

The GEO analysis uses:

- `Y`: expression of `CCND1`;
- `X`: expression of `ESR1`;
- `S`: expression of `PPFIA1`;
- `Z`: the remaining PAM50 genes after excluding `ESR1`.

The script restricts the data to luminal breast-cancer subtypes, removes anomalous observations using an isolation forest, standardizes the retained expression measurements, evaluates the predictive contribution of `S`, applies full-data GCM diagnostics, and then constructs an artificial semi-supervised experiment comparing GCM, PPGCM (`w = 1`), PPGCM (opt), CRT, and Maxway CRT.

### `real data/NHANES/`

- `PPGCM/real data/NHANES/label.csv`: processed labeled NHANES sample containing the outcome `Y` together with `X`, `S`, and `Z`.
- `PPGCM/real data/NHANES/unlabel.csv`: processed unlabeled NHANES sample containing `X`, `S`, and `Z` but not `Y`.
- `PPGCM/real data/NHANES/main.R`: reproduces the NHANES semi-supervised analysis.

The NHANES analysis uses:

- `Y`: two-hour OGTT plasma glucose;
- `X`: waist circumference;
- `S`: HbA1c;
- `Z`: BMI, age, systolic blood pressure, HDL cholesterol, total cholesterol, and creatinine.

The supplied processed files contain 2,239 labeled observations and 2,827 unlabeled observations. The script standardizes the variables, checks `X ⟂ S | Z`, assesses the predictive contribution of `S`, applies GCM and the two PPGCM variants, and then compares them with CRT and Maxway CRT.

# Simulation and real-data reproduction map

| Output / analysis | Purpose | Scripts to run | Order / notes |
|---|---|---|---|
| Main type-I error and power results in Section 4.2.1 (Table 1, Tables S.3--S.9) | Compare cross-fitted GCM, PPGCM (`w = 1`), and PPGCM (opt) under linear and nonlinear designs | `PPGCM/simulation/type-I_power/generate_data.R`, `PPGCM/simulation/type-I_power/main.R` | `main.R` sources `generate_data.R` automatically. Edit the active parameter block at the end of `main.R` to reproduce different `(n, N, c)` settings. |
| CRT and Maxway CRT results for the main simulation designs (Table 1, Tables S.3--S.5) | Compare CRT and Maxway CRT with the GCM/PPGCM procedures under the same data-generating mechanisms | `PPGCM/simulation/type-I_power/generate_data.R`, `PPGCM/simulation/type-I_power/CRT_MaxwayCRT/maxway.R`, `PPGCM/simulation/type-I_power/CRT_MaxwayCRT/main.R` | `main.R` loads the helper files. Linear and nonlinear nuisance models are handled separately inside the script. |
| Table S.1 | Illustrate the behavior of SGCM when `S` is conditionally associated with both `X` and `Y` although `X ⟂ Y \| Z` holds | `PPGCM/simulation/type-I_power/sgcm/sgcm_gcm_ppgcm_comparison.R` | The script runs the full simulation and prints `report_table`. |
| Figure 1 | Compare empirical power and estimated asymptotic variance of GCM, PPGCM (`w = 1`), and PPGCM (opt) under `c_n = t/sqrt(n)` | `PPGCM/simulation/hardness/ppgcm_hardness_correctly_specified.R`, `PPGCM/simulation/hardness/ppigcm_hardness_plot_functions.R`, `PPGCM/simulation/hardness/main.R` | This is the default active setting in `hardness/main.R`. The main outputs are `res_by_n`, `p1`, and `p2`. |
| Figure S.2 | Repeat the fixed-`t` experiment when the prediction function is deliberately misspecified | `PPGCM/simulation/hardness/ppgcm_hardness_misspecified.R`, `PPGCM/simulation/hardness/ppigcm_hardness_plot_functions.R`, `PPGCM/simulation/hardness/main.R` | In `hardness/main.R`, comment out the correctly specified source/run block and uncomment the misspecified block and its plotting ranges. |
| Table 2 | Gene-expression analysis using `CCND1`, `ESR1`, `PPFIA1`, and PAM50 genes; includes artificial labeled/unlabeled splits | `PPGCM/real data/GEO_breast_cancer/GEO_real_data.R`, with `GEO_data.csv` | The script reads `GEO_data.csv` from `PPGCM/real data/GEO_breast_cancer/`. Run from the parent directory containing `PPGCM/`. |
| Table S.2 | Genuine semi-supervised analysis using OGTT glucose, waist circumference, HbA1c, and cardiometabolic covariates | `PPGCM/real data/NHANES/main.R`, with `label.csv` and `unlabel.csv` | The processed NHANES data are included. Run from the parent directory containing `PPGCM/`. |

# Requirements

All supplied analysis code is written in R. The packages below are those explicitly required or loaded by the code.

## Main software

- R version: 4.3.1
- Operating system: Windows with Rtools

## Core PPGCM implementation

| Package | Version used |
|---|---|
| `glmnet` | 4.1-8 |
| `grf` | 2.3.2 |

## Simulation and plotting code

| Package | Version used |
|---|---|
| `glmnet` | 4.1-8 |
| `grf` | 2.3.2 |
| `randomForest` | 4.7-1.2 |
| `MASS` | 7.3-60 |
| `CompQuadForm` | 1.4.3 |
| `ggplot2` | 4.0.3 |
| `future` | 1.33.1 |
| `future.apply` | 1.11.1 |
| `RhpcBLASctl` | 0.23-42 |

## Real-data code

| Package | Version used |
|---|---|
| `dplyr` | 1.1.4 |
| `isotree` | 0.6.1-5 |
| `randomForest` | 4.7-1.2 |
| `glmnet` | 4.1-8 |
| `limma` | 3.56.2 |

# Data access

## NHANES

The processed NHANES data required by `PPGCM/real data/NHANES/main.R` are included in the repository:

- `PPGCM/real data/NHANES/label.csv`;
- `PPGCM/real data/NHANES/unlabel.csv`.

The labeled sample contains the two-hour OGTT outcome, whereas the unlabeled sample does not. Both contain waist circumference, HbA1c, and the conditioning covariates used by the analysis.

## GEO breast-cancer data

The GEO breast-cancer analysis uses `GEO_data.csv`, which can be downloaded from the
[latest GitHub Release](https://github.com/StatStudent888/PPGCM/releases/latest/download/GEO_data.zip).

After downloading the file, place it at:

- `PPGCM/real data/GEO_breast_cancer/GEO_data.csv`.

The data file contains the processed gene-expression measurements used by `GEO_real_data.R`, including a subtype column and gene-expression columns containing `CCND1`, `ESR1`, `PPFIA1`, and the `PAM50` genes used as conditioning variables. The script restricts the analysis to luminal breast-cancer subtypes, performs anomaly detection and standardization, and constructs the labeled/unlabeled splits used for the semi-supervised comparison.



