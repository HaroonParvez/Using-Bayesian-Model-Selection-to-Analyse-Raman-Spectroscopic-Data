# Bayesian Model Selection for Raman Spectroscopy of *E. coli*

This repository contains the code and data for my third-year final project, where I applied data analysis and Bayesian statistical methods to Raman spectroscopic data from the bacteria E.coli.

The project focuses on analysing spectral data to identify the biomolecular components that best explain each measurement using the Bayes Factor Integral (BFI). By selecting only statistically supported components, the approach improves model reliability while reducing overfitting and demonstrating a practical data-driven decision-making workflow.

**📄 Read the full report:** [FinalProjectReport_HaroonParvez_220475901.pdf](FinalProjectReport_HaroonParvez_220475901.pdf) — recommended starting point, includes full theory, methodology, results and discussion.

**💻 Analysis notebook:** [HaroonParvezFinalProject.ipynb](HaroonParvezFinalProject.ipynb) — full code implementation.

---

## Why this project

Raman spectra from bacterial cells contain overlapping signals from dozens of biomolecules, making it hard to know how many are actually detectable versus how many would just be fitting noise. This project adapts a Bayesian model selection method — originally developed for X-ray spectroscopy (EXAFS) — to a completely different domain (biological Raman spectroscopy), and benchmarks it against standard statistical criteria (AIC, BIC, χ², R²) to demonstrate a measurable advantage in handling correlated model parameters.

This is fundamentally a **model selection and statistical inference problem**: given a signal that can be explained by many overlapping candidate variables, how do you objectively decide how many to include without overfitting?

## Project Overview

*E. coli* Raman spectra (log and stationary growth phase) are modelled as linear combinations of 15 reference biomolecular spectra, plus a polynomial baseline to absorb residual background. A forward selection algorithm builds up the model one component at a time, using the Bayes Factor Integral to decide which additions are statistically justified — accounting for:

- goodness of fit
- model complexity
- parameter uncertainty (via the covariance matrix)
- correlations between fitted parameters
- prior parameter ranges

This is compared directly against conventional model selection criteria (AIC, BIC, χ², reduced χ², R²) to show where accounting for parameter correlation changes the outcome.

## Report

The full report (31 pages) covers:

- **Theory** — Raman spectroscopy fundamentals, linear combination fitting, the overfitting problem, and the Bayesian derivation of the BFI
- **Methodology** — data preprocessing, design matrix construction, bounded non-negative least squares fitting, and the forward selection procedure
- **Results** — optimal model complexity for each growth phase (8 components for log phase, 10 for stationary), fit quality assessment, and a full comparison against conventional criteria
- **Discussion** — biological interpretation of the selected components, limitations of the approach, and comparison with existing literature

Reading the report first gives the full context before diving into the code.

## Repository Contents

```text
.
├── Copy_of_Ecol_Raman.xlsx
├── FinalProjectReport_HaroonParvez_220475901.pdf
├── HaroonParvezFinalProject.ipynb
├── README.md
└── .gitattributes
```

## Files

### `HaroonParvezFinalProject.ipynb`

The full analysis pipeline, including:

- Data loading and preprocessing (normalisation, alignment of measured and reference spectra)
- Design matrix construction (baseline + reference spectra)
- Bounded non-negative least squares fitting (`scipy.optimize.lsq_linear`)
- A from-scratch implementation of the Bayes Factor Integral, including the Occam factor and covariance-based parameter correlation penalty
- Forward selection algorithm for model comparison
- Independent benchmarking against AIC, BIC, χ², reduced χ², and R²
- All plots and results referenced in the report

### `Copy_of_Ecol_Raman.xlsx`

Raw Raman spectroscopy data: measured *E. coli* spectra (log and stationary growth phase) and the 15 reference biomolecular spectra used for fitting.

## Method Summary

The measured spectrum is modelled as:

```text
y = Aβ + ε
```

where `A` is the design matrix (baseline terms + selected reference spectra), `β` is the vector of fitted coefficients, and `ε` is residual noise. Fitting is constrained to non-negative biomolecular coefficients, since a biomolecule can only contribute signal, never subtract it.

## Model Selection

Six criteria are implemented and compared, each running its own independent forward selection:

- **Bayes Factor Integral (BFI)** — the primary metric, incorporating an Occam factor that penalises correlated or poorly-constrained parameters
- Akaike Information Criterion (AIC)
- Bayesian Information Criterion (BIC)
- χ² and reduced χ²
- Coefficient of determination (R²)

The BFI consistently selects the smallest fully-supported model, while conventional criteria tend to include additional components that improve fit only marginally — a direct illustration of why accounting for parameter correlation matters in overlapped spectral data.

## Skills Demonstrated

- Statistical model selection and Bayesian inference
- Constrained optimisation (bounded non-negative least squares)
- Custom implementation of a published statistical method (BFI) in a novel domain
- Benchmarking against standard model selection criteria (AIC, BIC, χ²)
- Data preprocessing and signal alignment
- Scientific data visualisation (Matplotlib)
- Working with real experimental biological data

## Requirements

```bash
pip install numpy pandas scipy matplotlib openpyxl
```

## How to Run

1. Clone or download this repository.
2. Open `HaroonParvezFinalProject.ipynb` in Jupyter Notebook, JupyterLab, VS Code, or Google Colab.
3. Ensure `Copy_of_Ecol_Raman.xlsx` is in the same folder as the notebook.
4. Run the cells in order.

## Author

Haroon Parvez
BSc Physics, Queen Mary University of London
