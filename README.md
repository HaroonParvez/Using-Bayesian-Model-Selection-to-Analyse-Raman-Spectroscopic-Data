# Bayesian Model Selection for Raman Spectroscopy of *E. coli*

This repository contains the code and data for my third-year final project, where I applied data analysis and Bayesian statistical methods to Raman spectroscopic data from the bacteria E.coli.

The project focuses on analysing spectral data to identify the biomolecular components that best explain each measurement using the Bayes Factor Integral (BFI). By selecting only statistically supported components, the approach improves model reliability while reducing overfitting and demonstrating a practical data-driven decision-making workflow.

## Project Overview

Raman spectroscopy provides a non-invasive and label-free way to study biological samples by measuring molecular vibrational information. However, Raman spectra from bacterial cells are difficult to interpret because the measured signal is usually a mixture of many overlapping biomolecular contributions.

In this project, measured *E. coli* Raman spectra are modelled as a linear combination of reference biomolecular spectra, with additional baseline terms included to account for background effects. Different combinations of reference spectra are tested, and model selection metrics are used to decide which model is most appropriate.

The main focus is the Bayes Factor Integral, which combines information about:

- quality of fit
- model complexity
- parameter uncertainty
- prior parameter volume
- correlations between fitted parameters

This allows the project to compare models more carefully than using fit quality alone.

## Repository Contents

```text
.
├── Copy_of_Ecol_Raman.xlsx
├── HaroonParvezFinalProject.ipynb
├── README.md
└── .gitattributes
```

## Files

### `HaroonParvezFinalProject.ipynb`

This is the main Jupyter Notebook for the project. It contains the full analysis workflow, including:

- loading the Raman spectroscopy data
- preparing the measured and reference spectra
- fitting linear combination models
- applying non-negative least squares fitting
- calculating model comparison statistics
- implementing the Bayes Factor Integral
- comparing BFI with other metrics such as AIC, BIC, chi-squared, reduced chi-squared and R²
- producing plots and results used in the final project

### `Copy_of_Ecol_Raman.xlsx`

This spreadsheet contains the Raman spectroscopy data used in the project, including the measured *E. coli* spectra and reference biomolecular spectra.

## Method Summary

The measured Raman spectrum is treated as a combination of baseline terms and selected reference spectra. In matrix form, this can be written as:

```text
y = Aβ + ε
```

where:

- `y` is the measured Raman spectrum
- `A` is the design matrix containing baseline terms and selected reference spectra
- `β` is the vector of fitted coefficients
- `ε` is the residual noise

The project uses non-negative fitting because biomolecular contributions should not have negative concentrations. Models are then compared by adding reference spectra step by step and evaluating whether each added component is statistically justified.

## Model Selection

Several model selection and fit-quality metrics are compared:

- Bayes Factor Integral (BFI)
- Akaike Information Criterion (AIC)
- Bayesian Information Criterion (BIC)
- chi-squared
- reduced chi-squared
- coefficient of determination (R²)

The BFI is the main metric because it includes an Occam factor, meaning that unnecessary extra parameters are penalised. This is important because adding more spectra will usually improve the raw fit, but that does not mean every added component is genuinely supported by the data.

## Main Aim

The aim of the project is to improve the interpretation of Raman spectra from *E. coli* by selecting the most justified biomolecular components in the spectrum. This helps avoid overfitting and gives a more reliable interpretation of the biological information contained in the Raman data.

## Requirements

The notebook uses Python and common scientific computing libraries, including:

- NumPy
- pandas
- SciPy
- Matplotlib
- scikit-learn
- openpyxl

These can be installed using:

```bash
pip install numpy pandas scipy matplotlib scikit-learn openpyxl
```

## How to Run

1. Clone or download this repository.
2. Open `HaroonParvezFinalProject.ipynb` in Jupyter Notebook, JupyterLab, VS Code, or Google Colab.
3. Make sure `Copy_of_Ecol_Raman.xlsx` is in the same folder as the notebook.
4. Run the notebook cells in order.

## Author

Haroon Parvez  
Third Year Project  
Queen Mary University of London  
Physics and Data Science
