# third-year-project-raman-bfi
# Raman Spectroscopy – Bayesian Factor of Information (BFI) Analysis

This repository contains a snapshot of a third-year undergraduate project analysing the
**information content of Raman spectroscopy data** using a **Bayesian evidence / BFI approach**.

The purpose of this snapshot is **code review and validation** (not deployment or reuse).

---

## Project aim

Raman spectra often contain many visible peaks, but not all of them represent
**independent, statistically supported information**.

The core question of this project is:

> *How many independent spectral features are supported by the data, once model
complexity is properly penalised?*

This is addressed using a **Bayesian Factor of Information (BFI)** framework, rather than
heuristic peak counting.

---

## Data overview

The input CSV contains three conceptually different data components:

1. **Log spectrum**
   - A single representative Raman spectrum.
   - Used to estimate the *maximum* apparent spectral information content.

2. **STA spectrum**
   - A statistically processed version of the Raman data.
   - Used to test whether the information limit is noise-driven or structural.

3. **Biological data**
   - Non-spectral variables related to cells / biological conditions.
   - These are *not* fitted with peaks and are not treated as spectra.
   - They provide biological context rather than spectral information content.

This snapshot focuses on the **log spectrum** (and later STA) for BFI analysis.

---

## Analysis pipeline (current approach)

The analysis is implemented in a single Jupyter notebook (`main.ipynb`) and proceeds as follows:

1. **Baseline estimation**
   - Asymmetric least-squares (AsLS) baseline correction.
   - Baseline is used for peak detection but not over-subtracted in fitting.

2. **Peak detection**
   - SciPy `find_peaks` with physically motivated thresholds:
     - minimum separation (cm⁻¹)
     - prominence relative to signal scale
   - Produces a list of candidate peak centres.

3. **Spectral model**
   - Polynomial baseline (degree 2) + fixed-shape Voigt profiles.
   - Voigt widths (σ, γ) are fixed to physically reasonable values.
   - Only **amplitudes** are fitted (linear least squares).

4. **Bayesian evidence / ln(BFI)**
   - Evidence computed using a Laplace approximation.
   - Includes:
     - residual variance estimate
     - covariance determinant
     - explicit parameter priors
     - Occam penalty for model complexity

5. **Forward model selection**
   - Peaks are added sequentially.
   - At each step, ln(BFI) is recomputed.
   - The optimal number of peaks N* is identified as the **maximum of ln(BFI)**.
   - Δln(BFI) is also inspected to identify the onset of overfitting.

---

## Key issue encountered (and resolved)

### Problem
Initial implementations produced **monotonically increasing ln(BFI)** values,
sometimes reaching extremely large magnitudes.  
This indicated that model complexity was not being penalised correctly.

Root causes included:
- inconsistent noise variance definitions
- missing or mis-scaled prior volumes
- overly permissive peak selection (effectively double-counting structure)
- hidden state from notebook execution order

### Resolution
The current working solution fixes this by:
- using a consistent residual variance definition tied to degrees of freedom
- explicitly including prior widths in ln(BFI)
- fixing Voigt shapes to reduce non-identifiability
- enforcing a physically motivated minimum peak separation
- running forward selection with explicit stopping logic

As a result:
- ln(BFI) now **rises, peaks, and then decreases**
- Δln(BFI) becomes negative beyond the optimum
- the model exhibits a clear Occam turnover

This behaviour is the expected Bayesian outcome.

---

## Current results (log spectrum)

With physically reasonable settings (e.g. minimum separation ≈ 15 cm⁻¹):

- ln(BFI) peaks at **N* ≈ 12 peaks**
- Adding further peaks decreases evidence
- The result is stable to reasonable changes in constraints
  (e.g. stricter separation gives N* ≈ 11)

This demonstrates that the Raman spectrum contains a **finite, quantifiable number
of independent spectral degrees of freedom**, despite many visually identifiable peaks.

---

## Files in this snapshot

- `main.ipynb`
  - Source of truth.
  - Full analysis with plots and diagnostics.
  - Intended to be run top-to-bottom.

- `main.py`
  - Auto-exported linear script (`nbconvert`).
  - Provided for code review and static analysis.
  - Not edited directly.

- `README.md`
  - Project context for reviewers.

---

## Notes for code reviewers (Codex)

- The notebook is stateful; `main.py` represents the linearised execution.
- The primary focus of review should be:
  - correctness of the ln(BFI) calculation
  - consistency of noise variance handling
  - model selection logic and stopping criteria
  - separation of spectral vs biological data usage
- This repository is **not** intended as a reusable software package.

---

## Status

The analysis pipeline is currently **stable and working as intended**.
Further work applies the same framework to the STA spectrum and compares
spectral information content to biological variability.


