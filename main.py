#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
get_ipython().system('pip install lmfit')
from lmfit.models import GaussianModel


# In[2]:


log_df = pd.read_excel("Copy_of_Ecol_Raman.xlsx", sheet_name='ec_log_mean raman')
sta_df = pd.read_excel("Copy_of_Ecol_Raman.xlsx", sheet_name='ec_sta_mean raman')


# In[3]:


# Load each Raman spectrum from its sheet
log_df = pd.read_excel("Copy_of_Ecol_Raman.xlsx",
                       sheet_name="ec_log_mean raman")

sta_df = pd.read_excel("Copy_of_Ecol_Raman.xlsx",
                       sheet_name="ec_sta_mean raman")

# Standardise column names (first two columns are wavenumber + intensity)
log_df = log_df.iloc[:, :2]
sta_df = sta_df.iloc[:, :2]

log_df.columns = ["wavenumber", "intensity"]
sta_df.columns = ["wavenumber", "intensity"]

# Sort just in case
log_df = log_df.sort_values("wavenumber")
sta_df = sta_df.sort_values("wavenumber")

# Plot both on one graph
plt.figure(figsize=(10,6))
plt.plot(log_df["wavenumber"], log_df["intensity"], label="Log phase", linewidth=1.2)
plt.plot(sta_df["wavenumber"], sta_df["intensity"], label="Stationary phase", linewidth=1.2)

# Reverse x-axis (Raman convention)
plt.gca().invert_xaxis()

plt.xlabel("Raman Shift (cm⁻¹)")
plt.ylabel("Intensity (a.u.)")
plt.title("E. coli Raman Spectra: Log vs Stationary Phase")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


# In[4]:


x_log = log_df["wavenumber"].values
y_log = log_df["intensity"].values

# sampling step Δν
dnu_log = np.diff(x_log)
dnu_log_mean = np.mean(dnu_log)

fs_log = 1.0 / dnu_log_mean      # sampling "frequency" in 1/(cm⁻¹)
f_nyq_log = fs_log / 2.0         # Nyquist frequency
N_log = len(x_log)
max_params_log = N_log // 2      # upper bound on independent real parameters

print("LOG PHASE:")
print(f"  N points            = {N_log}")
print(f"  mean Δν             = {dnu_log_mean:.6f} cm⁻¹")
print(f"  Nyquist frequency   = {f_nyq_log:.6f} (1/cm)")
print(f"  Nyquist upper bound ≈ {max_params_log} independent real parameters\n")


# In[5]:


x_sta = sta_df["wavenumber"].values
y_sta = sta_df["intensity"].values

dnu_sta = np.diff(x_sta)
dnu_sta_mean = np.mean(dnu_sta)

fs_sta = 1.0 / dnu_sta_mean
f_nyq_sta = fs_sta / 2.0
N_sta = len(x_sta)
max_params_sta = N_sta // 2

print("STATIONARY PHASE:")
print(f"  N points            = {N_sta}")
print(f"  mean Δν             = {dnu_sta_mean:.6f} cm⁻¹")
print(f"  Nyquist frequency   = {f_nyq_sta:.6f} (1/cm)")
print(f"  Nyquist upper bound ≈ {max_params_sta} independent real parameters")


# In[6]:


# Transforming data into reciprocal space using fourier transform
Y_log = np.fft.rfft(y_log)
freqs_log = np.fft.rfftfreq(N_log, d=dnu_log_mean)
print("LOG PHASE: number of non-negative frequency bins =", len(freqs_log))
print("Max frequency from FFT (should be ~Nyquist):", freqs_log[-1])


# ## The value of 319 parameters using reciprocal space is consistent with upper bound using nyquist shannon theorem
# 
# ## Now we normalise so that the peak fitting is numerically stable

# In[7]:


# Normalise intensities so max = 1
log_df["intensity_norm"] = log_df["intensity"] / log_df["intensity"].max()
sta_df["intensity_norm"] = sta_df["intensity"] / sta_df["intensity"].max()


# In[8]:


# Plotting with normalised data

plt.figure(figsize=(10, 6))
plt.plot(log_df["wavenumber"], log_df["intensity_norm"], label="Log phase (norm)", linewidth=1.2)
plt.plot(sta_df["wavenumber"], sta_df["intensity_norm"], label="Stationary phase (norm)", linewidth=1.2)

plt.gca().invert_xaxis()
plt.xlabel("Raman Shift (cm⁻¹)")
plt.ylabel("Normalised intensity (a.u.)")
plt.title("E. coli Raman spectra (normalised)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


# # The peaks are noticeably more defined now, especially for the log phase data.

# In[9]:


x = log_df["wavenumber"].values
y = log_df["intensity_norm"].values


# In[10]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy import sparse
from scipy.sparse.linalg import spsolve

def baseline_asls(y, lam=1e6, p=0.01, niter=10):
    """Asymmetric least-squares baseline estimate."""
    y = np.asarray(y)
    L = y.size
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L-2))
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.T)
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z

# =============================
# EASY-TO-ADJUST SETTINGS
# =============================
xmin, xmax = 600, 1800          # analysis window

use_smoothing = True
sg_window = 11                  # must be odd; 11/15/21...
sg_poly   = 3

# "Strict" = reject noise bumps; if you want *more* small bumps, lower prom/height
min_sep_cm1  = 1.0              # strict separation (cm^-1)
height_frac  = 0.0001           # strict height floor (fraction of max)
prom_frac    = 0.0005           # strict prominence (fraction of max)

# Baseline fit (shown on plot; used for detection signal only)
use_baseline_for_detection = True
lam = 1e5                        # increase -> straighter baseline; decrease -> more curved
p   = 0.005                      # smaller -> baseline sits more under peaks

# =============================
# 1) Window the spectrum
# =============================
mask = (x >= xmin) & (x <= xmax)
x_use = x[mask]
y_use = y[mask]

# Ensure x increases for robust distance conversion
if x_use[0] > x_use[-1]:
    orderx = np.argsort(x_use)
    x_use = x_use[orderx]
    y_use = y_use[orderx]

# =============================
# 2) Baseline fit (not subtracted in plot)
# =============================
base = baseline_asls(y_use, lam=lam, p=p, niter=10)

# =============================
# 3) Detection signal (baseline-detrended, optional smoothing)
# =============================
y_det = (y_use - base) if use_baseline_for_detection else y_use
y_det = np.clip(y_det, 0, None)

if use_smoothing and len(y_det) >= sg_window:
    w = sg_window if sg_window % 2 == 1 else sg_window + 1
    w = min(w, len(y_det) if len(y_det) % 2 == 1 else len(y_det) - 1)
    y_det = savgol_filter(y_det, window_length=w, polyorder=min(sg_poly, w-2))

# =============================
# 4) Convert min separation (cm^-1) -> points
# =============================
dx = float(np.median(np.diff(x_use)))
distance_pts = int(np.ceil(min_sep_cm1 / dx))

# =============================
# 5) Peak finding (scaled thresholds)
# =============================
ymax = float(np.max(y_det))
height = height_frac * ymax
prominence = prom_frac * ymax

peaks, props = find_peaks(
    y_det,
    height=height,
    prominence=prominence,
    distance=distance_pts
)

peak_pos  = x_use[peaks]
peak_y    = y_use[peaks]               # markers on original spectrum
peak_prom = props["prominences"]

print(f"Detected peaks: {len(peaks)}")
print("Settings used:")
print(f"  window: {xmin}-{xmax} cm^-1")
print(f"  baseline: AsLS (lam={lam:.2g}, p={p})  | used_for_detection={use_baseline_for_detection}")
print(f"  smoothing: {use_smoothing} (sg_window={sg_window}, sg_poly={sg_poly})")
print(f"  min_sep_cm1={min_sep_cm1} -> distance_pts={distance_pts} (dx≈{dx:.3f} cm^-1)")
print(f"  height_frac={height_frac} -> height={height:.3g}")
print(f"  prom_frac={prom_frac} -> prominence={prominence:.3g}")

# =============================
# 6) Plot (spectrum + baseline + peaks)
# =============================
plt.figure(figsize=(10, 4))
plt.plot(x_use, y_use, linewidth=1, label="Log phase (norm)")
plt.plot(x_use, base, "--", linewidth=1.2, label="Baseline (AsLS)")
plt.plot(peak_pos, peak_y, "ro", markersize=4, label="Detected peaks")
plt.gca().invert_xaxis()
plt.xlabel("Raman shift (cm$^{-1}$)")
plt.ylabel("Normalised intensity (a.u.)")
plt.title("Peak finding result (600–1800 cm$^{-1}$)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# =============================
# 7) Sort peaks by prominence (for initial centres later)
# =============================
order = np.argsort(peak_prom)[::-1]
peak_pos_sorted  = peak_pos[order]
peak_prom_sorted = peak_prom[order]
peak_y_sorted    = peak_y[order]

top_n = min(20, len(peak_pos_sorted))
print(f"\nTop {top_n} peaks by prominence:")
for i in range(top_n):
    print(f"{i+1:>2}. {peak_pos_sorted[i]:8.2f} cm^-1   prominence={peak_prom_sorted[i]:.4f}")


# # Fitting and BFI

# In[11]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import voigt_profile
from scipy.optimize import lsq_linear

# ===========================
# SETTINGS
# ===========================
sigma_fixed = 6.0
gamma_fixed = 6.0
baseline_deg = 2

min_sep_fit_cm1 = 15.0

Amax = 2.0
bmax = 2.0

max_steps = 20          # force plotting up to this many peaks
force_full_curve = True # <-- key: keep going even if ΔlnBFI < 0

# ===========================
x_use = np.asarray(x_use, float)
y_use = np.asarray(y_use, float)
cands = np.asarray(peak_pos_sorted, float)

if x_use[0] > x_use[-1]:
    o = np.argsort(x_use)
    x_use, y_use = x_use[o], y_use[o]

n = len(x_use)
n_base = baseline_deg + 1

def x_scaled(x):
    return 2*(x - x.min())/(x.max() - x.min()) - 1.0

def make_voigt_col(x, center, sigma, gamma):
    v = voigt_profile(x - center, sigma, gamma)
    vmax = np.max(v)
    return v / vmax if (vmax > 0 and np.isfinite(vmax)) else np.zeros_like(x)

def build_X(x, centers, deg, sigma, gamma):
    xs = x_scaled(x)
    B = np.vstack([xs**i for i in range(deg+1)]).T
    V = np.vstack([make_voigt_col(x, c, sigma, gamma) for c in centers]).T if len(centers) else np.zeros((len(x),0))
    return np.hstack([B, V])

def fit_nonneg_amps(X, y, n_base):
    k = X.shape[1]
    lb = np.full(k, -np.inf)
    ub = np.full(k,  np.inf)
    lb[n_base:] = 0.0  # amplitudes >= 0

    res = lsq_linear(X, y, bounds=(lb, ub), lsmr_tol='auto', verbose=0)
    beta = res.x
    yhat = X @ beta
    r = y - yhat
    RSS = float(np.sum(r**2))

    dof = max(n - k, 1)
    s2 = RSS / dof

    XtX = X.T @ X
    try:
        XtX_inv = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        XtX_inv = np.linalg.pinv(XtX)

    cov = s2 * XtX_inv
    condX = float(np.linalg.cond(X))
    rank = np.linalg.matrix_rank(X)
    return beta, yhat, RSS, s2, cov, rank, condX

def ln_bfi(s2, cov, prior_widths, n):
    sign, logdet = np.linalg.slogdet(cov)
    if sign <= 0 or not np.isfinite(logdet) or s2 <= 0 or not np.isfinite(s2):
        return np.nan
    return float(
        -0.5*n*np.log(2*np.pi)
        -0.5*n*np.log(s2)
        -0.5*n
        +0.5*logdet
        -np.sum(np.log(prior_widths))
    )

def can_add(existing, cand, min_sep):
    return all(abs(cand - e) >= min_sep for e in existing)

# ---------- baseline-only start ----------
selected = []
history = []

X0 = build_X(x_use, [], baseline_deg, sigma_fixed, gamma_fixed)
beta0, yhat0, RSS0, s20, cov0, rank0, cond0 = fit_nonneg_amps(X0, y_use, n_base)
prior0 = np.array([2*bmax]*n_base, float)
lnB0 = ln_bfi(s20, cov0, prior0, n=n)

prev_lnB = lnB0

print("Forward selection (forced full curve)" if force_full_curve else "Forward selection (stops at saturation)")
print(f"n={n}, baseline_deg={baseline_deg}, sigma,gamma=({sigma_fixed},{gamma_fixed}), min_sep={min_sep_fit_cm1}")
print(f"Priors: Δb={2*bmax}, ΔA={Amax}, amps>=0")
print("")
print("step  center(cm^-1)   k  rank  cond(X)        RSS         s^2        lnBFI       ΔlnBFI")
print("----------------------------------------------------------------------------------------")
print(f"{0:>4d}  {'(baseline)':>12s}  {X0.shape[1]:>2d}  {rank0:>4d}  {cond0:>9.2e}  {RSS0:>10.3e}  {s20:>9.3e}  {lnB0:>10.4f}  {'-':>9s}")

for step in range(1, min(max_steps, len(cands)) + 1):
    best = None

    for c in cands:
        if c in selected:
            continue
        if not can_add(selected, float(c), min_sep_fit_cm1):
            continue

        centers_try = selected + [float(c)]
        X = build_X(x_use, centers_try, baseline_deg, sigma_fixed, gamma_fixed)
        beta, yhat, RSS, s2, cov, rank, condX = fit_nonneg_amps(X, y_use, n_base)

        prior = np.array([2*bmax]*n_base + [Amax]*len(centers_try), float)
        lnB = ln_bfi(s2, cov, prior, n=n)
        if not np.isfinite(lnB):
            continue

        dlnB = lnB - prev_lnB
        if (best is None) or (dlnB > best[1]):
            best = (lnB, dlnB, float(c), X.shape[1], rank, condX, RSS, s2, yhat)

    if best is None:
        print("\nNo remaining candidates satisfy min_sep. Stopping.")
        break

    lnB, dlnB, c_best, k_best, rank_best, cond_best, RSS_best, s2_best, yhat_best = best
    selected.append(c_best)
    history.append((step, c_best, lnB, dlnB, k_best, rank_best, cond_best, RSS_best, s2_best, yhat_best))

    print(f"{step:>4d}  {c_best:>12.2f}  {k_best:>2d}  {rank_best:>4d}  {cond_best:>9.2e}  {RSS_best:>10.3e}  {s2_best:>9.3e}  {lnB:>10.4f}  {dlnB:>9.4f}")

    prev_lnB = lnB

# ---- arrays for plotting
steps = np.array([h[0] for h in history], int)
lnBFI_arr = np.array([h[2] for h in history], float)
dln_arr = np.array([h[3] for h in history], float)

# argmax peak (this is N*)
best_idx = int(np.nanargmax(lnBFI_arr))
N_star = int(steps[best_idx])
ln_star = float(lnBFI_arr[best_idx])

print("\nArgmax along curve:")
print(f"  N* = {N_star} at lnBFI = {ln_star:.4f}")

# plots
plt.figure(figsize=(8,4))
plt.plot(steps, lnBFI_arr, marker="o")
plt.axvline(N_star, linestyle="--")
plt.scatter([N_star], [ln_star], zorder=5)
plt.xlabel("Step (number of selected peaks)")
plt.ylabel("ln(BFI)")
plt.title(f"ln(BFI) up to step {steps.max()} (peak at N*={N_star})")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8,4))
plt.plot(steps, dln_arr, marker="o")
plt.axhline(0.0, linestyle="--")
plt.xlabel("Step")
plt.ylabel("Δ ln(BFI)")
plt.title("Increment per added peak (can go negative after saturation)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()


# In[ ]:




