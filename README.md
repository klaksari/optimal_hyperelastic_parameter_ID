# Optimal Hyperelastic Parameter Identification  
MATLAB toolbox for stress-material Jacobian analysis, experiment design,  
and cross-validation of classical incompressible hyperelastic models  
(Neo-Hookean, Mooney-Rivlin, Ogden).



## 📦 What’s in here?

| File | Purpose | Key outputs |
|------|---------|-------------|
| **`MooneyRivlin.m`** | √det / √κ maps for Mooney-Rivlin | Condition/Determinant contour plots |
| **`NeoHookean.m`** | Same for Neo-Hookean | Single-parameter sanity check |
| **`Ogden.m`** | Parametric sweep over α-pairs (2-term Ogden) | Heat-maps + critical-line overlays |
| **`OptimalvsRandomvsLinear.m`** | Monte-Carlo comparison (optimal vs random vs linear λ choices) | KDE of fitted *C1*, *C2* |
| **`Shaded_Figure.m`** | Cross-prediction envelopes (UT → ET/PS) | 3×3 scatter + shaded-band dashboard |
| **`Percentage_Table_MR.m`** | **Brute-force grid audit** for all four loading modes, optional DEIM acceleration, CSV summary + 2-panel dashboard | • tables (`*_results.csv`) <br>• determinant/condition curves |
| **`deim()`** *(sub-function)* | Discrete Empirical Interpolation used when `n > 2` | Row selection indices |

**Why so many scripts?**  
> The repo mirrors the pipeline described in our paper/review response – from exhaustive
> enumeration of λ-space to synthetic-data cross-validation.

---

## 🚀 Quick start

% clone & cd
git clone https://github.com/klaksari/optimal_hyperelastic_parameter_ID.git
cd optimal_hyperelastic_parameter_ID

% 1. generate Mooney-Rivlin λ-maps
MooneyRivlin              % creates ./figures/MR_condition_* etc.

% 2. run a 10 000-iteration Monte-Carlo comparison (≈ 2 min)
OptimalvsRandomvsLinear

% 3. exhaustive grid audit (⚠️ large – adjust num_points first!)
Percentage_Table_MR


Interpreting the metrics
√|det(JᵀJ)| → information volume (D-optimality surrogate).
values < 1 ⇒ rank-deficient neighbourhood

√κ(JᵀJ) → noise amplification.
values > 20 ⇒ lose > 1 significant digit

The green lines (Ogden / Mooney-Rivlin scripts) mark the recommended
operating region: √κ < 20 and √det ≥ 1.

If this code helped your research, please cite:


title   = {Stress–Material Jacobian Framework for Repeatable Hyperelastic Characterisation},
author  = {Asadi, A. and Laksari, K.},
year    = {2025},
note    = {GitHub repository \url{https://github.com/klaksari/optimal_hyperelastic_parameter_ID}}


Happy Stretching!
