# HASOD: Hybrid Adaptive Screening-Optimization Design

A three-phase sequential experimental design framework that unifies factor screening and process optimization into a single adaptive workflow. HASOD achieves **97.08% factor detection accuracy** -- a 13.75 percentage point improvement over the best competing method -- while adaptively allocating experimental runs based on observed factor structure.

---

## Table of Contents

- [Problem Statement](#problem-statement)
- [Methodology](#methodology)
  - [Phase 1: Intelligent Screening via Modified DSD](#phase-1-intelligent-screening-via-modified-dsd)
  - [Phase 2: Adaptive Selective Augmentation](#phase-2-adaptive-selective-augmentation)
  - [Phase 3: GP-Based Optimization with Uncertainty-Guided Refinement](#phase-3-gp-based-optimization-with-uncertainty-guided-refinement)
- [Novelty and Contributions](#novelty-and-contributions)
- [Benchmark Results](#benchmark-results)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [References](#references)

---

## Problem Statement

In experimental science and engineering, practitioners face a common two-stage challenge:

1. **Factor screening** -- identifying which of many controllable factors actually influence the response.
2. **Process optimization** -- finding the factor settings that maximize (or minimize) the response.

Existing approaches force a choice between these goals. Screening designs (Plackett-Burman, Definitive Screening Designs) identify important factors but do not optimize. Optimization methods (Latin Hypercube Sampling, Bayesian Optimization, Box-Behnken) fit response surfaces but provide no factor classification. The traditional workaround -- running a screening study first, then a separate optimization study -- uses a fixed sequence that does not adapt to the structure revealed by the data. This leads to wasted experimental runs when the system is sparse, and missed factors when the system is dense or interaction-heavy.

HASOD fills this gap by integrating adaptive factor screening with GP-based optimization in a single principled sequential framework.

---

## Methodology

HASOD operates in three sequential phases. Each phase uses the information gained from previous phases to decide what to do next.

### Phase 1: Intelligent Screening via Modified DSD

**Goal:** Identify which factors are critical, moderate, or negligible, and detect significant two-factor interactions.

**Design Construction.** For *k* factors, a Modified Definitive Screening Design (M-DSD) is generated with *n_1 = 2k + 3* runs. The design includes a center point, *k* pairs of fold-over runs (where each pair isolates one factor at +1/-1 while the remaining factors take random levels from {-1, 0, +1}), and two corner points at the all-high and all-low settings. This ensures every factor spans all three levels and main effects are estimable with minimal aliasing.

**CWESS Statistic.** Factor importance is quantified using the Cumulative Weighted Effect Screening Statistic:

```
CWESS_i = |beta_hat_i| / (SE(beta_hat_i) + epsilon) * sqrt(SNR)
```

where coefficients are estimated via ElasticNet regression (lambda = 0.01, alpha = 0.5), the standard error uses the corrected formula `SE = sqrt(MSE * diag(inv(X'X + lambda*I)))`, and the signal-to-noise ratio is `SNR = Var(X * beta_hat) / MSE`. The CWESS combines effect magnitude, estimation precision, and overall signal quality into a single screening metric.

**Interaction Scoring.** Two-factor interactions are scored as `IS_ij = |beta_hat_ij| * sqrt(SNR) * w_int`, where the interaction weight `w_int = sqrt((k-1)/2)` balances the relative number of interaction terms against main effects.

**Hybrid Classification.** Factors are classified using dual thresholds -- a percentile-based threshold (60th percentile of CWESS) and an absolute threshold (0.8 * median of non-zero coefficients). A factor is classified as:
- **Critical** if its CWESS exceeds the hybrid threshold or its coefficient exceeds the absolute threshold.
- **Moderate** if it falls between the lower (20th percentile) and upper thresholds.
- **Negligible** if it falls below both.

A failsafe promotes the top 2 factors to critical status if fewer than 2 are classified as critical but the maximum CWESS exceeds 0.5.

### Phase 2: Adaptive Selective Augmentation

**Goal:** Collect additional data tailored to the specific factor structure discovered in Phase 1.

This is where HASOD fundamentally differs from traditional sequential approaches. Instead of following a fixed protocol, Phase 2 selects its design strategy based on two quantities observed in Phase 1 -- the number of critical factors (*k_c*) and the number of significant interactions (*n_int*):

| Observed Structure | Strategy | Phase 2 Runs |
|---|---|---|
| Interactions detected, k_c <= 5 | Full Factorial on critical factors | 2^k_c |
| Interactions detected, k_c > 5 | Resolution V Fractional Factorial | 2^(k_c - 1) |
| No interactions, k_c <= 3 | RSM Star Points (alpha = sqrt(k_c)) | 2*k_c + 2 |
| Otherwise | Fold-Over of Phase 1 design | 2k + 3 |

When interactions are present, a factorial design cleanly estimates all two-factor interactions among critical factors. When interactions are absent but curvature may exist, axial (star) points with rotatability parameter `alpha = sqrt(k_c)` enable quadratic model fitting. The fold-over fallback provides de-aliasing when the factor structure is complex.

### Phase 3: GP-Based Optimization with Uncertainty-Guided Refinement

**Goal:** Locate the process optimum and reduce prediction uncertainty around it.

A Gaussian Process with a Matern 5/2 kernel is fitted to the combined Phase 1 + Phase 2 data. The predicted optimum is found via differential evolution (200 iterations, population size 25).

Rather than placing refinement points at fixed radial offsets from the predicted optimum (as in v3), HASOD v4 uses **uncertainty-guided sampling**: for each of the 6 refinement runs, 500 candidate points are generated around the optimum at progressively increasing radii, and the candidate with the highest GP posterior variance is selected. This ensures maximum information gain where the surrogate model is least certain, directly reducing prediction variance at the optimum.

The GP model is then refitted on all data from all three phases, yielding final predictions and factor classifications.

---

## Novelty and Contributions

### The Gap HASOD Fills

Existing DOE methods occupy distinct quadrants of a screening-optimization capability space:

```
                      Factor Screening Capability
                      Low                    High
                +----------------------------------+
        High    | LHS, Sobol, BO,        |         |
  Optimization  | Box-Behnken            |  HASOD  |
  Capability    +------------------------+---------+
        Low     | Random designs         | PB, DSD |
                |                        | CCD     |
                +----------------------------------+
```

No existing published method integrates adaptive factor screening with GP-based optimization in a principled sequential framework. HASOD is the first to occupy the high-screening + high-optimization quadrant.

### Why HASOD is Novel Compared to Existing Methods

**1. Adaptive three-phase framework (vs. fixed sequential protocols).**
Traditional sequential DOE (e.g., Plackett-Burman followed by Central Composite Design) uses a fixed two-stage protocol regardless of what the screening stage reveals. If the system turns out to be sparse, the CCD stage wastes runs on negligible factors. If interactions dominate, the fixed CCD may not resolve them. HASOD's Phase 2 adapts its strategy -- choosing between full factorial, fractional factorial, RSM star points, or fold-over -- based on the actual number of critical factors and detected interactions. This factor-structure-aware augmentation is not present in the standard DOE literature.

**2. CWESS screening statistic (vs. simple t-tests or correlation-based screening).**
Traditional screening relies on individual t-statistics or correlation coefficients, which treat each factor in isolation. CWESS integrates three sources of information -- regularized coefficient magnitude (via ElasticNet), estimation precision (via corrected standard errors), and overall signal quality (via SNR) -- into a single metric. The hybrid classification using both percentile and absolute thresholds provides robustness across diverse factor structures. Theorems 1 and 2 establish that CWESS provides asymptotic separation between active (O(sqrt(n))) and inactive (O(1)) factors, guaranteeing classification consistency.

**3. Uncertainty-guided optimization refinement (vs. fixed-grid or random refinement).**
Standard response surface methodology places refinement points on fixed geometric patterns (star points, face-centered designs). Bayesian Optimization adapts point placement but provides no factor screening. HASOD's Phase 3 uses GP posterior variance to select refinement points where the model is most uncertain, combining the interpretability of DOE with the sample-efficiency of Bayesian methods. Theorem 4 proves that this approach monotonically reduces posterior variance at the predicted optimum.

**4. Rigorous theoretical guarantees (vs. heuristic-only methods).**
HASOD is supported by four theorems with complete proofs:
- **Theorem 1 (CWESS Separation):** Active factors produce CWESS = O(sqrt(n)); inactive factors produce CWESS = O(1).
- **Theorem 2 (Classification Consistency):** Under threshold conditions tau_n -> infinity and tau_n = o(sqrt(n)), misclassification probability goes to zero.
- **Theorem 3 (Adaptive Run Efficiency):** Total runs adapt to factor structure, achieving 64-88% reduction over traditional CCD for sparse systems.
- **Theorem 4 (GP Variance Reduction):** Phase 3 refinement provably reduces posterior variance at the predicted optimum via the information-never-hurts property of Bayesian conditioning.

Most screening methods (PB, DSD) lack formal detection guarantees, and most optimization methods (BO, LHS) provide no screening theory. HASOD provides both.

**5. Unified workflow (vs. separate screening and optimization studies).**
In practice, experimenters must currently run a screening study, analyze results, design an optimization study, and then run it -- often with weeks or months between stages and with the risk of losing context. HASOD integrates both goals into a single automated workflow that typically requires 35-51 total runs, compared to 81+ for a traditional CCD on 6 factors.

### Comparison to Specific Methods

| Aspect | Traditional (PB+CCD) | DSD | Bayesian Optimization | HASOD |
|---|---|---|---|---|
| Factor screening | Yes (fixed) | Yes (limited) | No | Yes (adaptive) |
| Interaction estimation | Partially | Some | No | Yes (via Phase 2 strategy) |
| Optimization capability | Yes (fixed model) | Limited | Yes (GP-based) | Yes (GP-based) |
| Adapts to observed structure | No | No | Yes (acquisition function) | Yes (strategy selection) |
| Theoretical guarantees | Partial | Partial | Regret bounds | Detection + optimization guarantees |
| Detection accuracy | 83.33% | 73.61% | 0% | 97.08% |

---

## Benchmark Results

HASOD v4 was benchmarked against 8 competitor methods across 6 test scenarios (sparse, dense, moderate, interaction-heavy, quadratic-heavy) with 10 replications each (540 total experimental configurations).

### Overall Performance

| Method | Detection Accuracy | Prediction Error | Avg Runs | Time (s) |
|---|---|---|---|---|
| **HASOD v4** | **97.08%** | 3.611 | 41.5 | 2.374 |
| Traditional (PB+CCD) | 83.33% | 4.312 | 29.0 | 0.005 |
| Augmented DSD | 82.22% | 3.496 | 25.0 | 0.003 |
| DSD | 73.61% | 3.618 | 13.0 | 0.002 |
| Box-Behnken | 0.00% | 2.124 | 63.0 | 0.001 |
| LHS | 0.00% | 0.834 | 17.0 | 0.113 |
| Sobol | 0.00% | 1.450 | 17.0 | 0.131 |
| Seq-D-Optimal | 0.00% | 2.456 | 18.0 | 0.013 |
| Bayesian Optimization | 0.00% | 20.133 | 30.0 | 10.166 |

Methods without factor classification (LHS, Sobol, Box-Behnken) can achieve lower prediction error but provide zero factor detection capability. Among methods that perform screening, HASOD v4 outperforms the best competitor by 13.75 percentage points.

### Robustness Across Scenarios

| Scenario | Detection Accuracy | Prediction Error | Total Runs |
|---|---|---|---|
| Sparse (few factors) | 100.0% | 3.526 | 44.2 |
| Sparse (many factors) | 100.0% | 4.130 | 37.0 |
| Quadratic-heavy | 100.0% | 7.236 | 40.2 |
| Moderate | 97.5% | 2.707 | 41.8 |
| Dense | 95.0% | 2.988 | 51.4 |
| Interaction-heavy | 90.0% | 1.077 | 34.6 |

HASOD maintains >= 90% detection accuracy across all six scenarios, demonstrating robustness to unknown factor structures. Run counts adapt automatically (34.6 to 51.4) based on detected complexity.

### Statistical Significance

All detection accuracy improvements are statistically significant at alpha = 0.05:
- vs Traditional: t = 5.51, p < 0.001
- vs DSD: t = 7.71, p < 0.001
- vs Augmented DSD: t = 6.23, p < 0.001

---

## Repository Structure

```
HASOD_DOE-main/
  HASOD_v4.ipynb                           # Full implementation and benchmark (Python)
  HASOD_Methodology_v4_Corrected.md        # Complete mathematical methodology with proofs
  HASOD_Strong_vs_Defensible_Claims.md     # Theorem correction rationale for peer review
  HASOD_v4_result.md                       # Benchmark results summary
  HASOD_v4_Benchmark_Results.csv           # Raw benchmark data (540 configurations)
  HASOD_v4_Complete_Benchmark_Visualization.png  # 9-subplot benchmark visualization
  novelty_analysis.md                      # Novelty assessment and industrial impact analysis
  hasod_flowchart_one_slide_from_tex.html  # HASOD framework flowchart (interactive HTML)
  HASOD v1.docx                            # Early draft manuscript
  KP Notes.docx                            # Working notes
  README.md                                # This file
```

---

## Getting Started

### Requirements

```
numpy
pandas
matplotlib
seaborn
scipy
scikit-learn
pyDOE3
scikit-optimize
```

### Running the Benchmark

Open `HASOD_v4.ipynb` in Jupyter Notebook or JupyterLab and run all cells sequentially. The notebook:

1. Loads all required libraries.
2. Defines the `HASOD_v4` class and 8 competitor method classes.
3. Creates 6 test scenarios with known ground-truth factor structures.
4. Runs the full benchmark (6 scenarios x 10 replications x 9 methods).
5. Generates results tables and a comprehensive 9-subplot visualization.

### Using HASOD on Your Own Problem

```python
from HASOD_v4 import HASOD_v4  # or copy the class from the notebook

def my_experiment(X):
    """Your response function. X is an (n_runs x k) array with values in [-1, 1]."""
    # Replace with actual experimental measurements
    return y  # array of n_runs responses

model = HASOD_v4(n_factors=6, noise_std=1.0, random_state=42)
results = model.fit(my_experiment)

print("Critical factors:", results['critical_factors'])
print("Total runs used:", results['total_runs'])
# results['gp_model'] can be used for prediction at new points
```

---

## References

1. Jones, B. and Nachtsheim, C.J. (2011). A class of three-level designs for definitive screening in the presence of second-order effects. *Journal of Quality Technology*, 43(1), 1-15.
2. Montgomery, D.C. (2017). *Design and Analysis of Experiments* (9th ed.). Wiley.
3. Rasmussen, C.E. and Williams, C.K.I. (2006). *Gaussian Processes for Machine Learning*. MIT Press.
4. Box, G.E.P. and Wilson, K.B. (1951). On the experimental attainment of optimum conditions. *Journal of the Royal Statistical Society B*, 13(1), 1-45.
5. Zou, H. and Hastie, T. (2005). Regularization and variable selection via the elastic net. *Journal of the Royal Statistical Society B*, 67(2), 301-320.
