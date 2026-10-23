# Plan: analytic gradient and Hessian for exact-likelihood models

**Status:** Phases 1–2 implemented · **Branch:** `plan/defm-gradient-hessian` · **Date:** 2026-10-21

## Motivation

Users fitting DEFMs through the `defm` R package report fits of up to 4 hours
for models with 5 binary outcomes, Markov order 1, and each outcome intercept
interacted with ~22 baseline covariates (K ≈ 122 parameters, ~8,800 rows).

A synthetic benchmark with that shape (2,000 individuals, 8,994 rows,
21 binary covariates, random data) gives:

| Spec                                  |   K | Unique supports / arrays | Time per `loglike_defm` |
|---------------------------------------|----:|-------------------------:|------------------------:|
| logit intercepts + motifs × site      |  17 |                 96 / 6,994 |                0.05 ms |
| + intercepts × site                   |  27 |                 96 / 6,994 |                0.10 ms |
| + intercepts × all covariates         | 122 |          **5,789** / 6,994 |               **23 ms** |

Two things compound:

1. **No derivatives.** `defm_mle()` hands only the objective to
   `stats4::mle()`, so `optim()` (L-BFGS-B) approximates each gradient with
   finite differences (≈ 2K = 244 evaluations, ~6 s here), and `mle()` then
   builds the Hessian numerically (`optim(..., hessian = TRUE)`, ≈ 4K² ≈ 60,000
   evaluations, ~23 min here).
2. **Expensive evaluations.** Covariate-interaction counters add the
   covariate values to the support hash, so almost every array gets its own
   support (83 % above), and each support stores all 2⁵ configurations × K
   statistics.

This plan addresses (1), which is the largest and most general win: it applies
to every barry model whose likelihood is computed by full enumeration of the
support (DEFM, and the ERGM-style models used by `ergmito`). Item (2) is listed
as follow-up work at the end.

## Prior art: `ergmito`

`ergmito` (same author) already implements the exact gradient and Hessian for
the same likelihood form:

- `ergmito/src/likelihood.h`, `ergmito_ptr::exact_gradient()`: observed minus
  expected statistics, reusing the cached `exp(statmat * params)` from the
  normalizing-constant update.
- `ergmito/src/ergmito_ptr.cpp`, `exact_hessiani()`: per-network
  `-(Σ w z t tᵀ · Z − (Σ w z t)(Σ w z t)ᵀ) / Z²`, summed over networks.

We port the math, not the code: `ergmito` is Armadillo-based, while barry is
dependency-free and stores supports as a flat `std::vector<double>`. Three
differences matter:

- barry de-duplicates supports, so each support's contribution must be
  weighted by the number of arrays that share it (`stats_support_n_arrays`).
- `exact_hessiani()` exponentiates without any shift, so it can overflow.
  barry should shift by the per-support maximum instead.
- `ergmito` supports offsets and barry does not, so we can drop them.

## The math

The reference version of this section, with implementation details, is
[`doc/likelihood-derivatives.md`](../doc/likelihood-derivatives.md).

### Notation

- Arrays (observations) $a = 1, \dots, A$, with observed sufficient statistics
  $s_a \in \mathbb{R}^K$ (`stats_target[a]`).
- Unique supports $u = 1, \dots, U$. Denote support $u$ by $\mathcal{Y}_u$; it is a list of $n_u$ distinct
  statistic vectors $t_{u1}, \dots, t_{un_u} \in \mathbb{R}^K$ with
  multiplicities (weights) $w_{u1}, \dots, w_{un_u}$, stored row-wise in
  `stats_support` starting at `stats_support_sizes_acc[u] * (K + 1)`, each row
  being $(w_{ui}, t_{ui}^\top)$.
- $m_u$ is the number of arrays mapped to support $u$
  (`stats_support_n_arrays[u]`). Then $\sum_u m_u = A$.
- Parameters $\theta \in \mathbb{R}^K$.

### Log-likelihood (already implemented in `likelihood_total`)

$$
Z_u(\theta) = \sum_{i=1}^{n_u} w_{ui} \exp\!\left(\theta^\top t_{ui}\right),
\qquad
\ell(\theta) = \sum_{a=1}^{A} \theta^\top s_a \;-\; \sum_{u=1}^{U} m_u \log Z_u(\theta).
$$

### Gradient (score)

Let $p_{ui}(\theta) = \Pr_{\mathcal{Y}_u}(Y = t_{ui}; \theta) = w_{ui}\exp(\theta^\top t_{ui}) / Z_u(\theta)$, the model
probability of configuration $i$ within support $u$, and
$\mu_u(\theta) = \sum_i p_{ui}\, t_{ui} = \mathbb{E}_{\mathcal{Y}_u}[t; \theta]$.
Because $\partial \log Z_u / \partial \theta = \mu_u$,

$$
\nabla \ell(\theta) = \sum_{a} s_a \;-\; \sum_{u} m_u\, \mu_u(\theta).
$$

The first term is constant in $\theta$ and can be precomputed once at
`init()`: $S = \sum_a s_a$.

### Hessian (observed = expected information)

Because $\partial \mu_u / \partial \theta^\top = \operatorname{Cov}_{\mathcal{Y}_u}[t; \theta]$,

$$
\nabla^2 \ell(\theta) = -\sum_{u} m_u\, \Sigma_u(\theta),
\qquad
\Sigma_u(\theta) = \sum_i p_{ui}\, (t_{ui} - \mu_u)(t_{ui} - \mu_u)^\top .
$$

Properties that we rely on:

- $\nabla^2 \ell$ is negative semi-definite, so $\ell$ is concave and the MLE
  (when it exists) is unique. Newton–Raphson is a valid optimizer.
- In a canonical exponential family the observed and expected information
  coincide: $I(\hat\theta) = -\nabla^2 \ell(\hat\theta)$ and
  $\widehat{\operatorname{Var}}(\hat\theta) = I(\hat\theta)^{-1}$.
- Using the centered form $(t - \mu)(t - \mu)^\top$ rather than
  $\mathbb{E}[tt^\top] - \mu\mu^\top$ avoids catastrophic cancellation
  (ergmito uses the uncentered form).

### Numerically stable evaluation

Per support, with $\eta_{ui} = \theta^\top t_{ui}$ and
$c_u = \max_i \eta_{ui}$:

$$
\tilde z_{ui} = w_{ui} \exp(\eta_{ui} - c_u), \quad
\tilde Z_u = \sum_i \tilde z_{ui}, \quad
\log Z_u = c_u + \log \tilde Z_u, \quad
p_{ui} = \tilde z_{ui} / \tilde Z_u .
$$

This replaces the fixed `BARRY_SAFE_EXP` (−100) shift for the new code paths
and is exact (the shift cancels in $p_{ui}$ and is added back in $\log Z_u$).

### Cost

For one pass over all supports, with $N = \sum_u n_u$ support rows:

| Quantity              | Cost            | Benchmark above (N ≈ 185k, K = 122) |
|-----------------------|-----------------|-------------------------------------|
| $\ell$                | $O(NK)$         | ≈ 23 ms (measured)                  |
| $\ell$ and $\nabla\ell$ (fused) | $O(NK)$ | ≈ 2× $\ell$                     |
| $\nabla^2\ell$        | $O(NK^2)$       | ≈ 1–3 s, computed once at the MLE    |

Compared with finite differences, each gradient goes from ≈ 2K evaluations to
≈ 2, and the final Hessian from ≈ 4K² evaluations to one $O(NK^2)$ pass. For
the K = 122 model that is roughly 6 s → 50 ms per gradient and 23 min → ~2 s
for the Hessian.

### Sanity checks the implementation must pass

- $\nabla\ell$ and $\nabla^2\ell$ match central finite differences of
  `likelihood_total` (relative error < 1e-6) on several random $\theta$.
- A DEFM with only `counter_logit_intercept` terms (with and without
  covariates) and order 0 is a set of independent logistic regressions:
  gradient and Hessian must match the closed-form logistic score
  $X^\top(y - \hat p)$ and $-X^\top \operatorname{diag}(\hat p(1-\hat p)) X$.
- The gradient is zero at the MLE, and $-\nabla^2\ell$ is symmetric positive
  definite there.
- Results are identical with `ncores = 1` and `ncores > 1`, and with
  `transform_model()` active (the derivatives are with respect to the
  parameters of the transformed statistics, which are what is stored).

## Implementation plan

### Phase 1: barry core (`include/barry/model-*.hpp`): **implemented**

1. **Free helpers** next to `update_normalizing_constant()` in
   `model-meat.hpp`:

   ```cpp
   // Fills probs (p_i) and mu (E[t | support]); returns log Z (-Inf if empty).
   inline double update_expected_stats(
       const std::vector< double > & params,
       const double * support, size_t k, size_t n,
       std::vector< double > & probs, double * mu
   );

   // Adds -weight * Cov[t | support] to the lower triangle of H (row-major).
   inline void add_covariance(
       const double * support, size_t k, size_t n,
       const std::vector< double > & probs, const double * mu,
       double weight, double * H, double * dev
   );
   ```

2. **`Model` methods** (declared in `model-bones.hpp` beside
   `likelihood_total`):

   ```cpp
   double likelihood_total_and_grad(
       const std::vector< double > & params,
       std::vector< double > & grad, BARRY_NCORES_ARG(=1)
   ) const;

   std::vector< double > likelihood_total_grad(
       const std::vector< double > & params, BARRY_NCORES_ARG(=1)
   ) const;

   // Row-major K x K.
   std::vector< double > likelihood_total_hess(
       const std::vector< double > & params, BARRY_NCORES_ARG(=1)
   ) const;
   ```

   - The methods are parallelized over supports with OpenMP, using
     per-thread accumulators reduced in a critical section. They fall back to
     one thread when there are fewer than 128 supports.
   - Supports with no rows contribute $-\infty$ to $\ell$ and nothing to the
     derivatives.
   - Differences from the original proposal:
     - The methods are `const` and do **not** touch the
       `normalizing_constants` / `params_last` cache. This keeps them thread
       safe and avoids mixing the stable $\log Z$ with the `BARRY_SAFE_EXP`
       convention. The fused function already returns $\ell$, so the cache
       isn't needed.
     - $S = \sum_a s_a$ is recomputed on each call rather than cached. It costs
       $O(AK)$, negligible next to the support pass, and needs no
       invalidation logic.

3. **Tests** (`tests/22-likelihood-gradient-hessian.cpp`, in
   `tests/main.cpp`). Results on macOS clang, no OpenMP:

   | Check | Result |
   |---|---|
   | DEFM (K = 11, 180 supports) with covariate interactions, a motif, and a covariate-interacted transition term, vs. central differences | gradient 4e-7, Hessian 4e-9 |
   | Order-0 logit-intercept DEFM vs. closed-form logistic score/Hessian | ~1e-15 |
   | Newton–Raphson from zero | converges to $\lVert\nabla\ell\rVert_\infty$ ≈ 6e-14 |
   | ERGM (3 networks) vs. central differences, before and after `set_transform_model()` | ~1e-9 |
   | Other properties | Hessian symmetric, $-\nabla^2\ell$ positive definite; `ncores = 2` matches `ncores = 1`; wrong parameter length throws |

4. **Benchmark** (`benchmarks/01-defm-derivatives/`, simulated data).
   Model: 8,979 rows, 6,979 arrays, 6,979 unique supports, K = 122,
   1 core, macOS clang `-O3`.

   | Quantity | Time |
   |---|---:|
   | `likelihood_total` | 22.3 ms |
   | `likelihood_total_and_grad` | **16.4 ms** |
   | Gradient by central differences (2K evals) | 5,439 ms |
   | `likelihood_total_hess` | **234 ms** |
   | Hessian by finite differences (4K² evals) | ~22 min |

   The fused likelihood + gradient is faster than the existing likelihood
   alone because it walks each support row by row (contiguous) instead of
   column by column.

### Phase 2: `defm` R package: **implemented in [UofUEpiBio/defm#14](https://github.com/UofUEpiBio/defm/pull/14)**

defm gets the barry headers from the barry R package (`LinkingTo: barry`,
installed from USCbiostats/barryr). So the merge order is: this PR, then a
header sync into barryr, then defm#14.

1. Rcpp exports `loglike_defm_grad(m, par, ncores = 1)` and
   `loglike_defm_hess(m, par, ncores = 1)`, plus an `ncores` argument for
   `loglike_defm()`. Its default stays at 2, the current behavior.
2. `defm_mle()` calls `stats::optim()` directly with `fn` and `gr`, then
   builds the `mle` / `defm_mle` object with
   `vcov = solve(-loglike_defm_hess(m, coef))`. The class contract
   (`coef()`, `vcov()`, `summary()`, `logLik()`, `AIC()`) is unchanged.
   - New arguments: `method`, `fixed`, `ncores`, and `control`.
   - The stored call is `defm_mle(...)` itself, so `profile()` /
     `confint()` refits work. They failed before.
   - A singular Hessian gives an `NA` vcov with a warning.
3. `method = "newton"`: Newton–Raphson with step halving, falling back to
   steepest descent if the Hessian is singular.
4. `inst/tinytest/test_gradient.R`: 30 expectations covering finite
   differences, agreement with the previous fit and with `glm()`, Newton,
   BFGS, fixed parameters, and refits.

Result on the simulated K = 122 model: `defm_mle()` takes 2.0 s (L-BFGS-B,
54 iterations) or 1.3 s (Newton, 6 iterations). The previous approach is
estimated at ~23 min.

### Phase 3: follow-ups that reduce the cost per evaluation

These are independent of the derivatives and can be scheduled separately.

1. **Drop lagged rows from the support hash when no counter looks at them.**
   `keygen_defm` and `defm_hasher_factory` always include rows
   $0, \dots, M-1$. With only current-time terms, that splits each support up
   to $2^{J \cdot M}$ ways for nothing (≤ 4.5× fewer supports in the
   benchmark). A counter could declare whether it reads past rows.
2. **Exploit the DEFM structure.** For logit-intercept-type terms the linear
   predictor is $\sum_j y_j(\beta_j + x^\top\gamma_j)$. Computing the $J$
   per-outcome log-odds once per array and enumerating only the motif terms
   over the $2^J$ configurations would avoid storing and multiplying the full
   $2^J \times K$ support matrix.
3. **`store_psets()` is always on in the `DEFM` constructor.** The benchmark
   model keeps 185k powerset arrays in memory. Make it opt-in (only `simulate`
   and `sample` need it).

## Out of scope

- Changing the likelihood itself, the counters, or the R-level formula
  interface.
- MCMC-MLE or other approximations. The exact likelihood stays exact.
