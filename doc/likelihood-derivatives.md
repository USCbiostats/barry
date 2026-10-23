# Likelihood, gradient, and Hessian of exact models

barry fits discrete exponential-family models, including DEFMs and small
ERGMs, by *full enumeration*: for every observation it lists all the
configurations the observation could have taken (its **support**), so the
normalizing constant, and hence the likelihood, is computed exactly. This page
derives the log-likelihood, its gradient, and its Hessian, and explains how
`barry::Model` computes them.

## The model

Each observation (an *array*) $y_a$, $a = 1, \dots, A$, has sufficient
statistics $s(y_a) \in \mathbb{R}^K$. Its support, denoted \(\mathcal{Y}\), is a set of distinct
configurations $y_1, \dots, y_n$ with multiplicities (weights)
$w_1, \dots, w_n$, where $w_i$ counts the configurations sharing the
sufficient-statistic vector $s(y_i)$. Given parameters $\theta \in \mathbb{R}^K$, the probability of
the observed array is

```math
\mathrm{Pr}_{\mathcal{Y}}(Y = y_a; \theta) = \frac{\exp(\theta^{\mathbf{t}}s(y_a))}{\eta(\theta)},
\qquad
\eta(\theta) = \sum_{i=1}^{n} w_i \exp(\theta^{\mathbf{t}}s(y_i)).
```

Arrays that have the same support (as decided by the counters' hashers) share
it, so barry stores and evaluates each **unique support** $u = 1, \dots, U$
only once. Denote support $u$ by $\mathcal{Y}_u$. Let $m_u$ be the number of arrays that map to support $u$
($\sum_u m_u = A$).

## Log-likelihood

Assuming independent arrays,

```math
\ell(\theta) = \sum_{a=1}^{A} \theta^{\mathbf{t}}s(y_a) - \sum_{u=1}^{U} m_u \log \eta_u(\theta).
```

## Gradient

Define the probability of row $i$ of support $u$ and the expected statistics
within that support:

```math
p_{ui}(\theta) = \mathrm{Pr}_{\mathcal{Y}_u}(Y = y_{ui}; \theta) = \frac{w_{ui}\exp(\theta^{\mathbf{t}}s(y_{ui}))}{\eta_u(\theta)},
\qquad
\mu_u(\theta) = \sum_i p_{ui}\, s(y_{ui}) = \mathbb{E}_{\mathcal{Y}_u}[s(Y); \theta].
```

Since $\partial \log \eta_u / \partial\theta = \mu_u$, the gradient (score) is
the observed minus the expected statistics:

```math
\nabla \ell(\theta) = \sum_{a} s(y_a) - \sum_{u} m_u\, \mu_u(\theta).
```

## Hessian

Since $\partial \mu_u / \partial\theta^{\mathbf{t}} = \operatorname{Cov}_{\mathcal{Y}_u}[s(Y); \theta]$,

```math
\nabla^2 \ell(\theta) = -\sum_{u} m_u\, \Sigma_u(\theta),
\qquad
\Sigma_u(\theta) = \sum_i p_{ui}\,(s(y_{ui}) - \mu_u)(s(y_{ui}) - \mu_u)^{\mathbf{t}}.
```

The Hessian does not depend on the observed statistics, and it is negative
semi-definite. Consequences:

- $\ell$ is concave, so the MLE is unique when it exists, and Newton–Raphson
  ($\theta \leftarrow \theta - [\nabla^2\ell]^{-1}\nabla\ell$) is a natural
  optimizer.
- In a canonical exponential family the observed and expected (Fisher)
  information coincide, so the variance of the MLE is estimated by
  $[-\nabla^2 \ell(\hat\theta)]^{-1}$.
- The Hessian is singular when a statistic does not vary within any support,
  or when the MLE is on the boundary, e.g., under separation.

## Numerically stable evaluation

$\exp(\theta^{\mathbf{t}}s(y))$ easily overflows. Within each support, the derivative
code subtracts the largest linear predictor, $c_u = \max_i \theta^{\mathbf{t}}s(y_{ui})$
(over rows with positive weight):

```math
\tilde z_{ui} = w_{ui}\exp(\theta^{\mathbf{t}}s(y_{ui}) - c_u),
\qquad
p_{ui} = \frac{\tilde z_{ui}}{\sum_i \tilde z_{ui}},
\qquad
\log \eta_u = c_u + \log \sum_i \tilde z_{ui}.
```

The shift cancels exactly in $p_{ui}$ and is added back to $\log \eta_u$. The
Hessian uses the centered form $(s(Y) - \mu)(s(Y) - \mu)^{\mathbf{t}}$, which avoids the
cancellation that $\mathbb{E}[s(Y)s(Y)^{\mathbf{t}}] - \mu\mu^{\mathbf{t}}$ suffers from.
(`likelihood_total()` instead uses the fixed shift `BARRY_SAFE_EXP`.)

## Implementation

The supports live in `Model::stats_support` as one flat, row-major vector.
Row $i$ of support $u$ has $K + 1$ entries, the weight followed by the
statistics, $(w_{ui}, s(y_{ui})^{\mathbf{t}})$:

| Member | Content |
|---|---|
| `stats_support` | All support rows, $(w, s(y)^{\mathbf{t}})$, concatenated |
| `stats_support_sizes[u]` | Number of rows $n_u$ of support $u$ |
| `stats_support_sizes_acc[u]` | Index of the first row of support $u$ |
| `stats_support_n_arrays[u]` | $m_u$, the number of arrays sharing support $u$ |
| `stats_target[a]` | Observed statistics $s(y_a)$ |

If the model is transformed (`Model::set_transform_model()`), both the
supports and the targets are stored transformed. The derivatives are then
with respect to the parameters of the transformed statistics.

| Function | Computes | Cost |
|---|---|---|
| `Model::likelihood_total()` | $\ell(\theta)$ | $O(NK)$ |
| `Model::likelihood_total_and_grad()` | $\ell(\theta)$ and $\nabla\ell(\theta)$ in one pass | $O(NK)$ |
| `Model::likelihood_total_grad()` | $\nabla\ell(\theta)$ | $O(NK)$ |
| `Model::likelihood_total_hess()` | $\nabla^2\ell(\theta)$, row-major $K \times K$ | $O(NK^2)$ |

Here $N = \sum_u n_u$ is the total number of support rows. Two free functions
in `model-meat.hpp` do the per-support work:

- `update_expected_stats()` computes $p_{ui}$ and $\mu_u$ for one support,
  with the max shift, and returns $\log \eta_u$ (or $-\infty$ for an empty
  support).
- `add_covariance()` adds $-m_u\,\Sigma_u$ to the lower triangle of the
  Hessian.

The derivative methods are `const`: they neither read nor update the cache of
normalizing constants that `likelihood_total()` keeps. With OpenMP, they split
the supports across `ncores` threads, each with its own accumulator, and
merge the results at the end. They run on one thread when there are fewer
than 128 supports.

## Example

```cpp
#include "barry/models/defm.hpp"

// ... build a defm::DEFM `model`, add counters, and call model.init() ...
auto & m = model.get_model();

std::vector< double > theta(m.nterms(), 0.0), grad;
double ll = m.likelihood_total_and_grad(theta, grad); // log-likelihood + score
std::vector< double > H = m.likelihood_total_hess(theta); // row-major K x K
```

In R, the defm package exposes these as `loglike_defm_grad()` and
`loglike_defm_hess()`. Its `defm_mle()` passes the gradient to `optim()`
(or runs Newton–Raphson with `method = "newton"`), and computes the
variance-covariance matrix from the exact Hessian.

## Checks

`tests/22-likelihood-gradient-hessian.cpp` verifies:

- the gradient and Hessian against central finite differences, for a DEFM
  with covariate interactions and transition terms, and for an ERGM before
  and after a model transformation;
- an order-0 DEFM with only logit intercepts against the closed-form
  logistic-regression score $X^{\mathbf{t}}(y - \hat p)$ and Hessian
  $-X^{\mathbf{t}} \operatorname{diag}\{\hat p(1 - \hat p)\} X$;
- that Newton–Raphson converges to a zero gradient, and that the Hessian is
  symmetric, negative definite, and the same with one or more threads.

`benchmarks/01-defm-derivatives/` times these functions on a DEFM with
$K = 122$ terms.

## References

- Vega Yon, G. G., Pugh, M. J., & Valente, T. W. (2022). Discrete
  Exponential-Family Models for Multivariate Binary Outcomes.
  [arXiv:2211.00627](https://arxiv.org/abs/2211.00627).
- Vega Yon, G. G., Slaughter, A., & de la Haye, K. (2021). Exponential random
  graph models for little networks. *Social Networks*, 64, 225–238. (The
  `ergmito` R package implements the same exact gradient and Hessian.)
