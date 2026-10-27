# Approximate sampling (MCMC)

`Model::sample()` draws arrays exactly. It needs the full support of the
array (and `store_psets()`), so it only works for small arrays.
`Model::sample_approx()` draws arrays with Markov chain Monte Carlo (MCMC)
instead, so it never computes the support. It is also the building block for
simulation-based estimation (MCMC-MLE and stochastic approximation), so the
chain state persists across calls and the sampler can stream sufficient
statistics instead of arrays.

## Model

For an array $y$ in the support $\mathcal{Y}$ with sufficient statistics
$s(y)$,

```math
\mathrm{Pr}_{\mathcal{Y}}(Y = y; \theta) = \frac{\exp\left(\theta^{\mathbf{t}}s(y)\right)}{\eta(\theta)},
```

where $\eta(\theta)$ is the normalizing constant. MCMC only needs ratios of
probabilities, so $\eta(\theta)$ cancels:

```math
\frac{\mathrm{Pr}_{\mathcal{Y}}(Y = y'; \theta)}{\mathrm{Pr}_{\mathcal{Y}}(Y = y; \theta)}
  = \exp\left(\theta^{\mathbf{t}}\left[s(y') - s(y)\right]\right).
```

When $y'$ and $y$ differ in a single cell, $s(y') - s(y)$ is the change
statistic of that cell. If the model has a transformation $T$
(`set_transform_model()`), the exponent is
$\theta^{\mathbf{t}}\left[T(s(y')) - T(s(y))\right]$.

## The tie/no-tie (TNT) sampler

The default method, `SampleApproxMethod::tnt`, is Metropolis–Hastings with the
tie/no-tie proposal, as in `ergm`'s `MH_TNT`. Let $D$ be the number of free
cells (cells the rules allow to change) and $E$ the number of free cells that
are currently non-zero ("ties"). Each step:

1. If $E > 0$, with probability $1/2$ picks a random tie. Otherwise (and
   always when $E = 0$), picks a random free cell.
2. If the cell is a tie, proposes removing it. If it is empty, proposes adding
   it. The change statistic is computed with the cell set to one, the same
   convention as the support enumeration.
3. Accepts the proposal $y \to y'$ with probability

```math
\min\left\{1,\;
  \exp\left(\theta^{\mathbf{t}}\left[s(y') - s(y)\right]\right)
  \frac{q(y' \to y)}{q(y \to y')}\right\},
```

where the probability of proposing cell $c$ from a state with $e$ ties is

```math
q_e(c) = \begin{cases}
  \frac{1}{2}\frac{\mathbb{1}\{c \text{ is a tie}\}}{e} + \frac{1}{2D} & e > 0, \\
  \frac{1}{D} & e = 0.
\end{cases}
```

The correction $q(y' \to y)/q(y \to y')$ is needed because a tie can be
proposed through both branches but an empty cell only through the second.
Without it, the chain would converge to the wrong distribution unless the
array happens to be about half full.

Proposals leading to states that fail the dynamic rules (`get_rules_dyn()`)
are rejected. Dynamic rules are evaluated on the raw statistics of the
proposed state (as DEFM and GEESE rules do); rules that inspect the array
itself are not supported.

## State and statistics stream

```cpp
auto state = model.sample_approx_init(y0);   // array + stats + free cells

std::vector< double > stats;                 // n x nterms(), row-major
model.sample_approx(state, theta, 1000u);    // burn-in: 1000 proposals
model.sample_approx(state, theta, n, thin, &stats);
```

- `SampleApproxState` holds the current array, its raw sufficient statistics,
  the free cells, the list of current ties and the number of proposals and
  acceptances. The statistics are only updated with change statistics; they
  are never recounted.
- `sample_approx(state, params, n, thin, stats_out)` makes `n * thin`
  proposals. When `stats_out` is not `nullptr`, it appends the statistics of
  the current state every `thin` proposals. These are on the same scale as
  `get_stats_target()` (after the transformation, if any).
- `params` can change between calls, and later calls resume from the current
  state.
- `sample_approx(y0, params, nsteps, stats_out)` is a shortcut that returns a
  single array after `nsteps` proposals.
- `DEFM::simulate_approx()` is the counterpart of `DEFM::simulate()`: it draws
  each row with `sample_approx()`, starting from its observed values.

## Towards simulation-based estimation

The stream of statistics is what simulation-based estimators need. With
$x_i = s(Y_i) - s(y_{\mathrm{obs}})$ sampled at $\theta_0$:

- **MCMC-MLE** approximates the log-likelihood ratio by importance sampling,

  ```math
  \ell(\theta) - \ell(\theta_0) \approx
    -\log \frac{1}{n}\sum_{i=1}^{n} \exp\left((\theta - \theta_0)^{\mathbf{t}}x_i\right),
  ```

  or by its lognormal approximation, and maximizes it over $\theta$.
- **Stochastic approximation** (Robbins–Monro) runs a few steps at
  $\theta_t$, then updates
  $\theta_{t+1} = \theta_t - a_t D^{-1}\left(s(Y_t) - s(y_{\mathrm{obs}})\right)$,
  and continues the same chain at $\theta_{t+1}$.

Both only need the statistics, a persistent state and the ability to change
$\theta$ between calls, which is what `sample_approx()` provides.
