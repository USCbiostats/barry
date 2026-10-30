# GEESE correctness audit and fix plan

**Status:** Proposed; no fixes in this PR.

**Audited Barry revision:** `70901def12dd6a991b1f640c182467edc4e165c2`
(2026-10-29).

Paths are relative to `include/barry/`. The [reproducer](geese-audit/reproduce.cpp)
and its [output](geese-audit/results.txt) cover findings 1–9; a second build
on 2026-10-30 gave identical output. Findings 10–12 were added in review from
code inspection and have not been executed yet.

```sh
c++ -std=c++17 -O0 -g -Iinclude plans/geese-audit/reproduce.cpp -o /tmp/barry-geese-audit
for mode in basic cache underflow hash internal indices counters loo-check; do
    /tmp/barry-geese-audit "$mode"
done
```

The fixture uses IDs `0..4`, parents `-1,0,0,1,1`, one function, all
duplications, parameters `(-2,-1,0.3)`, and 9 for missing annotations.

## Findings

| # | Problem | Where | Evidence (mode) |
|---|---|---|---|
| 1 | Annotating a tip in a previously all-missing subtree is ignored by the reduced likelihood, because `reduced_sequence` is never rebuilt | `geese-meat.hpp`, `update_annotations` | reduced 1 vs. full and fresh 0.4072 (`basic`) |
| 2 | An all-missing tree returns `-inf` on the log scale instead of 0 | `geese-meat-likelihood.hpp:219` | (`basic`) |
| 3 | `simulate` reuses powerset probabilities computed by `likelihood` at other parameters, because the global `update_pset_probs` does not update `params_last` | `model-meat.hpp:354`, `Model::sample` | 4279 vs. 5718 out of 10,000 (`cache`) |
| 4 | Copied nodes keep parent pointers into the source; updating the copy corrupts the source | `geese-meat-constructors.hpp:230` | source throws afterwards (`basic`) |
| 5 | Root logistic `exp(p)/(exp(p)+1)` overflows to NaN | likelihood, likelihood_exhaust, predict, predict_exhaust, simulate | logit 1000 gives NaN (`basic`) |
| 6 | Log likelihood underflows on large observed trees because pruning messages are unscaled | `geese-meat-likelihood.hpp` | 1,100 tips give `-inf`, expected −762.46 (`underflow`) |
| 7 | `keygen_full` packs parent states as decimal digits in a `double`, so keys collide at 17+ functions | `geese-bones.hpp:36` | equal keys (`hash`) |
| 8 | `likelihood_exhaust` indexes columns by node `id` instead of `ord` | `geese-meat-likelihood_exhaust.hpp:103` | throws with IDs 10, 20, 30 (`indices`) |
| 9 | Pruning ignores annotations on internal nodes; the exhaustive functions condition on them | pruning vs. exhaust | 0.0781 vs. 0.0722 (`internal`) |
| 10 | Each `update_annotations` call runs `init_node(*parent)`, which appends new arrays to the (shared) model; LOO does this twice per annotated tip, so the model grows on every LOO call | `geese-meat.hpp`, `update_annotations` | inspection |
| 11 | `Model`'s copy constructor and assignment drop `with_pset` and `rengine`, so `simulate()` on a deep-copied `Geese` throws | `model-meat.hpp` ~518–620 | inspection |
| 12 | Sampling loops can read past the end on roundoff: the root draw in `Geese::simulate` and both `Model::sample` overloads; `predict(only_annotated = true)` on an all-missing tree reads `preorder[0]` of an empty vector | `geese-meat-simulate.hpp`, `model-meat.hpp`, `geese-meat-predict.hpp` | inspection |

The `counters` mode (384 configurations) and `loo-check` mode found no problems.

## Fix plan

Two Barry PRs, then one GEESE PR. Every fix gets one Catch regression test,
built from its reproducer mode, in the existing `tests/06-geese.cpp` …
`tests/10-geese-predict.cpp` files. The fixed `likelihood_exhaust` is the
reference for small trees; no separate enumerator is needed.

### PR 1 — small fixes (all findings except 6)

Each item is a few lines.

- **1 and 10, `update_annotations`:** reject non-leaf nodes and values other
  than 0, 1, or 9. Replace the `init_node(*parent)` call with: set
  `annotations`, clear `reduced_sequence`, reset `subtree_prob` to 1 for every
  node, and call `calc_reduced_sequence()`. The rebuild is safe to drop because
  `keygen_full` ignores cell values and pruning reads leaf annotations directly.
  Test: an updated model matches a freshly built one (missing→observed and
  observed→missing), and the model's array count is unchanged after LOO.
- **2:** return `as_log ? 0.0 : 1.0`; make `predict` return early on an empty
  `preorder`.
- **3:** add a `pset_params_last` vector to `Model`, written only inside
  `update_pset_probs` for every block it recomputes (including `i == -1`) and
  read only by the two `sample` overloads. This is separate from `params_last`,
  which also tracks normalizing constants. Test: `simulate(a)`,
  `likelihood(b)`, `simulate(a)` with a fixed seed matches a fresh model.
- **4 and 11:** delete the `Geese` copy constructor (`= delete`, as the
  assignment operators already are). Barry's tests, `Flock` (which uses move),
  and the GEESE R package do not use it. If something turns out to need it,
  fall back to: rebind `parent` pointers, copy `with_pset` in `Model`, and call
  `model->set_rengine`.
- **5:** add `inline double logistic(double p)` (`1/(1+exp(-p))` for
  `p >= 0`, otherwise `exp(p)/(1+exp(p))`) and use it at the five call sites.
- **7:** in `keygen_full`, `push_back` one element per parent-state bit instead
  of summing powers of 10.
- **8:** `offspring[o]->id` → `offspring[o]->ord`.
- **9:** reject non-missing internal annotations in the constructor (and, via
  item 1, in `update_annotations`). Conditioning on ancestral annotations can be
  a separate feature later if it is needed.
- **12:** stop each cumulative-probability loop at the last state.

### PR 2 — numerical underflow (finding 6)

This only matters once a tree has on the order of a thousand informative tips,
so it can follow PR 1. Per node, divide `subtree_prob` by its maximum and store
that factor's log in a new `log_scale` member. The log likelihood adds the
scales of all included nodes. `predict_backend` normalizes the root posterior
directly and divides the `above` term by the parent's own scale factor (all
other scales cancel). Test: 1,100 fair-coin tips give `-1100 * log(2)`, and
small trees match `likelihood_exhaust` to 1e-12.

### GEESE R package

Vendor each Barry fix as it merges. Add tinytest integration tests in
`inst/tinytest/` for: `predict_geese(..., leave_one_out = TRUE)` with default
arguments, a Flock with one all-missing tree, likelihood followed by
simulation, and the error on internal annotations. Rerun any parameter-recovery
results that used simulation after likelihood evaluation (finding 3).

## Not in scope

Expanding counter checks, Monte Carlo frequency tests, sanitizer matrices,
thread-safety documentation, and changes to Barry's shared normalization code.
None of these address a reproduced failure; open separate issues if needed.
