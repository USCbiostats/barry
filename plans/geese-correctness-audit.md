# GEESE correctness audit and fix plan

**Status:** Fixed in this PR.

**Audited Barry revision:** `70901def12dd6a991b1f640c182467edc4e165c2`
(2026-10-29).

Paths are relative to `include/barry/`. The [reproducer](geese-audit/reproduce.cpp)
and its [output](geese-audit/results.txt) cover findings 1–9; a second build
on 2026-10-30 gave identical output. Findings 10–12 were added in review from
code inspection. The regression tests cover 10 and 11; the sampling-loop
bounds in 12 are not tested directly (they only matter on roundoff).

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

## Fixes (this PR)

All fixes are in this PR. Regression tests are in
`tests/24-geese-audit-fixes.cpp`; they compare against freshly built models,
`likelihood_exhaust`, and values computed before the fixes
(`results.txt`). The full suite passes (`make test-portable`).

- **1 and 10, `update_annotations`:** no longer calls `init_node(*parent)`, so
  updates do not add arrays to the model. It validates values (0, 1, or 9),
  updates the annotation counts, resets `subtree_prob` (and `log_scale`), and
  rebuilds `reduced_sequence`.
- **2:** an empty traversal returns `as_log ? 0.0 : 1.0`; `predict` returns
  early on an empty `preorder`.
- **3 and 12:** `Model` has a new `pset_params_last`, written only by
  `update_pset_probs` and checked by both `sample` overloads. Cumulative
  sampling loops, including the root draw in `Geese::simulate`, stop at the
  last state.
- **4 and 11:** the `Geese` copy constructor is deleted (Barry's tests, `Flock`,
  and the GEESE R package do not use it). `Model` copies keep `with_pset`.
- **5:** a stable `logistic()` helper is used at all five call sites.
- **6:** each node's `subtree_prob` is divided by its maximum, and
  `Node::log_scale` accumulates the log factors over the subtree. The log
  likelihood adds the root's `log_scale`; `predict_backend` normalizes the root
  posterior directly and divides by each node's own factor;
  `get_probabilities()` undoes the scaling.
- **7:** `keygen_full` stores one key element per parent-state bit.
- **8:** `likelihood_exhaust` uses `ord`.
- **9:** internal annotations are **ignored** everywhere: the exhaustive
  likelihood and prediction now treat them as missing, as pruning always did.
  Rejecting them was considered, but existing callers (for example,
  `tests/07-geese-flock.cpp`) pass internal annotations.

The reproducer documents the audited revision. It does not compile against
the fixed code because mode `basic` copies a `Geese` object.

### GEESE R package (follow-up)

- `src/geese-common.cpp` (`get_probabilities`) reads `node.subtree_prob`
  directly; with scaling it must call `Geese::get_probabilities()` instead.
- Add tinytest integration tests in `inst/tinytest/` for
  `predict_geese(..., leave_one_out = TRUE)` with default arguments, a Flock
  with one all-missing tree, likelihood followed by simulation, and internal
  annotations being ignored.
- Rerun any parameter-recovery results that used simulation after likelihood
  evaluation (finding 3).

## Not in scope

Expanding counter checks, Monte Carlo frequency tests, sanitizer matrices,
thread-safety documentation, and changes to Barry's shared normalization code.
None of these address a reproduced failure; open separate issues if needed.
