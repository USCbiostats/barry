# GEESE correctness audit and remediation plan

**Status:** Proposed; no algorithm fixes in this PR.

**Audited Barry revision:** `70901def12dd6a991b1f640c182467edc4e165c2`.

**Audit date:** 2026-10-29.

Eight correctness bugs and one public likelihood inconsistency were reproduced
against this Barry checkout. They affect annotation updates, missing-data log
likelihoods, probability-cache validity, object ownership, numerical stability,
support identity, and exhaustive indexing. The implementation plan below orders
repairs around independent reference calculations and regression tests.

## Evidence and reproduction

All finding paths below are relative to `include/barry/` and line numbers refer
to the audited revision. The standalone [reproducer](geese-audit/reproduce.cpp)
and [captured output](geese-audit/results.txt) accompany this report. It prints
observations and discrepancies; a successful exit is not an assertion that the
bugs are fixed. Counts from seeded simulations may vary between standard-library
implementations; deterministic cache-state tests are part of the proposed plan.

From the Barry repository root:

```sh
c++ -std=c++17 -O0 -g -Iinclude plans/geese-audit/reproduce.cpp -o /tmp/barry-geese-audit
for mode in basic cache underflow hash internal indices counters loo-check; do
    /tmp/barry-geese-audit "$mode"
done
```

The standard fixture uses IDs `0,1,2,3,4`, parents `-1,0,0,1,1`, one function,
all duplication events, gains/losses/root parameters `(-2,-1,0.3)`, and missing
annotations encoded as 9. Individual modes vary the fixture as described below.

The existing six test cases in `tests/06-geese.cpp` through
`tests/11-phylo-counts.cpp` passed **15 assertions** without shortening their
simulation loops. Separately, 384 small-state counter configurations and three
LOO-versus-fresh-model comparisons passed. These checks establish limited
coverage, not overall correctness. No biochemical interpretation or literature
novelty claims are evaluated in this report.

## Confirmed bugs

### 1. Annotation updates leave reduced traversal stale — high priority

**Location:** `models/geese/geese-meat.hpp:288`, `Geese::update_annotations`,
especially line 311; `Geese::calc_reduced_sequence`, line 361.

**Trigger:** Initialize an entirely unannotated tree, then annotate tip 3 as 1.
Run mode `basic`.

**Observed:** Reduced likelihood is **1**, full likelihood and a fresh equivalent
model both give **0.4071904237454682**.

**Cause and consequences:** Updating the tip rebuilds its parent's arrays but
never rebuilds `reduced_sequence`. An initially omitted subtree stays omitted.
This silently discards new evidence in default likelihood and can affect default
predictions. The all-missing example makes the discrepancy particularly clear;
adding evidence to an omitted branch of a partially observed tree has the same
structural problem.

**Correction:** Clear and recompute the reduced sequence when annotations change.
Reset omitted subtree messages to 1 so previously computed messages cannot leak
into a newly reduced traversal. `calc_reduced_sequence` currently appends, so
calling it without clearing is insufficient. Refresh annotation summary counts
as well. Regression: compare updated and reconstructed objects after transitions
between all missing, partially annotated, and fully annotated tips.

### 2. Entirely missing trees have incorrect log likelihood — high priority

**Location:** `models/geese/geese-meat-likelihood.hpp:219–220`,
`Geese::likelihood`.

**Trigger:** All tips missing; default reduced traversal. Mode `basic`.

**Observed:** Probability likelihood **1**, log likelihood **-inf**. Full traversal
returns approximately **0** on the log scale.

**Cause and consequences:** The empty-sequence branch explicitly returns negative
infinity for `as_log`. Integrating all latent states gives probability 1 and log
likelihood 0. In `Flock::likelihood_joint`, one such tree therefore changes an
otherwise finite joint log likelihood to negative infinity.

**Correction:** Return `as_log ? 0.0 : 1.0`. Test standalone and pooled inference.

### 3. Global probability updates invalidate sampling's cache bookkeeping — high priority

**Location:** `model-meat.hpp:354–363`, `Model::update_pset_probs`;
`model-meat.hpp:1750`, `Model::sample`.

**Trigger:** Simulate at parameters a, evaluate likelihood at different parameters
b, then simulate again at a. Mode `cache` uses a=(-8,-8,0.3), b=(8,8,0.3),
100 warmup simulations, then 10,000 draws with seed 123.

**Observed:** Tip 2 is present **4279/10000** times after the intervening likelihood;
a fresh model simulated at a with the same seed gives **5718/10000**.

**Cause and consequences:** The all-support update (`i=-1`) overwrites every
`pset_probs` block but only the single-support update records `params_last`.
Sampling sees the old a label and reuses probabilities computed at b. In this
example low gain/loss probabilities preserve the root state at a, whereas b
nearly flips it on a one-edge path. The discrepancy is not Monte Carlo evidence
alone: the cache label and data are explicitly inconsistent in the source.
Inference followed by simulation can therefore simulate a different model.
Shared supports make this relevant across trees too, though the executable
reproducer exercises one tree.

**Correction:** Associate each probability block with its actual parameter vector,
updating every affected cache label. Prefer separate cache validity for
normalizing constants and powerset probabilities: other Model methods also use
`params_last`. Regression should interleave likelihood, prediction, and simulation
across parameter vectors and verify probabilities or identical seeded draws.

### 4. Copied nodes retain parent pointers into the source object — high priority

**Location:** `models/geese/geese-meat-constructors.hpp:230–303`, copy constructor;
`nodes(model_.nodes)` at line 237; only offspring pointers are rewired at line 293.

**Trigger:** Copy an initialized model, then change a copied leaf annotation.
Mode `basic` uses tip annotations 2=1, 3=0, 4=1 and changes copied tip 3 to 1.

**Observed:** `copy.nodes[3].parent == &source.nodes[1]` is **true**. After the
copy update, evaluating the source throws **The requested support is out of range**.

**Cause and consequences:** `update_annotations` follows the copied node's
parent pointer into the source and rebuilds that source node using the copy's
support model. The source receives support indices belonging to another model.
Destroying the source would also leave dangling parent pointers in the copy;
that lifetime case was not separately executed under a sanitizer.

**Correction:** Rebind every non-null parent pointer by ID into the destination
node map, as well as offspring pointers. Audit RNG wiring separately. Regression:
updates in a copy must leave the original unchanged, and the copy must remain
usable after destruction of the source.

### 5. Root logistic transform overflows — medium/high priority

**Location:** `models/geese/geese-meat-likelihood.hpp:118–120`, `Geese::likelihood`;
same expression in `geese-meat-simulate.hpp:20`,
`geese-meat-predict.hpp:19`, and `geese-meat-likelihood_exhaust.hpp:18`.

**Trigger:** Finite root logit 1000, otherwise ordinary parameters and observed
tips. Mode `basic`.

**Observed:** Likelihood is **NaN**.

**Cause and consequences:** `exp(p)/(exp(p)+1)` becomes infinity/infinity. The
correct root probability is effectively 1. This can disrupt optimization and
posterior calculation. The executable directly tests likelihood; the shared
expression establishes the same arithmetic defect in the other named methods.

**Correction:** Use a sign-stable logistic function: for nonnegative p,
`1/(1+exp(-p))`; otherwise `exp(p)/(1+exp(p))`. Share the helper across methods.

### 6. Log likelihood still underflows through unscaled pruning — high priority

**Location:** `models/geese/geese-meat-likelihood.hpp:68`, `pset_loop`;
`Geese::likelihood`, lines 187–223 (message sum and final logarithm).

**Trigger:** A binary tree with 1,100 observed tips, one function, all duplication,
all transition and root parameters 0. Mode `underflow` constructs 2,199 nodes.

**Observed:** Log likelihood **-inf**; exact answer **-762.4618986159398**.

**Why incorrect:** Here each tip is an independent fair coin under the transition
model, so log likelihood is exactly -1100 log(2). Products underflow before
`as_log` is applied. This is a numerical correctness failure for a finite,
well-defined log likelihood, not evidence of impossible data. Flock's sum of log
likelihoods cannot repair underflow inside individual trees.

**Correction:** Normalize pruning messages and accumulate log scaling factors,
or implement log messages with log-sum-exp. Posterior propagation must use the
same scaling convention. Test against the analytic fair-coin result.

### 7. Support keys collide for distinct parent states — lower practical priority

**Location:** `models/geese/geese-bones.hpp:36–57`, `keygen_full`.

**Trigger:** 17 functions; compare parent state with only function index 16 present
against the same state plus function index 0. Mode `hash` calls the actual key
function directly, avoiding enumeration of an enormous support.

**Observed:** Distinct states produce **equal keys**.

**Cause and consequences:** Parent bits are encoded as decimal digits in a
`double`. At 10^16 the additional unit cannot be represented. Support reuse then
cannot distinguish parent states that may have different gain/loss statistics.
This confirms a key collision; a complete 17-function likelihood was not run,
and exponential support growth makes the trigger less practical today.

**Correction:** Encode individual Boolean states as separate key elements or use
an exact packed integer/bit-vector key with equality checks. Do not fix this by
changing the decimal scale alone.

### 8. Exhaustive likelihood confuses node IDs with column ordinals — medium priority

**Location:** `models/geese/geese-meat-likelihood_exhaust.hpp:103`,
`Geese::likelihood_exhaust`.

**Trigger:** Three-node tree with IDs 10,20,30; parents -1,10,10; annotations 9,1,0.
Mode `indices`.

**Observed:** Pruning gives **0.1576230597216791**. Exhaustive evaluation throws
**The column is out of range: 20 >= 3**.

**Cause and consequences:** `base` stores nodes in ordinal columns but transition
extraction reads `offspring[o]->id`. Even contiguous IDs in a different input
order can silently select the wrong column. This also compromises exhaustive
likelihood as a validation oracle.

**Correction:** Read `offspring[o]->ord`. Test sparse IDs and permuted row orders.

### 9. Pruning and exhaustive likelihood disagree about ancestral observations

**Location:** `models/geese/geese-meat-likelihood.hpp:42–68`, `pset_loop` only
checks leaf annotations; `Geese::likelihood` sums unrestricted internal states.
Contrast `geese-meat-likelihood_exhaust.hpp:35–42`, which fixes all annotated
nodes in its enumeration.

**Trigger:** Five-node model with annotations root=1, node1=9, tip2=1, tip3=0,
tip4=1. Mode `internal`.

**Observed:** Pruning **0.07809990947097682**, exhaustive **0.0722198208145861**.

**Interpretation:** The two public likelihood implementations compute different
events on identical accepted input. This is a confirmed API inconsistency.
Whether ancestral observations should be supported is a modeling/API decision;
if they are intended as observations, pruning silently ignores evidence.

**Correction:** Either consistently reject nonmissing internal annotations with
a clear error, or mask incompatible internal states in pruning, reduced traversal,
and posterior propagation. Exhaustive enumeration must follow the same contract.

## Checks that did not identify bugs

Mode `counters` independently counts gains, losses, complementary offspring
function partitioning, parent-gated neofunctionalization (symmetric and A-to-B),
and joint gains for two functions and two offspring. It checks **384** combinations
(4 parent states × 16 offspring states × 2 event types × 3 event selectors),
covering eight statistic outputs per combination. **Zero mismatches.** This
supports these specific binary-offspring formulas and duplication/speciation
selectors. It does not test every counter, all insertion paths, polytomies,
dynamic support restrictions, or statistical identifiability.

Mode `loo-check`: leave-one-out predictions for each of three annotated tips
match newly constructed models with the entire corresponding tip masked:

| Tip | LOO and fresh masked probability |
|---|---:|
| 2 | 0.5272873804128089 |
| 3 | 0.6846375092305400 |
| 4 | 0.3501804795632941 |

No standalone LOO bug was established by this check. It does not establish
correctness with dynamic support constraints, extreme parameters, or ownership
failures. Missing tips are explicitly marginalized in the inspected pruning code;
there is no basis here for claiming they are generally treated as absent.

## Suspected issues and remaining coverage

- `predict_backend` divides by `parent.subtree_prob[s]` without a zero-message
  guard (`geese-meat-predict.hpp`, the `above.push_back` expression). Structural
  zeros from restricted support could cause 0/0; no minimal executable case was
  established in this audit. Treat as suspected, not a confirmed prediction bug.
- `Model::sample` uses an unbounded cumulative-probability loop
  (`model-meat.hpp:1756` and 1766). A probability sum slightly below a draw could
  read past the block. No actual failing draw was established.
- Copy/move RNG ownership, repeated initialization, invalid topology handling,
  dynamic rules, and concurrent Flock operations need dedicated follow-up tests.
  The Geese copy constructor creates a new RNG but does not explicitly bind the
  copied Barry model to it; this was inspected but not reproduced as a failure.
- The audit directly tested single-tree cache corruption; it did not run a
  comprehensive pooled, parallel, multi-tree regression matrix.
- No claim of complete implementation correctness follows from the passing
  small-state checks. All other evolutionary counters remain unaudited here.

## Modeling limitations established from the implementation

`Geese::init_node` (`geese-meat.hpp:15`) builds branch-length vectors of ones.
The supplied Geese constructor accepts no branch-length argument. Thus this
interface implements event-based transition probabilities, not a branch-length
rate model. `Geese::likelihood` constructs the root distribution as a product of
function-specific Bernoulli probabilities; root dependencies are not represented
by that prior. These are model choices, not bugs by themselves.

The tested partitioning and innovation statistics count particular endpoint
configurations. They do not establish biochemical mechanisms, unique histories,
or identifiability. The literature and proposed publication framing were not
reviewed in this implementation audit.

Prioritize fixes 1–6 before relying on annotation updates, likelihood optimization,
or simulation-based parameter recovery. Repair exhaustive indexing before using
that implementation as a reference. Then expand exact small-tree checks to
polytomies, dynamic constraints, event mixtures, sparse annotations, and Flock.

## Remediation plan

This is a proposed sequence of implementation PRs. Checkboxes denote future work,
not fixes delivered by this documentation PR. Each implementation PR should add
regressions with the fix; the standalone reproducer deliberately prints baseline
failures rather than serving as a passing CI test. Resolve the internal-annotation
contract in Phase 1 before implementing its behavior.

### Phase 1 — establish reference behavior and observation contracts

- [ ] Convert the existing reproductions into focused Catch tests in
  `tests/06-geese.cpp`, `tests/07-geese-flock.cpp`, and `tests/10-geese-predict.cpp`,
  with additional dedicated test files if needed. Register new files in
  `tests/main.cpp`. Keep a failing test beside each subsequent fix rather than
  merging a permanently failing default suite.
- [ ] Fix exhaustive likelihood's `id` versus `ord` access (finding 8). Test
  nonconsecutive IDs, input row permutations, and shuffled ID-to-row assignments.
- [ ] Build an independent tiny-tree enumerator for validation. Enumerate node
  states, form each transition's statistics with direct formulas, normalize its
  support, multiply transitions and the root prior, and sum states compatible
  with observations. Do not reuse pruning messages or its support lookup.
  Keep fixtures within tractable sizes: one/two functions, two to five nodes.
- [ ] Decide whether ancestral annotations are observations or unsupported input
  (finding 9). If supported, enforce compatibility at every internal/root state
  in likelihood and prediction and include annotated ancestors in reduced
  traversal. If unsupported, validate both construction and updates consistently;
  reject such input rather than silently ignoring it. Document the decision in
  the public API and in the GEESE R interface on integration.

**Acceptance:** The repaired exhaustive algorithm, independent oracle, and full
pruning agree on supported input within absolute tolerance 1e-12 for moderate
small-model parameters. ID and row-order permutations leave results unchanged.
Both algorithms consistently condition on or reject ancestral observations.

### Phase 2 — repair annotation state and graph ownership

- [ ] Recompute reduced traversal after annotation changes (finding 1). Clear
  `reduced_sequence` before rebuilding; reset messages for newly excluded nodes
  to their integrated value of 1. Recompute annotation summary counts. Ensure
  repeated updates do not append duplicate traversal entries.
- [ ] Make updating annotations atomic with respect to validation and support
  construction: validate values and dimensions first, and avoid leaving partially
  updated annotations if rebuilding throws. Confirm whether rebuilding parent
  arrays is actually necessary when only the observation mask changed.
- [ ] Rebind every copied parent and offspring pointer into the copied node map
  (finding 4). Audit Geese copy/move/assignment and Barry model RNG binding;
  verify deep copies own their RNG/model while borrowed Flock objects do not.
  Keep shared support ownership explicit and avoid changing Flock sharing into
  accidental deep copying.
- [ ] Test updated versus reconstructed models for missing-to-observed,
  observed-to-missing, 0-to-1, and partial multi-function updates. Exercise a
  previously excluded branch of a partially observed tree, not only a completely
  unobserved tree. Repeat both full and reduced calls in alternating order.
- [ ] Test source isolation after copy updates, copy survival after source
  destruction, and Flock vector growth/moves. Run these under AddressSanitizer
  and UndefinedBehaviorSanitizer.

**Acceptance:** Updated and fresh likelihoods/posteriors match; source objects
remain unchanged after copy updates; no stale/dangling pointers are detected.

### Phase 3 — make cache validity explicit across all shared-model callers

- [ ] Separate validity of normalizing constants from validity of powerset
  probabilities (finding 3). For each support block, record the parameter vector
  or generation that actually produced each cached artifact. A normalizer update
  alone must not certify an older probability block.
- [ ] Update labels for every block touched by a global update; invalidate or
  initialize labels when adding support. Inspect all `params_last` readers and
  writers, including likelihood, gradient/Hessian, exact sampling, approximate
  sampling, copy/move, and model transformations. Do not patch only the Geese
  simulation caller.
- [ ] Add a deterministic regression inspecting the actual powerset probabilities
  after interleaving `simulate(a)`, `likelihood(b)`, and `simulate(a)`. Supplement
  it with identical-seed comparisons to a fresh model. Verify every support
  block, not just the last one sampled.
- [ ] Repeat across two Flock trees sharing support and two trees introducing
  distinct support shapes or event types. Include prediction between sampling
  calls. Establish and document whether simultaneous external calls on a shared
  model are supported; avoid assuming sequential cache fixes make them safe.

**Acceptance:** Each returned probability block matches direct normalization at
its requested parameters, regardless of call order or which tree requested it.
Existing gradient/Hessian and sampling tests pass without numerical changes
outside affected cache-reuse paths.

### Phase 4 — stabilize likelihood and propagate the convention to predictions

- [ ] Return log likelihood 0 for an entirely missing tree (finding 2). Test a
  Flock containing one observed tree and one entirely missing tree: adding the
  latter must not change the joint log likelihood.
- [ ] Introduce one stable logistic helper for root probabilities (finding 5).
  Use it in likelihood, exhaustive evaluation, prediction, and simulation.
  Test logits -1000, 0, and 1000 as well as ordinary values.
- [ ] Choose and document scaled probability messages or log-space messages for
  pruning (finding 6). With scaled messages, record each subtree's accumulated
  log scale, include child scales in the final log likelihood, and keep omitted
  subtrees at unit message and zero log scale. With log messages, use
  log-sum-exp and handle structurally impossible states explicitly.
- [ ] Update posterior propagation at the same time: its current ratios use raw
  subtree messages. Derive the corresponding normalized/log expression and
  avoid 0/0 for impossible parent states. Define behavior when the entire
  observed event has probability zero, rather than returning NaN posteriors.
- [ ] Audit transition normalization in shared Barry code. Replace fixed
  exponential shifts with per-support maximum shifts where necessary, retaining
  a log normalizer so finite log probabilities remain representable. Check
  compatibility with existing likelihood derivatives before changing stored
  normalizer semantics. This is follow-up robustness work, not an additional
  reproduced transition-overflow finding.
- [ ] Bound cumulative sampling loops and distinguish roundoff in a valid
  normalized distribution from an invalid distribution. Never silently send a
  materially missing probability mass to the last state.

**Acceptance:** The 1,100-tip fair-coin fixture returns finite log likelihood
`-1100 * log(2)` within 1e-10. Small-model results match the independent oracle.
Root logits of magnitude 1000 produce valid limits. Posterior state masses sum
to 1 and observed tips have their observed marginal states. Sampling produces
only valid states. Probability-scale likelihood may legitimately underflow to 0;
the log-scale result must remain finite for the analytic fixture.

### Phase 5 — use exact support identity and expand evolutionary-statistic checks

- [ ] Replace decimal parent-state packing in `keygen_full` with an exact key
  (finding 7). One straightforward option within the existing vector-key API is
  separate elements for dimensions, event type, and each state bit. Verify the
  downstream hash/equality implementation compares the complete key.
- [ ] Test the demonstrated 17-function collision directly without enumerating
  its huge support. Add distinct event-type and shape key fixtures. Determine
  whether dimensions, parent states, event type, and applicable rule/counter data
  exhaust the support identity needed for each supported configuration.
- [ ] Extend the 384 passing statistic fixtures to three/four offspring,
  mixed duplication/speciation trees, additional function pairs, and all relevant
  counters. Verify direct counts and incremental change statistics along
  multiple insertion/removal paths, especially for neo/subfunctionalization and
  joint dependencies. These are validation extensions, not established counter
  bugs.
- [ ] Add structural-zero/dynamic-rule fixtures and partial annotation fixtures.
  Compare full/reduced/exhaustive likelihoods and exact posteriors where feasible.
  Compare simulated small-tree state frequencies against exact probabilities
  with predeclared Monte Carlo tolerances and fixed seeds.

**Acceptance:** Distinct model-defining configurations do not alias support keys;
pooled and independently initialized trees agree. Counter checks match their
explicit endpoint-statistic definitions across the expanded fixtures.

### Phase 6 — release and downstream integration

- [ ] Run `make test-portable` plus the supported OpenMP test configuration;
  run sanitizer builds for ownership and indexing fixtures. Preserve the existing
  six GEESE tests and broader Barry tests for unrelated consumers.
- [ ] Record numerical compatibility, changed edge-case behavior, and any
  ancestral-annotation contract decision in release notes. State explicitly that
  endpoint configuration statistics are not biochemical mechanisms and this
  constructor does not model branch-length-dependent rates.
- [ ] Regenerate the distributable `barry.hpp` through the repository's existing
  `make barry.hpp` process after implementation changes, and verify modular and
  amalgamated headers behave identically.
- [ ] Vendor the fixed Barry revision into GEESE in a separate integration PR;
  rerun R-level likelihood, annotation-update, Flock, LOO, and simulation checks.
- [ ] Rerun parameter-recovery studies after fixes. Retain seeds, parameter
  vectors, tree fixtures, and library revisions. Results generated through the
  confirmed stale-cache/update paths require reassessment; this report does not
  establish that any specific published analysis was affected.

**Release gate:** All confirmed cases have automated regression coverage;
suspected issues are either reproduced and fixed or explicitly documented as
unresolved; the full supported test matrix passes; the GEESE integration pins the
corrected Barry commit. Track performance after correctness is established,
particularly memory use and full/reduced traversal costs.
