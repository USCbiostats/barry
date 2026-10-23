# Mathematical notation

When writing or revising mathematical documentation, use support-indexed
probability notation: `\mathrm{Pr}_{\mathcal{Y}}` denotes probability under the
support \(\mathcal{Y}\). Use `\mid` (or `|` in inline notation) to separate
random conditioning variables, and separate fixed model parameters with a
semicolon. For example:

\[
\mathrm{Pr}_{\mathcal{Y}}(Y = y \mid X = x; \theta).
\]

For support-specific quantities, use the relevant support index, such as
`\mathrm{Pr}_{\mathcal{Y}_u}(Y = y; \theta)`. Apply the same convention to related
expectation and covariance notation where appropriate.

Use `\mathbf{t}` as the transpose marker (for example,
`\theta^{\mathbf{t}}`). Denote normalizing constants by `\eta`, not `Z`.
Write probabilities in terms of the observed outcome and its sufficient
statistics, for example
`\mathrm{Pr}_{\mathcal{Y}}(Y = y; \theta) = \exp(\theta^{\mathbf{t}}s(y)) / \eta(\theta)`.
