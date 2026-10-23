# Mathematical notation

When writing or revising mathematical documentation, use support-indexed
probability notation: `\Pr_{\mathcal{Y}}` denotes probability under the
support \(\mathcal{Y}\). Use `\mid` (or `|` in inline notation) to separate
random conditioning variables, and separate fixed model parameters with a
semicolon. For example:

\[
\Pr_{\mathcal{Y}}(Y = y \mid X = x; \theta).
\]

For support-specific quantities, use the relevant support index, such as
`\Pr_{\mathcal{Y}_u}(Y = y; \theta)`. Apply the same convention to related
expectation and covariance notation where appropriate.
