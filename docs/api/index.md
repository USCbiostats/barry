# API Reference

This section is generated from the Doxygen comments in the headers under
`include/barry/` by [mkdoxy](https://github.com/JakubAndrysek/MkDoxy).

## Main components

- **`barry::BArray` and `barry::BArrayDense`**: sparse and dense binary arrays.
- **`barry::Counters` and `barry::Rules`**: user-defined sufficient statistics
  and support constraints.
- **`barry::Support` and `barry::PowerSet`**: enumeration of the support of
  an array.
- **`barry::Model`**: discrete exponential-family models with exact
  likelihoods.
- **Modules**: network models (`barry::counters::network`), phylogenetic
  models (`geese`), and Discrete Exponential-Family Models (`defm`).

Use the sidebar to browse classes, namespaces, and files.

## Likelihood derivatives

`barry::Model` provides
[`likelihood_total_and_grad()`](../barry/classModel.md#function-likelihood_total_and_grad),
[`likelihood_total_grad()`](../barry/classModel.md#function-likelihood_total_grad),
and [`likelihood_total_hess()`](../barry/classModel.md#function-likelihood_total_hess).
See [Likelihood, gradient, and Hessian](../impl/likelihood-derivatives.md) for
the mathematical derivation and usage.
