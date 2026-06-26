# Contributing to invasimapr

First off, thanks for taking the time to contribute! All contributions
are welcome, from reporting a bug to proposing a new feature or
improving the documentation.

## Report bugs or request features

Report bugs and request features at
<https://github.com/b-cubed-eu/invasimapr/issues>.

When reporting a bug, please include:

- Your operating system name and version, and your R session info
  ([`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)).
- A minimal, self-contained, reproducible example (a “reprex”).
- Detailed steps to reproduce the problem and what you expected to
  happen.

## Contribute code

We use the standard [GitHub
flow](https://docs.github.com/en/get-started/quickstart/github-flow):

1.  Fork the repository and create a branch from `main`.
2.  Make your changes, following the [B-Cubed software development
    guide](https://docs.b-cubed.eu/guides/software-development/) and the
    [tidyverse style guide](https://style.tidyverse.org/).
3.  Add or update unit tests (`testthat`) and roxygen2 documentation
    (`@return` and `@examples` for exported functions).
4.  Run
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html),
    [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
    and
    [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
    locally and make sure they pass without errors.
5.  Open a pull request describing your changes and referencing any
    related issues.

## Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://b-cubed-eu.github.io/CODE_OF_CONDUCT.md). By
participating in this project you agree to abide by its terms.
