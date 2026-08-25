# brmplus

`brmplus` implements extended binary regression models for risk-ratio and
risk-difference inference in small-sample and rare-event studies.

## Installation

Install the package from the repository root:

```sh
R CMD INSTALL .
```

The package requires R 3.5.0 or later, Rcpp, and RcppArmadillo.

## Basic use

```r
library(brmplus)
help(package = "brmplus")
```

## Simulation studies

The `compare/` directory contains method-comparison and simulation code.
The Cluster entry points (`run_*.R` and `.sh`) remain available in the GitHub
repository but are excluded from the built R package.

See [`compare/README.md`](compare/README.md) for dependencies and commands.
