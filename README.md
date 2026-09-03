# brmplus

`brmplus` implements extended binary regression models for risk-ratio and
risk-difference inference in small-sample and rare-event studies.

## Installation

Install the package from the repository root:

```sh
R CMD INSTALL .
```

The package requires R 3.5.0 or later, Rcpp, RcppArmadillo, and a C/C++
toolchain suitable for compiling R packages.

## Basic use

```r
library(brmplus)
help(package = "brmplus")
```

## Simulation studies

The Git repository contains the research workflows used for the simulation and
method-comparison studies. The standalone entry points are located in
`compare/`:

- `run_simulation.R`: the main method-comparison simulation;
- `run_rarediff.R`: rare-event stress scenarios;
- `run_brmfirth.R`: BRM, Firth, and exact-interval comparisons;
- `run_simulation_time.R`: computation-time comparisons.

These workflows can be rerun from a Git checkout after the required R packages
have been installed. They source supporting implementations from the same
directory and use `brmplus` installed from the same checkout.

The simulation entry points (`compare/run_*.R`) and cluster wrappers
(`compare/*.sh`) are excluded from the source package produced by
`R CMD build`. They are therefore available from the Git repository, but not
from an installed copy of `brmplus`.

The standalone workflows are research code rather than exported package
functions. Individual helper functions sourced by these workflows should not
be treated as a stable public API.

## Reproducibility

The simulation scripts use fixed replication seeds or registered random-number
streams. To identify the exact version of the study code, record the Git commit
used for an analysis:

```sh
git rev-parse HEAD
```

Third-party R package versions are not currently pinned. Exact numerical
agreement across systems may therefore depend on the R version, package
versions, compiler, and numerical libraries.
