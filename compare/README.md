# Simulation and comparison support code

This directory contains supporting implementations used to evaluate and
compare the statistical methods developed for `brmplus`. These files are
research support code and are not part of the package's public API.

Cluster launchers, scheduler configuration, and simulation entry-point scripts
are intentionally excluded from the distributed source package. Consequently,
the files retained here are not intended to be executed as standalone
workflows.

Package users should install `brmplus` and use the documented functions exposed
by its namespace. The code in this directory is retained only to document the
method implementations underlying the comparisons.
