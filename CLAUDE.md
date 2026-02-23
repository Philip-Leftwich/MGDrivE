# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

MGDrivE (Mosquito Gene Drive Explorer) is a simulation modelling framework for evaluating gene drive interventions for mosquito-borne disease control. This repository contains two R packages:

- **MGDrivE** (v1.6.0): Original package using tensor-oriented approach for population dynamics and genetic inheritance
- **MGDrivE2** (v1.0.1): Extended framework based on stochastic Petri nets, adding epidemiological dynamics and time-varying parameters

Website: https://marshalllab.github.io/MGDrivE/

## Scope of Work

**We are working exclusively with MGDrivE (v1). Ignore MGDrivE2 entirely unless explicitly instructed otherwise.**

All new code will be written as standalone scripts or functions that sit outside the package source. Do not modify any files under `MGDrivE/R/` or `MGDrivE/src/`. The package itself is to remain untouched. New functionality will be developed as supplementary code that calls into the existing package rather than extending it internally.

## Collaboration Style

This is an exploratory, discussion-first project. Before writing any code:
- Explain where relevant functionality currently lives in the codebase
- Describe how it works conceptually
- Discuss the approach and any trade-offs before implementing anything

Only proceed to write code when explicitly asked to do so.

## Current Objective

Understand how mating is currently implemented in MGDrivE and explore how multiple mating (i.e. females mating more than once) could be implemented as external code bolted onto the existing framework without modifying the package source.

## Building and Development

```R
# Install dependencies
install.packages(c("Rcpp", "R6", "data.table", "Rdpack", "roxygen2", "knitr", "rmarkdown"))  # MGDrivE
install.packages(c("Matrix", "deSolve", "ggplot2"))  # MGDrivE2

# Install packages from local source
devtools::install("./MGDrivE")
devtools::install("./MGDrivE2")

# Or install from CRAN
install.packages("MGDrivE")
install.packages("MGDrivE2")

# Rebuild documentation after changing roxygen comments
devtools::document("./MGDrivE")
devtools::document("./MGDrivE2")

# Build vignettes
devtools::build_vignettes("./MGDrivE2")

# Run R CMD check
devtools::check("./MGDrivE")
devtools::check("./MGDrivE2")
```

## Package Architecture

### MGDrivE - Gene Drive Cubes and Spatial Dynamics

Key file patterns:
- `MGDrivE/R/Cube-*.R`: Gene drive inheritance systems (CRISPR, Wolbachia, RIDL, MEDEA, ClvR, etc.)
- `MGDrivE/R/Network-*.R`: R6 class for landscape/metapopulation simulation
- `MGDrivE/R/Patch-*.R`: R6 class for individual patch dynamics
- `MGDrivE/src/*.cpp`: C++ distance functions and spatial kernels (via Rcpp)

The "cube" abstraction represents 3D inheritance tensors mapping (female genotype × male genotype → offspring genotypes).

### MGDrivE2 - Petri Net Framework

Key file patterns:
- `MGDrivE2/R/PN-*.R`: Petri net model definitions (places `-P.R`, transitions `-T.R`)
  - `lifecycle-*`: Mosquito-only dynamics
  - `epiSIS-*`: SIS epidemiological model
  - `epiSEIR-*`: SEIR epidemiological model
  - `-node-*`: Single patch
  - `-network-*`: Multi-patch with migration
- `MGDrivE2/R/equilibrium-*.R`: Equilibrium state calculations
- `MGDrivE2/R/sampling-*.R`: Numerical methods
  - `ODE`: Deterministic (uses deSolve)
  - `CLE`: Chemical Langevin Equation (stochastic approximation)
  - `DM`: Direct Method (Gillespie algorithm)
  - `PTS`: Poisson Time Stepping
- `MGDrivE2/R/hazard-*.R`: Event rate/hazard functions

The Petri net design separates model specification (places and transitions) from numerical methods, allowing independent extension of either.

### Vignettes (recommended reading order)

1. lifecycle-node: Single node mosquito dynamics
2. lifecycle-network: Metapopulation mosquito dynamics
3. epi-node: Single node with epidemiology (SIS)
4. epi-network: Network epidemiology
5. epi-SEIR: SEIR human dynamics
6. output-storage: CSV output handling
7. inhomogeneous: Time-varying parameters
8. advanced_topics: Custom numerical methods

## Key Design Notes

- MGDrivE2 uses MGDrivE's cube structures for genetic inheritance - both packages work together
- Petri nets are bipartite graphs: Places (P) hold tokens/state, Transitions (T) fire to move tokens
- Documentation uses roxygen2 with markdown enabled
- Both packages target CRAN submission (follow CRAN policies)

## Examples

- `Examples/SoftwarePaper/`: Gene replacement and suppression demos
- `Examples/SoftwarePaper2/`: Additional simulation scenarios

## Analysis Companion

Python package [MoNeT_MGDrivE](https://pypi.org/project/MoNeT-MGDrivE/) provides data analysis tools for MGDrivE output.
