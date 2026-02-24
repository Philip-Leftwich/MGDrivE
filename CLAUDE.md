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

An intial multi_mating.R script has been created with a preliminary implementation. The next steps are to review this code, identify any issues or improvements, and discuss how it integrates with the existing MGDrivE workflow. Evaluate whether the current multimating script allows differential remating rates that are conditional on mated male genotypes, and if not, discuss how this could be implemented.

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
# Multiple Mating Extension for MGDrivE

This extension adds multiple mating capabilities to MGDrivE simulations. Females can remate over their lifespan, with last-male sperm precedence.

## Features

- **Refractory period**: Configurable delay after mating before females become eligible to remate
- **Encounter-rate dependent remating**: Remating probability scales with male density
- **Dynamic choosiness**: Female preferences can adjust based on male genotype frequencies
- **Last-male precedence**: New matings completely replace previous sperm

## Files

- `multi_mating.R` - Core implementation (source this file)
- `example_simulation.R` - Working example with MGDrivE integration

## Quick Start

```r
library(MGDrivE)
source("multi_mating.R")

# After setting up MGDrivE network...

# Configure parameters
mm_params <- list(
  T_refractory = 5,        # days before eligible to remate
  rho = 0.05,              # base remating propensity (0-1)
  lambda = 0.001,          # encounter rate parameter
  stochastic = FALSE,      # TRUE for stochastic, FALSE for deterministic
  modulate_fn = modulate_choosiness_linear,
  modulate_params = list(strictness = 0.2)
)

# Initialize state from existing population
mm_state <- init_multiMating_fromPopulation(popFemale, mm_params$T_refractory)

# In daily simulation loop, after standard mating:
mm_state <- oneDay_multiMating(
  state = mm_state,
  popMale = popMale,
  base_eta = driveCube$eta,
  newly_mated = newly_mated,
  params = mm_params
)

# Get total mated females for analysis
total_mated <- get_total_mated(mm_state)
```

## Parameters

### Core Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `T_refractory` | Integer | Days after mating before female can remate |
| `rho` | Numeric (0-1) | Base daily probability of remating (when males abundant) |
| `lambda` | Numeric | Encounter rate; higher = remating less dependent on male density |
| `stochastic` | Logical | Use stochastic sampling (TRUE) or deterministic (FALSE) |

### Remating Probability

The daily probability that an eligible female remates:

```
p_remate = rho × (1 - exp(-lambda × total_males))
```

- At low male density: `p_remate ≈ rho × lambda × total_males`
- At high male density: `p_remate → rho`

### Choosiness Modulation Functions

| Function | Behaviour |
|----------|-----------|
| `modulate_choosiness_none` | Static preferences (use base eta unchanged) |
| `modulate_choosiness_linear` | Gradually relax toward uniform as preferred males become rare |
| `modulate_choosiness_threshold` | Binary: strict preferences if preferred males above threshold, else uniform |
| `modulate_choosiness_sigmoid` | Smooth S-curve transition between strict and relaxed |

#### Linear Relaxation Parameters
- `strictness`: Frequency threshold below which relaxation begins (default 0.2)

#### Threshold Parameters
- `threshold`: Minimum frequency of preferred males to maintain strict preferences (default 0.1)

#### Sigmoid Parameters
- `midpoint`: Frequency at which relaxation is 50% (default 0.15)
- `steepness`: How sharp the transition is (default 20)

## Data Structures

### MultiMating State

The state object tracks females in two pools:

1. **Eligible pool** (`popFemale_eligible`): Matrix of females past refractory period, can remate
2. **Refractory buffer** (`popFemale_buffer`): 3D array tracking recently-mated females by days since mating

### Daily Cycle

1. Eligible females may remate (redistributed by male availability × preferences)
2. Buffer advances (females age)
3. Females completing refractory period graduate to eligible pool
4. Newly mated females enter buffer at day 1

## Analysis Functions

```r
# Summary statistics
summary <- summarise_multiMating(mm_state)
summary$total_females
summary$n_eligible_to_remate
summary$n_in_refractory
summary$total_by_mate_genotype

# Total mated population (for oviposition calculations)
total_mated <- get_total_mated(mm_state)

# Compare before/after to measure switching
switching <- calc_switching_rate(state_before, state_after)
switching$switch_rate
```

## Notes

- This extension operates externally to MGDrivE; no package modifications required
- Both deterministic and stochastic modes are supported
- The `get_total_mated()` function returns the combined population for use in offspring calculations
