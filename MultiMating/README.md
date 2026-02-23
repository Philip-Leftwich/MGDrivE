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
