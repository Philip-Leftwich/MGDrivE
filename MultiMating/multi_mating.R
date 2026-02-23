###############################################################################
#
#   Multiple Mating Extension for MGDrivE
#
#   This script provides functions to extend MGDrivE with multiple mating
#   capabilities. It operates externally to the package, wrapping the
#   existing simulation framework.
#
#   Features:
#     - Refractory period before females become eligible to remate
#     - Encounter-rate dependent remating probability
#     - Dynamic choosiness modulation based on male genotype frequencies
#     - Last-male sperm precedence
#
###############################################################################

#' Initialize Multiple Mating Data Structures
#'
#' Creates the data structures needed to track female mating history,
#' including the refractory buffer and eligible pool.
#'
#' @param genotypesN Integer, number of genotypes
#' @param genotypesID Character vector of genotype names
#' @param T_refractory Integer, days after mating before eligible to remate
#'
#' @return A list containing:
#'   - popFemale_eligible: matrix of females eligible to remate
#'   - popFemale_buffer: 3D array of recently-mated females (in refractory period)
#'   - T_refractory: the refractory period
#'   - genotypesN: number of genotypes
#'   - genotypesID: genotype names
#'
init_multiMating <- function(genotypesN, genotypesID, T_refractory) {


  # Eligible females (mating age >= T_refractory)
  popFemale_eligible <- matrix(
    data = 0,
    nrow = genotypesN,
    ncol = genotypesN,
    dimnames = list(genotypesID, genotypesID)
  )

# Buffer for recently mated females (days 1 to T_refractory)
  # Third dimension is days since mating: 1 = just mated today, T_refractory = about to become eligible
  popFemale_buffer <- array(
    data = 0,
    dim = c(genotypesN, genotypesN, T_refractory),
    dimnames = list(genotypesID, genotypesID, paste0("day_", 1:T_refractory))
  )

  return(list(
    popFemale_eligible = popFemale_eligible,
    popFemale_buffer = popFemale_buffer,
    T_refractory = T_refractory,
    genotypesN = genotypesN,
    genotypesID = genotypesID
  ))
}


#' Initialize Multiple Mating State from Existing Population
#'
#' Converts an existing popFemale matrix (from standard MGDrivE) into
#' the multiple mating data structure. Assumes all existing females
#' are past the refractory period (eligible to remate).
#'
#' @param popFemale Matrix (nGeno x nGeno) of mated females from MGDrivE
#' @param T_refractory Integer, days after mating before eligible to remate
#'
#' @return A multiMating state list (see init_multiMating)
#'
init_multiMating_fromPopulation <- function(popFemale, T_refractory) {

  genotypesN <- nrow(popFemale)
  genotypesID <- rownames(popFemale)

  # Initialize empty structures
  state <- init_multiMating(genotypesN, genotypesID, T_refractory)

  # Place all existing females in the eligible pool
  state$popFemale_eligible[] <- popFemale

  return(state)
}


###############################################################################
# Choosiness Modulation Functions
###############################################################################

#' No Choosiness Modulation (Default)
#'
#' Returns base eta unchanged. Use this when you don't want dynamic preferences.
#'
#' @param base_eta Matrix (nGeno x nGeno), base mating preferences
#' @param popMale Vector of male population by genotype
#' @param params List of additional parameters (unused)
#'
#' @return Matrix of effective eta (unchanged from base)
#'
modulate_choosiness_none <- function(base_eta, popMale, params = NULL) {
  return(base_eta)
}


#' Linear Relaxation of Choosiness
#'
#' As preferred males become rare, females relax their preferences toward
#' uniform (non-selective) mating. The degree of relaxation is proportional
#' to how rare preferred males are.
#'
#' @param base_eta Matrix (nGeno x nGeno), base mating preferences
#' @param popMale Vector of male population by genotype
#' @param params List containing:
#'   - strictness: threshold for "preferred males available" (0-1)
#'                 Below this frequency, choosiness starts to relax
#'
#' @return Matrix of effective eta (relaxed toward uniform)
#'
modulate_choosiness_linear <- function(base_eta, popMale, params) {

  strictness <- params$strictness
  if (is.null(strictness)) strictness <- 0.2

  # Calculate male frequencies
  M_total <- sum(popMale)
  if (M_total == 0) return(base_eta)

  male_freq <- popMale / M_total

  # For each female genotype, calculate relaxation
  effective_eta <- base_eta

  for (i in seq_len(nrow(base_eta))) {
    # Identify "preferred" male genotypes (above-average preference)
    preferred <- base_eta[i, ] > mean(base_eta[i, ])
    preferred_available <- sum(male_freq[preferred])

    # Relaxation factor: 0 when preferred males common, 1 when very rare
    relaxation <- max(0, 1 - preferred_available / strictness)
    relaxation <- min(1, relaxation)

    # Interpolate between base preferences and uniform
    uniform <- rep(1, ncol(base_eta))
    effective_eta[i, ] <- base_eta[i, ] * (1 - relaxation) + uniform * relaxation
  }

  return(effective_eta)
}


#' Threshold-Based Choosiness Switching
#'
#' Females use strict preferences if any preferred male genotype exceeds
#' a threshold frequency. Otherwise, they become non-selective.
#'
#' @param base_eta Matrix (nGeno x nGeno), base mating preferences
#' @param popMale Vector of male population by genotype
#' @param params List containing:
#'   - threshold: minimum frequency of preferred males to maintain choosiness
#'
#' @return Matrix of effective eta (either base or uniform)
#'
modulate_choosiness_threshold <- function(base_eta, popMale, params) {

  threshold <- params$threshold
  if (is.null(threshold)) threshold <- 0.1

  # Calculate male frequencies
  M_total <- sum(popMale)
  if (M_total == 0) return(base_eta)

  male_freq <- popMale / M_total

  effective_eta <- base_eta

  for (i in seq_len(nrow(base_eta))) {
    # Identify "preferred" male genotypes (preference > 1, i.e., above neutral)
    preferred <- base_eta[i, ] > 1

    # Check if any preferred genotype exceeds threshold
    if (any(preferred) && any(male_freq[preferred] > threshold)) {
      # Maintain strict preferences
      effective_eta[i, ] <- base_eta[i, ]
    } else {
      # Become non-selective
      effective_eta[i, ] <- rep(1, ncol(base_eta))
    }
  }

  return(effective_eta)
}


#' Sigmoid Choosiness Relaxation
#'
#' Choosiness relaxes sigmoidally as preferred males become rare.
#' Provides smooth transition with adjustable steepness.
#'
#' @param base_eta Matrix (nGeno x nGeno), base mating preferences
#' @param popMale Vector of male population by genotype
#' @param params List containing:
#'   - midpoint: frequency at which relaxation is 50%
#'   - steepness: how sharp the transition is
#'
#' @return Matrix of effective eta
#'
modulate_choosiness_sigmoid <- function(base_eta, popMale, params) {

  midpoint <- params$midpoint
  steepness <- params$steepness
  if (is.null(midpoint)) midpoint <- 0.15
  if (is.null(steepness)) steepness <- 20

  # Calculate male frequencies
  M_total <- sum(popMale)
  if (M_total == 0) return(base_eta)

  male_freq <- popMale / M_total

  effective_eta <- base_eta

  for (i in seq_len(nrow(base_eta))) {
    # Identify "preferred" male genotypes
    preferred <- base_eta[i, ] > mean(base_eta[i, ])
    preferred_available <- sum(male_freq[preferred])

    # Sigmoid relaxation: high when preferred rare, low when common
    # relaxation = 1 / (1 + exp(steepness * (preferred_available - midpoint)))
    relaxation <- 1 / (1 + exp(steepness * (preferred_available - midpoint)))

    # Interpolate
    uniform <- rep(1, ncol(base_eta))
    effective_eta[i, ] <- base_eta[i, ] * (1 - relaxation) + uniform * relaxation
  }

  return(effective_eta)
}


###############################################################################
# Core Remating Functions
###############################################################################

#' Normalise a Vector to Sum to One
#'
#' Utility function matching MGDrivE's internal normalise.
#'
#' @param vector Numeric vector
#'
#' @return Normalised vector (sums to 1), or original if all zeros
#'
normalise <- function(vector) {
  if (all(vector == 0)) {
    return(vector)
  } else {
    return(vector / sum(vector))
  }
}


#' Calculate Remating Probability
#'
#' Computes the daily probability that an eligible female remates,
#' based on encounter rate and base propensity.
#'
#' @param popMale Vector of male population by genotype
#' @param rho Base remating propensity (0-1)
#' @param lambda Encounter rate parameter
#'
#' @return Scalar probability of remating
#'
calc_remating_prob <- function(popMale, rho, lambda) {
  M_total <- sum(popMale)

  # Encounter probability increases with male density
  p_encounter <- 1 - exp(-lambda * M_total)

  # Overall remating probability
  p_remate <- rho * p_encounter

  return(p_remate)
}


#' Perform Remating Step (Deterministic)
#'
#' Executes one day of remating for eligible females.
#' Females are extracted proportionally and redistributed according to
#' current male availability and (optionally modulated) preferences.
#'
#' @param state MultiMating state list
#' @param popMale Vector of male population by genotype
#' @param base_eta Matrix of base mating preferences (from driveCube)
#' @param rho Base remating propensity
#' @param lambda Encounter rate parameter
#' @param modulate_fn Choosiness modulation function
#' @param modulate_params Parameters for choosiness modulation
#'
#' @return Updated state list
#'
remate_deterministic <- function(state, popMale, base_eta, rho, lambda,
                                  modulate_fn = modulate_choosiness_none,
                                  modulate_params = NULL) {

  # Calculate remating probability
  p_remate <- calc_remating_prob(popMale, rho, lambda)

  # Skip if no remating occurs
  if (p_remate == 0 || sum(state$popFemale_eligible) == 0) {
    return(state)
  }

  # Check if males are available
  if (sum(popMale) == 0) {
    return(state)
  }

  # Get effective preferences (with choosiness modulation)
  effective_eta <- modulate_fn(base_eta, popMale, modulate_params)

  nGeno <- state$genotypesN

  # For each female genotype
  for (i in seq_len(nGeno)) {

    # Total eligible females of this genotype (across all mate genotypes)
    total_eligible_i <- sum(state$popFemale_eligible[i, ])

    if (total_eligible_i == 0) next

    # Number remating
    n_remating <- total_eligible_i * p_remate

    # Remove from current positions (proportionally from each mate column)
    removal_fraction <- p_remate
    remating_females <- state$popFemale_eligible[i, ] * removal_fraction
    state$popFemale_eligible[i, ] <- state$popFemale_eligible[i, ] * (1 - removal_fraction)

    # Redistribute according to male availability and preferences
    mate_probs <- normalise(popMale * effective_eta[i, ])

    # Add to new positions (enters buffer at day 1)
    state$popFemale_buffer[i, , 1] <- state$popFemale_buffer[i, , 1] +
                                       n_remating * mate_probs
  }

  return(state)
}


#' Perform Remating Step (Stochastic)
#'
#' Stochastic version of remating. Uses binomial sampling for number
#' remating and multinomial for mate choice.
#'
#' @param state MultiMating state list
#' @param popMale Vector of male population by genotype
#' @param base_eta Matrix of base mating preferences (from driveCube)
#' @param rho Base remating propensity
#' @param lambda Encounter rate parameter
#' @param modulate_fn Choosiness modulation function
#' @param modulate_params Parameters for choosiness modulation
#'
#' @return Updated state list
#'
remate_stochastic <- function(state, popMale, base_eta, rho, lambda,
                               modulate_fn = modulate_choosiness_none,
                               modulate_params = NULL) {

  # Calculate remating probability
  p_remate <- calc_remating_prob(popMale, rho, lambda)

  # Skip if no remating occurs
  if (p_remate == 0 || sum(state$popFemale_eligible) == 0) {
    return(state)
  }

  # Check if males are available
  if (sum(popMale) == 0) {
    return(state)
  }

  # Get effective preferences (with choosiness modulation)
  effective_eta <- modulate_fn(base_eta, popMale, modulate_params)

  nGeno <- state$genotypesN

  # For each female genotype
  for (i in seq_len(nGeno)) {

    # For each current mate genotype
    for (j in seq_len(nGeno)) {

      n_eligible <- state$popFemale_eligible[i, j]

      if (n_eligible == 0) next

      # Sample number remating
      n_remating <- rbinom(n = 1, size = round(n_eligible), prob = p_remate)

      if (n_remating == 0) next

      # Remove from current position
      state$popFemale_eligible[i, j] <- state$popFemale_eligible[i, j] - n_remating

      # Choose new mates
      mate_probs <- normalise(popMale * effective_eta[i, ])

      if (sum(mate_probs) > 0) {
        new_mates <- rmultinom(n = 1, size = n_remating, prob = mate_probs)
        # Add to buffer at day 1
        state$popFemale_buffer[i, , 1] <- state$popFemale_buffer[i, , 1] + new_mates
      }
    }
  }

  return(state)
}


#' Advance Buffer by One Day
#'
#' Shifts the refractory buffer forward by one day:
#' - Females at day T_refractory graduate to eligible pool
#' - All other cohorts age by one day
#' - Day 1 position is cleared (will receive newly mated females)
#'
#' @param state MultiMating state list
#'
#' @return Updated state list
#'
advance_buffer <- function(state) {

  T_ref <- state$T_refractory

  # Graduate females completing refractory period
  state$popFemale_eligible <- state$popFemale_eligible +
                               state$popFemale_buffer[, , T_ref]

  # Shift buffer: day t -> day t+1 (working backwards to avoid overwriting)
  if (T_ref > 1) {
    for (t in T_ref:2) {
      state$popFemale_buffer[, , t] <- state$popFemale_buffer[, , t - 1]
    }
  }

  # Clear day 1 (will receive newly mated females)
  state$popFemale_buffer[, , 1] <- 0

  return(state)
}


#' Add Newly Mated Females to Buffer
#'
#' Takes newly mated females (from standard MGDrivE mating) and adds
#' them to the refractory buffer at day 1.
#'
#' @param state MultiMating state list
#' @param newly_mated Matrix (nGeno x nGeno) of newly mated females
#'
#' @return Updated state list
#'
add_newly_mated <- function(state, newly_mated) {
  state$popFemale_buffer[, , 1] <- state$popFemale_buffer[, , 1] + newly_mated
  return(state)
}


#' Get Total Mated Females (for Oviposition)
#'
#' Returns the total mated female population across both the eligible
#' pool and all refractory buffer slots. This is what should be used
#' for oviposition calculations.
#'
#' @param state MultiMating state list
#'
#' @return Matrix (nGeno x nGeno) of total mated females
#'
get_total_mated <- function(state) {
  # Sum across all buffer days
  buffer_total <- apply(state$popFemale_buffer, c(1, 2), sum)

  # Add eligible pool
  total <- state$popFemale_eligible + buffer_total

  return(total)
}


###############################################################################
# Daily Cycle Integration
###############################################################################

#' Run One Day of Multiple Mating Dynamics
#'
#' Executes a complete daily cycle for the multiple mating extension:
#' 1. Remating of eligible females
#' 2. Advance refractory buffer (aging)
#' 3. Add newly mated females from standard mating
#'
#' This should be called AFTER the standard MGDrivE daily dynamics,
#' with the newly_mated females extracted from popFemale.
#'
#' @param state MultiMating state list
#' @param popMale Vector of male population by genotype
#' @param base_eta Matrix of base mating preferences (from driveCube)
#' @param newly_mated Matrix of newly mated females from this day's mating
#' @param params List of parameters:
#'   - rho: base remating propensity
#'   - lambda: encounter rate parameter
#'   - stochastic: logical, use stochastic (TRUE) or deterministic (FALSE)
#'   - modulate_fn: choosiness modulation function (optional)
#'   - modulate_params: parameters for choosiness modulation (optional)
#'
#' @return Updated state list
#'
oneDay_multiMating <- function(state, popMale, base_eta, newly_mated, params) {

  # Extract parameters
  rho <- params$rho
  lambda <- params$lambda
  stochastic <- params$stochastic
  modulate_fn <- params$modulate_fn
  modulate_params <- params$modulate_params

  # Defaults
  if (is.null(stochastic)) stochastic <- FALSE
  if (is.null(modulate_fn)) modulate_fn <- modulate_choosiness_none
  if (is.null(modulate_params)) modulate_params <- list()

  # Step 1: Remating of eligible females
  if (stochastic) {
    state <- remate_stochastic(state, popMale, base_eta, rho, lambda,
                                modulate_fn, modulate_params)
  } else {
    state <- remate_deterministic(state, popMale, base_eta, rho, lambda,
                                   modulate_fn, modulate_params)
  }

  # Step 2: Advance buffer (aging)
  state <- advance_buffer(state)

  # Step 3: Add newly mated females to buffer
  state <- add_newly_mated(state, newly_mated)

  return(state)
}


###############################################################################
# Analysis and Reporting Functions
###############################################################################

#' Summarise Multiple Mating State
#'
#' Provides a summary of the current multiple mating population structure.
#'
#' @param state MultiMating state list
#'
#' @return List with summary statistics
#'
summarise_multiMating <- function(state) {

  total_mated <- get_total_mated(state)

  # Total by female genotype
  by_female <- rowSums(total_mated)

# Total by mate genotype
  by_mate <- colSums(total_mated)

  # Totals in each pool
  n_eligible <- sum(state$popFemale_eligible)
  n_refractory <- sum(state$popFemale_buffer)

  # Age distribution in buffer
  age_dist <- apply(state$popFemale_buffer, 3, sum)

  return(list(
    total_females = sum(total_mated),
    total_by_female_genotype = by_female,
    total_by_mate_genotype = by_mate,
    n_eligible_to_remate = n_eligible,
    n_in_refractory = n_refractory,
    refractory_age_distribution = age_dist,
    popFemale_combined = total_mated
  ))
}


#' Calculate Mate Switching Rate
#'
#' Compares mating distribution before and after remating to quantify
#' how much mate-switching occurred.
#'
#' @param state_before MultiMating state before remating
#' @param state_after MultiMating state after remating
#'
#' @return List with switching statistics
#'
calc_switching_rate <- function(state_before, state_after) {

  before <- get_total_mated(state_before)
  after <- get_total_mated(state_after)

  # Net change by mate genotype (for each female genotype)
  delta <- after - before

  # Total females that moved (half of absolute changes, since it's zero-sum)
  n_switched <- sum(abs(delta)) / 2

  # Fraction of population that switched
  total_pop <- sum(before)
  switch_rate <- if (total_pop > 0) n_switched / total_pop else 0

  return(list(
    n_switched = n_switched,
    switch_rate = switch_rate,
    delta_by_cross = delta
  ))
}


###############################################################################
# Example Usage / Integration Template
###############################################################################

#' Example: Run Simulation with Multiple Mating
#'
#' This is a template showing how to integrate multiple mating with
#' an MGDrivE simulation. NOT a complete working example - requires
#' an actual MGDrivE network setup.
#'
#' @examples
#' \dontrun{
#' # After setting up MGDrivE network...
#'
#' # Initialize multiple mating state
#' mm_params <- list(
#'   T_refractory = 3,        # 3 day refractory period
#'   rho = 0.1,               # 10% base remating propensity
#'   lambda = 0.01,           # encounter rate
#'   stochastic = TRUE,
#'   modulate_fn = modulate_choosiness_linear,
#'   modulate_params = list(strictness = 0.2)
#' )
#'
#' # Get initial population and create state
#' popFemale_init <- patch$get_femalePopulation()
#' mm_state <- init_multiMating_fromPopulation(popFemale_init, mm_params$T_refractory)
#'
#' # Get base eta from drive cube
#' base_eta <- driveCube$eta
#'
#' # In daily simulation loop:
#' for (day in 1:simTime) {
#'
#'   # Record female population before mating
#'   popFemale_before <- patch$get_femalePopulation()
#'
#'   # Run standard MGDrivE daily dynamics
#'   # ... (death, maturation, pupation, releases, mating, oviposition)
#'
#'   # Get population after mating
#'   popFemale_after <- patch$get_femalePopulation()
#'
#'   # Extract newly mated females (difference)
#'   newly_mated <- popFemale_after - popFemale_before
#'   newly_mated[newly_mated < 0] <- 0  # only additions
#'
#'   # Get current male population
#'   popMale <- patch$get_malePopulation()
#'
#'   # Run multiple mating dynamics
#'   mm_state <- oneDay_multiMating(
#'     state = mm_state,
#'     popMale = popMale,
#'     base_eta = base_eta,
#'     newly_mated = newly_mated,
#'     params = mm_params
#'   )
#'
#'   # The total mated population for analysis
#'   total_mated <- get_total_mated(mm_state)
#' }
#' }
#'
example_usage <- function() {
  cat("See the function documentation for usage template.\n")
  cat("This requires an actual MGDrivE network setup.\n")
}
