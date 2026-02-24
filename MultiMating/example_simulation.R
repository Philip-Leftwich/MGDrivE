###############################################################################
#
#   Example: MGDrivE Simulation with Multiple Mating
#
#   This script demonstrates how to run an MGDrivE simulation with the
#   multiple mating extension. It sets up a simple single-patch scenario
#   with Mendelian inheritance to illustrate the integration.
#
###############################################################################

# Load required packages
library(MGDrivE)

# Source the multiple mating extension
source("multi_mating.R")


###############################################################################
# Simulation Parameters
###############################################################################

# Standard MGDrivE parameters
simTime <- 365        # simulation days
sampTime <- 1         # output every day
nPatch <- 1           # single patch
adultPopEQ <- 500     # equilibrium adult population
muAd <- 0.09          # adult daily mortality rate (shared with mm_params)

# Multiple mating parameters
# rho is now a named vector: one value per male genotype (ordered to match cube$genotypesID).
# Females mated to preferred males (AA) have lower remating propensity than
# those mated to less-preferred males (aa). A scalar would apply uniformly.
mm_params <- list(
  T_refractory = 5,                          # 5 day refractory period after mating
  rho = c(AA = 0.02, Aa = 0.05, aa = 0.08), # genotype-conditional remating propensity
  lambda = 0.001,                            # encounter rate parameter
  stochastic = FALSE,                        # deterministic for this example
  survival_f = rep(1 - muAd, 3),            # daily survival per female genotype (uniform here)
  modulate_fn = modulate_choosiness_linear,
  modulate_params = list(strictness = 0.2)
)


###############################################################################
# Setup MGDrivE
###############################################################################

# Initialize for deterministic simulation
setupMGDrivE(stochasticityON = FALSE, verbose = TRUE)

# Create inheritance cube (simple Mendelian with 3 genotypes)
# AA = wild-type, Aa = heterozygote, aa = homozygote
cube <- cubeMendelian(
  gtype = c("AA", "Aa", "aa"),
  eta = NULL,   # default: all matings equally likely
  phi = NULL,   # default: 50% female
  omega = NULL, # default: no mortality differences
  xiF = NULL,   # default: equal female emergence
  xiM = NULL,   # default: equal male emergence
  s = NULL      # default: equal fertility
)

# Alternatively, set up preferences where females prefer AA males:
# cube <- cubeMendelian(
#   gtype = c("AA", "Aa", "aa"),
#   eta = list(
#     c("AA", 1.5),   # AA males 50% more attractive
#     c("aa", 0.5)    # aa males 50% less attractive
#   )
# )

# Generate parameters
params <- parameterizeMGDrivE(
  runID = 1,
  simTime = simTime,
  sampTime = sampTime,
  nPatch = nPatch,
  beta = 32,              # eggs per female per day
  muAd = muAd,            # adult mortality rate
  popGrowth = 1.1,        # population growth rate
  tEgg = 2,               # days in egg stage
  tLarva = 8,             # days in larva stage
  tPupa = 1,              # days in pupa stage
  AdssPopEQ = rep(adultPopEQ, nPatch),
  inheritanceCube = cube
)

# Create release schedule (no releases in this example)
patchReleases <- replicate(
  n = nPatch,
  expr = list(
    maleReleases = NULL,
    femaleReleases = NULL,
    eggReleases = NULL,
    matedFemaleReleases = NULL
  ),
  simplify = FALSE
)

# Migration matrices (single patch = no migration)
migrationMale <- matrix(1, nrow = nPatch, ncol = nPatch)
migrationFemale <- matrix(1, nrow = nPatch, ncol = nPatch)

# Output directory
outDir <- tempdir()

# Create network
network <- Network$new(
  params = params,
  driveCube = cube,
  patchReleases = patchReleases,
  migrationMale = migrationMale,
  migrationFemale = migrationFemale,
  directory = outDir,
  verbose = TRUE
)


###############################################################################
# Initialize Multiple Mating State
###############################################################################

# Get reference to patch
patch <- network$.__enclos_env__$private$patches[[1]]

# Get initial female population
popFemale_init <- patch$get_femalePopulation()

# Initialize multiple mating state
mm_state <- init_multiMating_fromPopulation(popFemale_init, mm_params$T_refractory)

# Get base eta from drive cube
base_eta <- cube$eta

cat("\n=== Initial State ===\n")
cat("Female population matrix:\n")
print(round(popFemale_init, 2))
cat("\nMultiple mating state initialized.\n")
cat("Refractory period:", mm_params$T_refractory, "days\n")
cat("Eligible to remate:", sum(mm_state$popFemale_eligible), "\n")


###############################################################################
# Run Simulation with Multiple Mating
###############################################################################

cat("\n=== Running Simulation ===\n")

# Storage for results
results <- data.frame(
  day = integer(),
  total_females = numeric(),
  eligible_to_remate = numeric(),
  in_refractory = numeric(),
  mated_to_AA = numeric(),
  mated_to_Aa = numeric(),
  mated_to_aa = numeric()
)

# Run simulation
for (day in 1:simTime) {

  # Record female population before mating step
  popFemale_before <- patch$get_femalePopulation()

  # Run standard MGDrivE daily dynamics
  # Note: In a full integration, you'd call the network's oneDay() method
  # Here we manually step through to demonstrate the integration point
  patch$oneDay_PopDynamics()

  # Get population after mating
  popFemale_after <- patch$get_femalePopulation()

  # Extract newly mated females (females that weren't mated before)
  # This is the increase in popFemale from the mating step
  newly_mated <- popFemale_after - popFemale_before
  newly_mated[newly_mated < 0] <- 0  # only count additions

  # Get current male population
  popMale <- patch$get_malePopulation()

  # Run multiple mating dynamics
  mm_state <- oneDay_multiMating(
    state = mm_state,
    popMale = popMale,
    base_eta = base_eta,
    newly_mated = newly_mated,
    params = mm_params
  )

  # Get summary for this day
  summary <- summarise_multiMating(mm_state)

  # Store results
  results <- rbind(results, data.frame(
    day = day,
    total_females = summary$total_females,
    eligible_to_remate = summary$n_eligible_to_remate,
    in_refractory = summary$n_in_refractory,
    mated_to_AA = summary$total_by_mate_genotype["AA"],
    mated_to_Aa = summary$total_by_mate_genotype["Aa"],
    mated_to_aa = summary$total_by_mate_genotype["aa"]
  ))

  # Progress update
  if (day %% 50 == 0) {
    cat("Day", day, "- Total mated females:", round(summary$total_females), "\n")
  }
}

cat("\n=== Simulation Complete ===\n")


###############################################################################
# Analysis and Plotting
###############################################################################

cat("\n=== Results Summary ===\n")
cat("Final total mated females:", round(tail(results$total_females, 1)), "\n")
cat("Final females eligible to remate:", round(tail(results$eligible_to_remate, 1)), "\n")
cat("Final females in refractory:", round(tail(results$in_refractory, 1)), "\n")

# Final mating distribution
final_state <- summarise_multiMating(mm_state)
cat("\nFinal mating distribution (by mate genotype):\n")
print(round(final_state$total_by_mate_genotype, 2))

cat("\nFinal combined popFemale matrix:\n")
print(round(final_state$popFemale_combined, 2))

# Plot results if running interactively
if (interactive()) {

  par(mfrow = c(2, 2))

  # Plot 1: Total females over time
  plot(results$day, results$total_females,
       type = "l", col = "black", lwd = 2,
       xlab = "Day", ylab = "Count",
       main = "Total Mated Females")

  # Plot 2: Eligible vs refractory
  plot(results$day, results$eligible_to_remate,
       type = "l", col = "blue", lwd = 2,
       xlab = "Day", ylab = "Count",
       main = "Eligible vs Refractory",
       ylim = c(0, max(results$eligible_to_remate, results$in_refractory)))
  lines(results$day, results$in_refractory, col = "red", lwd = 2)
  legend("right", legend = c("Eligible", "Refractory"),
         col = c("blue", "red"), lwd = 2)

  # Plot 3: Mate genotype distribution
  plot(results$day, results$mated_to_AA,
       type = "l", col = "forestgreen", lwd = 2,
       xlab = "Day", ylab = "Count",
       main = "Mate Genotype Distribution",
       ylim = c(0, max(c(results$mated_to_AA, results$mated_to_Aa, results$mated_to_aa))))
  lines(results$day, results$mated_to_Aa, col = "orange", lwd = 2)
  lines(results$day, results$mated_to_aa, col = "purple", lwd = 2)
  legend("right", legend = c("AA", "Aa", "aa"),
         col = c("forestgreen", "orange", "purple"), lwd = 2)

  # Plot 4: Proportion by mate genotype
  total <- results$mated_to_AA + results$mated_to_Aa + results$mated_to_aa
  plot(results$day, results$mated_to_AA / total,
       type = "l", col = "forestgreen", lwd = 2,
       xlab = "Day", ylab = "Proportion",
       main = "Mate Genotype Proportions",
       ylim = c(0, 1))
  lines(results$day, results$mated_to_Aa / total, col = "orange", lwd = 2)
  lines(results$day, results$mated_to_aa / total, col = "purple", lwd = 2)
  legend("right", legend = c("AA", "Aa", "aa"),
         col = c("forestgreen", "orange", "purple"), lwd = 2)

  par(mfrow = c(1, 1))
}


###############################################################################
# Parameter Sensitivity Example
###############################################################################

cat("\n=== Parameter Sensitivity Notes ===\n")
cat("
Key parameters and their effects:

T_refractory (refractory period):
  - Higher values: fewer females eligible to remate at any time
  - Lower values: faster population turnover in mate distribution

rho (base remating propensity):
  - Higher values: more remating events per day
  - Lower values: slower mate switching

lambda (encounter rate):
  - Higher values: remating approaches rho even at low male density

  - Lower values: remating strongly dependent on male density

Choosiness modulation:
  - modulate_choosiness_none: static preferences
  - modulate_choosiness_linear: gradual relaxation as preferred males become rare
  - modulate_choosiness_threshold: binary switch at frequency threshold
  - modulate_choosiness_sigmoid: smooth S-curve transition

To explore different scenarios, modify mm_params and re-run.
")
