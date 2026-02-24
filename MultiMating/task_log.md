# MultiMating Development Log

---

## Session 1 — Initial implementation

**Files created:** `multi_mating.R`, `example_simulation.R`, `README.md`

Built the core multiple mating extension as external code (no modifications to `MGDrivE/R/`).

### Data structures

- `popFemale_eligible` — matrix (nGeno × nGeno): females past the refractory period, eligible to remate
- `popFemale_buffer` — 3D array (nGeno × nGeno × T_refractory): recently-mated females indexed by [female genotype, mate genotype, days since mating]

### Core features implemented

- `init_multiMating` / `init_multiMating_fromPopulation` — initialise state from scratch or from an existing MGDrivE `popFemale` matrix
- `calc_remating_prob` — encounter-rate model: `p_remate = rho × (1 - exp(-lambda × N_males))`
- `remate_deterministic` / `remate_stochastic` — daily remating step
- `advance_buffer` — ages the refractory buffer, graduates females completing their refractory period into the eligible pool
- `add_newly_mated` — adds newly mated females (from standard MGDrivE mating) into the buffer at day 1
- `oneDay_multiMating` — orchestrates the full daily cycle
- `get_total_mated` / `summarise_multiMating` / `calc_switching_rate` — analysis utilities
- Four choosiness modulation functions: `none`, `linear`, `threshold`, `sigmoid`

---

## Session 2 — Conditional remating and mortality correction

**Files modified:** `multi_mating.R`, `example_simulation.R`, `README.md`

### 1. Vector-based conditional remating (`rho` per mate genotype)

**Problem:** `rho` was a scalar applied uniformly to all eligible females, regardless of which male genotype they were currently mated to. This prevented modelling mate-quality-dependent remating (e.g. females mated to preferred males less likely to remate).

**Change:** `rho` is now a numeric vector of length nGeno. Element `rho[j]` is the base daily remating propensity for females currently mated to a male of genotype `j`. A scalar is silently recycled to a uniform vector for backward compatibility.

**Implementation:**
- `remate_deterministic` restructured from a single loop over female genotype `i` to a double loop over `(i, j)`. Each cohort `popFemale_eligible[i, j]` draws its remating probability from `rho[j]`. `mate_probs` (which depends only on female genotype `i`) is now cached outside the inner `j` loop. Removed unused `remating_females` intermediate variable.
- `remate_stochastic` already had a double loop; updated to compute `p_remate_j = calc_remating_prob(popMale, rho[j], lambda)` inside the `j` loop and moved `mate_probs` outside it.

### 2. Mortality drift correction

**Problem:** MGDrivE applies daily adult mortality to its internal `popFemale` each day, but `mm_state` had no corresponding step. Over time, females that died inside MGDrivE remained counted in `mm_state`, causing upward population drift.

**Change:** New function `apply_mortality_multiMating(state, survival_f)`. Accepts a survival probability vector (one per female genotype, `survival_f[i] = 1 - muAd * omega_f[i]`). A scalar is recycled. Applied as **Step 0** at the start of `oneDay_multiMating`, before remating.

**Implementation detail:** For the eligible matrix, R's column-major vector recycling means `mat * vec` applies `vec[i]` to row `i` when `length(vec) == nrow(mat)`. For the 3D buffer, `sweep(..., 1, survival_f, "*")` is used to apply survival along the first dimension (female genotype).

If `survival_f` is `NULL` in params the mortality step is skipped (opt-in).

### Updated daily cycle order

```
0. Apply mortality    (apply_mortality_multiMating)
1. Remate eligible    (remate_deterministic / remate_stochastic)
2. Advance buffer     (advance_buffer)
3. Add newly mated    (add_newly_mated)
```

### Example simulation updates (`example_simulation.R`)

- `muAd` promoted to a shared variable used by both `parameterizeMGDrivE` and `mm_params`
- `rho` changed to `c(AA = 0.02, Aa = 0.05, aa = 0.08)` to demonstrate differential remating
- `survival_f = rep(1 - muAd, 3)` added to `mm_params`

---

## Session 3 — Fix newly_mated extraction

**Files modified:** `example_simulation.R`, `task_log.md`

### Problem

`newly_mated` was computed as the difference in `popFemale` bracketing a full `oneDay_PopDynamics()` call:

```r
popFemale_before <- patch$get_femalePopulation()
patch$oneDay_PopDynamics()
popFemale_after  <- patch$get_femalePopulation()
newly_mated <- popFemale_after - popFemale_before
newly_mated[newly_mated < 0] <- 0
```

`oneDay_PopDynamics()` applies death, maturation, pupation, releases, mating, and oviposition in one call. The difference therefore reflects all of those steps, not mating alone. Clipping negatives discards the mortality signal but the sign of the net change per cell depends on which effect dominates, so the result is unreliable.

### Solution

All individual lifecycle steps are public methods on the `Patch` R6 class. The mating function (`oneDay_mating_deterministic_Patch`) **only adds** to `popFemale` — it takes females from `popUnmated` and distributes them into `popFemale[i, ]` via a weighted outer product. It never subtracts.

Therefore, a snapshot bracketing only the mating call is exact and produces no negatives:

```r
# Run all pre-mating steps individually
patch$oneDay_adultD()
patch$oneDay_pupaDM()
patch$oneDay_larvaDM()
patch$oneDay_eggDM()
patch$oneDay_pupation()
patch$oneDay_releases()

# Exact newly mated females
popFemale_pre_mating <- patch$get_femalePopulation()
patch$oneDay_mating()
newly_mated <- patch$get_femalePopulation() - popFemale_pre_mating

# Complete the day
patch$oneDay_layEggs()
patch$oneDay_releaseEggs()
```

The `newly_mated[newly_mated < 0] <- 0` guard is no longer needed and was removed.

Note: if mated female releases are scheduled (via `matedFemaleReleases`), they are applied in `oneDay_releases()` — before the pre-mating snapshot — so they are correctly excluded from `newly_mated`. If such releases should also be registered in `mm_state`, a separate mechanism would be needed. The current example has no releases.

---

---

## Session 4 — Code audit and fixes

**Files modified:** `multi_mating.R`, `task_log.md`

Four issues identified by cross-checking `task_log.md` against the R source files.

### Fix 1 — `example_usage` docstring out of date (Issue 1, must fix)

The `\dontrun{}` example in the `example_usage` function had not been updated to reflect either Session 2 (vector `rho`, `survival_f`) or Session 3 (correct `newly_mated` extraction via individual lifecycle steps). The old pattern used `popFemale_after - popFemale_before` with negative clipping. Replaced with the correct pattern matching `example_simulation.R`.

### Fix 2 — `rmultinom` matrix coercion in `remate_stochastic` (Issue 2)

`rmultinom(n = 1, ...)` returns an `nGeno × 1` matrix. Adding this to `state$popFemale_buffer[i, , 1]` (a vector) produced a matrix via R's implicit coercion, which was then assigned back to the array slice relying on silent element-count matching. Wrapped with `drop()` to explicitly convert the result to a plain named vector before the addition.

### Fix 3 — Indentation inconsistencies (Issue 3)

Two comment lines were not indented to match their surrounding code block:
- Line 44 (`# Buffer for recently mated females`) inside `init_multiMating`
- Line 616 (`# Total by mate genotype`) inside `summarise_multiMating`

### Fix 4 — `calc_switching_rate` precondition documented (Issue 4)

Added a docstring note that the zero-sum assumption underlying `sum(abs(delta)) / 2` only holds when the two states differ solely due to remating. Calling the function around a full `oneDay_multiMating` call (which also applies mortality) violates this assumption.

---

## Known issues / outstanding items

- **`omega_f` handling.** If genotype-specific female mortality fitness (`omega_f`) differs across genotypes, `survival_f` should be computed as `1 - muAd * omega_f` (vector). The example uses a uniform scalar which is only exact when `omega_f = rep(1, nGeno)`.
