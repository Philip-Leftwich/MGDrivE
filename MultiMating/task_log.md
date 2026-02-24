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

## Known issues / outstanding items

- **`newly_mated` extraction is incorrect.** The example computes `newly_mated` as `popFemale_after - popFemale_before` where both snapshots bracket a full `oneDay_PopDynamics()` call. That call applies death, maturation, pupation, releases, mating, and oviposition — not just mating. The difference reflects all of these, not mating alone. Approach to fix: TBD (next session).
- **`omega_f` handling.** If genotype-specific female mortality fitness (`omega_f`) differs across genotypes, `survival_f` should be computed as `1 - muAd * omega_f` (vector). The example uses a uniform scalar which is only exact when `omega_f = rep(1, nGeno)`.
