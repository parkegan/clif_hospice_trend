#!/usr/bin/env Rscript
# =============================================================================
# aim3_chart_validation.R  —  Aim 3 chart-validation substudy (protocol v7.0 §3.6.3)
# -----------------------------------------------------------------------------
# Aim 3 is EXPLORATORY and contingent on validating two algorithm-inferred
# constructs against the medical record:
#   (1) the same-facility hospice-transition algorithm, and
#   (2) the AND -> comfort-care mapping.
# The protocol asks for 50-100 encounters reviewed at 2-3 sites, and reports the
# positive predictive value (PPV) of each construct.
#
# The chart review itself is MANUAL. This script does the two things code can:
#   MODE 1 (default) — draw a stratified random sample of algorithm-flagged
#     encounters and write a BLANK review worksheet for abstraction.
#   MODE 2 (env AIM3_CV_SCORE=/path/to/completed_worksheet.csv) — read a
#     completed worksheet and compute PPV with Wilson confidence intervals.
#
# Nothing here fabricates a gold standard: PPV is computed only from reviewer-
# entered columns, and MODE 2 refuses to run on a worksheet whose review columns
# are still blank.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
if (!exists("data_list")) source(here::here("code", "00_setup.R"))

set.seed(20260713)                      # reproducible sample
N_REVIEW      = 80L                     # target sample size (protocol: 50-100)
score_path    = Sys.getenv("AIM3_CV_SCORE", "")

# Wilson score interval for a proportion (better than normal at small n / extreme p)
wilson_ci = function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z = qnorm(1 - (1 - conf) / 2)
  phat = x / n
  denom = 1 + z^2 / n
  centre = (phat + z^2 / (2 * n)) / denom
  half = z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  round(c(centre - half, centre + half), 4)
}


# =============================================================================
# MODE 2: score a completed worksheet
# =============================================================================
if (nzchar(score_path)) {

  stopifnot("worksheet not found" = file.exists(score_path))
  ws = fread(score_path, colClasses = list(character = "hospitalization_id"))

  need = c("hospitalization_id", "algo_hospice_transition", "algo_comfort_and",
           "chart_hospice_confirmed", "chart_comfort_confirmed")
  miss = setdiff(need, names(ws))
  if (length(miss))
    stop("Worksheet missing required columns: ", paste(miss, collapse = ", "))

  # reviewer columns must be filled (0/1); refuse to score a blank worksheet
  norm01 = function(x) {
    v = tolower(trimws(as.character(x)))
    fcase(v %in% c("1", "y", "yes", "true"), 1L,
          v %in% c("0", "n", "no", "false"), 0L,
          default = NA_integer_)
  }
  ws[, chart_hospice_confirmed := norm01(chart_hospice_confirmed)]
  ws[, chart_comfort_confirmed := norm01(chart_comfort_confirmed)]

  # PPV of the hospice-transition algorithm: among algorithm-flagged transitions,
  # the share the chart confirms.
  hp = ws[algo_hospice_transition == 1L & !is.na(chart_hospice_confirmed)]
  # PPV of AND -> comfort: among AND-flagged encounters, the share with
  # chart-documented comfort-focused care.
  cf = ws[algo_comfort_and == 1L & !is.na(chart_comfort_confirmed)]

  if (nrow(hp) == 0 && nrow(cf) == 0)
    stop("No completed review rows found — fill chart_* columns before scoring.")

  ppv = function(dt, col, construct) {
    if (nrow(dt) == 0) {
      message(sprintf("  %s: no reviewed rows.", construct)); return(NULL)
    }
    x = sum(dt[[col]] == 1L); n = nrow(dt)
    ci = wilson_ci(x, n)
    message(sprintf("  %-38s PPV = %d/%d = %.1f%%  (95%% CI %.1f-%.1f)",
                    construct, x, n, 100 * x / n, 100 * ci[1], 100 * ci[2]))
    data.table(construct = construct, n_reviewed = n, n_confirmed = x,
               ppv = round(x / n, 4), ci_lower = ci[1], ci_upper = ci[2],
               site = site_lowercase)
  }

  message("\n===== AIM 3 CHART-VALIDATION PPV =====")
  out = rbindlist(list(
    ppv(hp, "chart_hospice_confirmed", "Same-facility hospice-transition algorithm"),
    ppv(cf, "chart_comfort_confirmed", "AND -> comfort-care mapping")
  ), fill = TRUE)
  message("======================================\n")

  fwrite(out, shareable_path(
    paste0("aim3_chart_validation_ppv_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim3_chart_validation_ppv.")
  quit(save = "no")
}


# =============================================================================
# MODE 1: draw the review sample and write a blank worksheet
# =============================================================================

# algorithm outputs: the cohort (hospice-transition determinations) joined to the
# inflection flags (comfort/AND detection).
cohort = fread(latest_output("cohort_all_hospitalizations"),
               colClasses = list(character = c("patient_id", "hospitalization_id")))
flags  = fread(latest_output("clinical_inflection_flags"),
               colClasses = list(character = "hospitalization_id"))

# candidate pool: encounters the algorithm FLAGGED as same-facility hospice
# transitions (these are what a PPV is defined over).
sf = cohort[outcome_category %in% c("same_facility_hospice_death",
                                    "same_facility_hospice_alive")]
if (nrow(sf) == 0) stop("No same-facility hospice transitions to validate at this site.")

# attach the comfort/AND determination if available
comfort_col = intersect(c("first_comfort_dttm", "flag_comfort", "comfort_focused"),
                        names(flags))
sf = merge(sf, flags[, c("hospitalization_id", comfort_col), with = FALSE],
           by = "hospitalization_id", all.x = TRUE)
sf[, algo_comfort_and := as.integer(
  if (length(comfort_col)) !is.na(get(comfort_col[1])) else FALSE)]

# stratify so the sample spans short-survival vs standard AND comfort vs not,
# so the PPV is not dominated by the majority stratum.
sf[, surv_stratum := fifelse(!is.na(hours_transition_to_death) &
                               hours_transition_to_death >= 0 &
                               hours_transition_to_death < 48,
                             "short", "standard")]
sf[, stratum := paste(surv_stratum, fifelse(algo_comfort_and == 1L, "comfort", "noncomfort"))]

# proportional allocation across strata, at least 1 per non-empty stratum
n_target = min(N_REVIEW, nrow(sf))
alloc = sf[, .(n = .N), by = stratum]
alloc[, take := pmax(1L, round(n_target * n / sum(n)))]
sample_rows = rbindlist(lapply(alloc$stratum, \(st) {
  pool = sf[stratum == st]
  k = min(alloc[stratum == st, take], nrow(pool))
  pool[sample(.N, k)]
}))
# trim any rounding overage back to the target
if (nrow(sample_rows) > n_target) sample_rows = sample_rows[sample(.N, n_target)]

worksheet = sample_rows[, .(
  hospitalization_id,
  stratum,
  algo_hospice_transition = 1L,        # every sampled row was algorithm-flagged
  algo_comfort_and,
  # --- reviewer fills these from the chart ---
  chart_hospice_confirmed = NA_integer_,   # 1 = chart confirms a hospice transition
  chart_comfort_confirmed = NA_integer_,   # 1 = chart documents comfort-focused care
  reviewer = NA_character_,
  review_date = NA_character_,
  notes = NA_character_
)]
setorder(worksheet, stratum, hospitalization_id)

# the worksheet contains PHI (encounter IDs) -> private, not shareable
out_path = file.path(private_dir,
                     paste0("aim3_chart_validation_worksheet_", site_lowercase, "_", today, ".csv"))
fwrite(worksheet, out_path)

message(sprintf("\n  Drew %d encounters for chart review across %d strata:",
                nrow(worksheet), uniqueN(worksheet$stratum)))
for (st in sort(unique(worksheet$stratum)))
  message(sprintf("    %-22s %d", st, worksheet[stratum == st, .N]))
message(sprintf("\n  Blank worksheet (PHI -> PRIVATE): %s", out_path))
message("  Reviewers fill chart_hospice_confirmed / chart_comfort_confirmed (1/0),")
message("  then re-run with AIM3_CV_SCORE=<completed worksheet path> to compute PPV.")
