# 04_aim2a3_risk_adjusted_variation.R
# Aim 2A: Hospital-level variation in discharge disposition
# Protocol v7.0, Sections 3.4.4 and 2A.1-2A.3

# Scope: this script covers all of Aim 2A.
#   2A.1  Unadjusted between-hospital variation in disposition rates
#   2A.2  Hospice-to-mortality ratios, including the 30-day death-linked ratio
#   2A.3  Risk-adjusted variation, the substantive analysis

# Approach (2A.3):
#   The protocol specifies a mixed-effects multinomial logistic regression with
#   hospital random intercepts. We implement this as a set of separate binary
#   GLMMs (lme4::glmer), each outcome vs. in_hospital_death_acute within the
#   END OF LIFE cohort, and fit the multinomial (mclogit) alongside as a
#   consistency check. See aim2a3_model_consistency for the agreement
#   diagnostics.

# Covariates, primary model (protocol 3.4.4):
#     - Age (continuous, centered)
#     - Sex
#     - Race/ethnicity (collapsed to 5-level)
#     - Admission diagnosis category (ICD-10 -> CCSR)
#     - Comorbidity burden (Elixhauser van Walraven index, PRESENT ON ADMISSION
#       diagnoses only where the site supplies a POA indicator)
#     - Admission severity (APACHE II, first 24h; plus a missingness indicator)
#     - Preexisting code status
#     - Admission source
#     - Pandemic period

# Analytic unit: ONE admission per patient, the LAST admission in the study
# period. Each patient appears once, so no model carries a patient term.

# Analytic cohorts (CHANGED 2026.10.14, see section 4):
#   eol        PRIMARY. End of life cohort: admissions ending in
#              in_hospital_death_acute, same_facility_hospice (death and alive
#              collapsed), or external_hospice_discharge. Outcomes are modeled
#              against in_hospital_death_acute, so each hospital effect is the
#              adjusted hospice to death ratio on the log odds scale.
#   full       One admission per patient, all dispositions. Two INCIDENCE fits
#              only: any hospice and in hospital death, each vs routine_discharge.
#   terminal   Death linked sites only. Died in hospital or within 30 days of
#              discharge, any disposition; adds died_30d_no_hospice. Same fits
#              as eol, plus a concordance table against eol.

# Random effects:
#     - (1 | hospital_unit)
#     - (1 | health_system)   only where a site spans 3 or more systems
#
# Fixed effects fallback (binary and multinomial; below the 3 hospital floor,
# fewer than 2 hospitals with events, or on a random effects error): hospital
# enters as a FIXED intercept, one dummy per hospital beyond the largest,
# wherever the analytic sample has 2 or more hospitals, so the fallback ORs
# stay within hospital. See add_hospital_fe in section 4. sigma^2, the MOR,
# the shrunken profiles and the risk standardized composition are GLMM only.

# Variation reporting: the protocol replaced the ICC with a profiling suite.
# The headline quantities are sigma^2 with an interval, the median odds ratio,
# risk-standardized hospital probabilities, reliability, and rank uncertainty.
# The ICC is retained as a secondary column only.

# Outputs:
#   2A.1 / 2A.2
#     - aim2a_hospital_characteristics_{site}.csv
#     - aim2a_hospital_rates_{site}.csv
#     - aim2a_hospital_ratios_{site}.csv
#     - aim2a_variation_stats_{site}.csv
#     - aim2a_unadjusted_variation_{site}.png
#   2A.3 primary
#     - aim2a3_analytic_summary_{site}.csv           — covariate distributions
#     - aim2a3_standardization_margins_{site}.csv    — for pooled-standard assembly (all admissions)
#     - aim2a3_standardization_margins_eol_{site}.csv — same, END OF LIFE cohort (for the eol standard population)
#     - aim2a3_model_coefficients_{site}.csv         — fixed effects, binary GLMMs
#     - aim2a3_variation_summary_{site}.csv          — sigma^2, MOR, p10-p90, reliability
#     - aim2a3_hospital_profile_{site}.csv           — RS probs, reliability, rank CIs
#     - aim2a3_hospital_ranef_{site}.csv             — hospital random intercepts (BLUPs)
#     - aim2a3_icc_{site}.csv                        — secondary
#     - aim2a3_rs_composition_{site}.csv             — multinomial composition
#     - aim2a3_model_consistency_{site}.csv          — binary vs. multinomial agreement
#     - aim2a3_eol_terminal_concordance_{site}.csv   — death linked sites: eol vs terminal hospital effects
#   2A.3 secondary (clinical course, separate estimand)
#     - aim2a3_model_coefficients_clinical_course_{site}.csv
#     - aim2a3_icc_clinical_course_{site}.csv
#   NOTE: the sigma^2 interval method comparison (exact vs fast, formerly the
#   TEMP-SIGMACI blocks here) now lives in validate_sigma_ci_ohsu.R, which is
#   run at OHSU only and is NOT part of the distributed consortium bundle.
#   Every site, including OHSU, reports the fast one-dimensional profile
#   interval (sigma_ci_fast below) in production output.

#   Added 2026.08.24 changes (folded 2026.08.27):
#     - aim2a3_timing_{site}.csv                     — per model fit timing (always on)
#     - aim2_multinomial_status_{site}.csv           — multinomial fit diagnostics (always on)
#     - clinical course outputs now carry course_window_hours + course_flag_provenance
#     - coefficient CSVs carry control_sampled, control_fraction, n_reference_used

# Last updated 26.10.20 (Emory, first site with ONE pandemic era; data begin
# October 2022, so every admission is Post-pandemic): the era factor was
# hard releveled to "Pre-pandemic", which stopped the stage with "'ref' must
# be an existing level". Three changes, version 04_v261020:
#   (a) pandemic_period is built from the eras OBSERVED in the analytic
#       sample, in calendar order, with the earliest observed era as the
#       reference (Pre-pandemic at full window sites, so their output is
#       unchanged).
#   (b) with only one observed era the term is constant, so it is left out
#       of base_rhs (and therefore of every primary, clinical course,
#       multinomial and fallback model), exactly as the other baseline
#       adjusters are when they do not vary. The era refits and the cross
#       era Spearman are skipped with a message (they need 2 or more eras).
#   (c) the standard population used for risk standardization is no longer
#       filtered down to the site's single era when the term is not in the
#       model (it was coerced to the site's era levels, which would have
#       silently dropped every Pre and Acute pandemic standard row).
#   run_config records era_levels_observed, era_reference and
#   era_term_in_model so the coordinating center knows which sites carry no
#   era coefficients. No change at sites with two or more eras beyond (a).
#
# Previous update 26.10.19: single hospital sites no longer crash stage 04.
# UChicago (one hospital_id) hit a stop in 2A.1: compute_var_stats returned
# NULL below 2 hospitals, the variation table came back with zero columns,
# and the summary print errored. It now returns a typed empty table, so a
# one hospital site writes header only variation stats and continues through
# the fixed effects fits; its hospital variation enters through the pooled
# coordinating center model. Builds on v261016b (hospital fixed intercept in
# the fallbacks, below); version bumped to 04_v261019.
#
# Previous update 26.10.16 (second change set): hospital FIXED intercept in the
# fixed effects fallbacks, binary and multinomial. Below the 3 hospital floor
# (or on a random effects error) the fallbacks used to drop hospital from the
# model entirely, so a two hospital site's ORs were population averaged
# rather than within hospital and their SEs ignored the clustering, and they
# were not comparable with the GLMM's conditional ORs pooled at the
# coordinating center. Hospital now enters every fallback as a fixed
# intercept wherever the analytic sample has 2 or more hospitals
# (add_hospital_fe, section 4; the reference is the largest hospital in the
# fit). The hospital contrasts are written as ordinary coefficient rows with
# term prefix hospital_fe, and two new columns, hospital_fixed_intercept and
# hospital_fe_reference, travel on every fallback coefficient row and on
# aim2_multinomial_status. model_type labels are unchanged (fixed_only,
# fixed_fallback, multinomial_fixed), so the compile census still works. No
# change to any GLMM, to sigma^2, or to the profiles. Version 04_v261016b.
#
# Earlier on 26.10.16: glmer smoke test added at startup. The UMN 261015
# run showed every glmer() call erroring instantly (environment level; most
# often a Matrix / lme4 binary mismatch after a package update), so all 13
# binary models silently degraded to fixed effects glm and the bundle shipped
# without sigma^2, MOR, or any hospital profile. An environment that cannot
# fit a toy mixed model now STOPS the run with reinstall instructions instead
# of quietly writing a bundle without the variation outputs. Version bumped
# to 04_v261016 (no statistical changes from v261014c).
#
# Previous update bpe 26.10.14 (third change set): Elixhauser from PRESENT ON
# ADMISSION diagnoses only (section 2a). The POA column is detected by name
# (poa_present in CLIF 2.1); Y / W / 1 count as present, N / U / 0 as not
# present and are dropped, exempt or blank values are kept as not stated. The
# principal diagnosis (2b) is untouched. The restriction's completeness is
# reported in aim2a3_run_config (five new rows) and a data quality flag fires
# when the indicator is absent, never populated, or more than half not
# stated. Preexisting code status is unchanged (admission plus 24 h window).
#
# Earlier on 26.10.14 (second change set, research in progress
# feedback): the PRIMARY variation model moves to the END OF LIFE cohort.
#   - eol cohort = admissions ending in in_hospital_death_acute,
#     same_facility_hospice, or external_hospice_discharge; outcomes are
#     modeled vs in_hospital_death_acute (same_facility_hospice at class (a)
#     hospitals only; external_hospice_discharge and any_hospice_disposition
#     everywhere). The hospital effect is the adjusted hospice to death ratio
#     on the log odds scale, the model version of the 2A.2 ratio. The LST and
#     ICU subgroups, the clinical course model, the multinomial (class (a)
#     hospitals, three levels), and the era refits all run inside this cohort.
#     subgroup labels: eol, eol_lst_any, eol_icu_only, eol_excl_death_48h.
#   - opportunity sensitivity eol_excl_death_48h drops in hospital deaths
#     within early_death_window_hours of admission.
#   - two INCIDENCE fits on the full one admission sample vs routine_discharge
#     (any_hospice_disposition, in_hospital_death_acute), subgroup
#     full_cohort_incidence. The only fits that still subsample controls.
#   - terminal cohort at death linked sites (died in hospital or within
#     terminal_death_window_days of discharge, any disposition; routine,
#     transfer and AMA deaths become died_30d_no_hospice), subgroup
#     terminal_30d, plus aim2a3_eol_terminal_concordance comparing its
#     hospital effects with eol.
#   - coefficient, ICC, variation, profile and ranef tables carry a
#     reference_category column; subgroup = "full_cohort" rows are no longer
#     produced.
#   - a second margins file (aim2a3_standardization_margins_eol) and a second
#     distributed standard population (config/network_standard_population_eol.csv)
#     for the end of life fits; site local until the coordinating center
#     distributes it. The all admissions population serves the incidence fits.
#   - admission_source rare level collapse now computed on the end of life rows.
#   - run_era_models switched ON (cheap inside the eol cohort).
#   - fit_family() wraps the model_base / ref_group swap the subgroup blocks
#     used to repeat by hand.
#
# Earlier on 26.10.14, after the research in progress meeting:
#   1. ONE same facility hospice outcome. same_facility_hospice_death and
#      same_facility_hospice_alive are collapsed to "same_facility_hospice"
#      in the modeling sample (model_base, top of section 4). The multinomial
#      therefore has five disposition categories against routine discharge,
#      and the per disposition binary GLMMs, the era refits, and the
#      consistency checks use the same five, so every multinomial category
#      still has its binary comparator. The split is NOT changed in 01;
#      2A.1, 2A.2, and the analytic summary still report it, and Aims 1 and
#      3 are untouched.
#   2. The patient level random effect is REMOVED entirely: no (1 | patient_id)
#      term anywhere, no estimability test, no all admissions sensitivity.
#      One admission per patient (the last) is the only analytic unit. Rows
#      with subgroup = "all_admissions_patient_re" are no longer produced and
#      the patient_re_sensitivity row leaves aim2a3_run_config.
#
# Previous update bpe 26.10.13: analytic unit changed to one admission per
# patient, the LAST admission in the study period (before this, all
# admissions with a patient random intercept where estimable, and a FIRST
# admission per patient refit as the sensitivity). subgroup = "full_cohort"
# now means the one admission per patient sample; subgroup =
# "one_admission_per_patient" rows are no longer produced. 2A.1, 2A.2, the
# analytic summary, and the standardization margins still describe ALL
# admissions. aim2a3_run_config carries the admission rule and the encounter
# and patient counts.
#
# Previous update bpe 26.10.07: control_subsample_threshold raised from 500,000
# to 1,000,000 after review of the first subsampled OHSU run (261007). OHSU
# (reference 602k) barely crossed the old threshold and took the approximation
# for a saving it did not need; at 1,000,000 a mid size site runs the exact
# protocol while truly large sites (UMN reference 2.81M) still sample. Also
# fixed aim2a3_run_config value serialization: values now format per element,
# so logicals write as TRUE/FALSE and the threshold writes as 1000000 rather
# than the 1/0/5e+05 produced by the old as.character(c(...)) coercion.
#
# Previous update bpe 26.08.28: runtime triage decision wired in. Control
# subsampling ENABLED (threshold gated, small sites unaffected); Package A
# reduction applied as switches (run_multinomial = FALSE, run_era_models =
# FALSE, run_clinical_course = TRUE); every run now writes
# aim2a3_run_config_{site}.csv so skipped model families are distinguishable
# from failed ones at the coordinating center. Skipped multinomial writes its
# status file with model_type = "disabled".
#
# Previous update bpe 26.08.27. That version =
#   26.08.18 consortium bundle (TEMP-SIGMACI removed)
#   + the 2026.08.24 and 2026.08.26 changes, folded from the project spec:
#       1. reference group control subsampling (DEFAULT OFF, switches in section 0;
#          intercept on sampled scale, RS outputs offset corrected via lp_offset)
#       2. per model timing instrumentation (always on)
#       3. windowed clinical course exposures (icu_early_flag + lst_early_flag,
#          first course_window_hours h, whole stay fallback with provenance)
#       4. exact row compression (DEFAULT OFF; weights = row_weight on every
#          binary fit, weighted sigma_ci_fast, weighted standardization draws,
#          payoff diagnostic always printed; NOT applied to multinomial or era)
#       5. locked network centering constants (always on: age 65, elix VW 0,
#          APACHE II 10, imputation 10; margins CSV still reports true site means)
#   + four BUGFIXES 2026.08.27 (search "BUGFIX 2026.08.27"):
#       a. era refits: pandemic_period constant within an era killed every era
#          fit silently -> era RHS now rebuilt per era; failures reported loudly
#       b. era result collection: nested rbindlist pattern errored the moment a
#          fit succeeded -> summary/ranef collected separately
#       c. consistency checks were dead code (all_var_summary did not exist yet)
#          -> binary profiling assembled before the multinomial block
#       d. mclogit 0.9.x has no VarCov() function -> read mn_fit$VarCov directly
#   + aim2_multinomial_status output so a missing multinomial file is
#     diagnosable from the shareable bundle.

# 0. setup ---------------------------------------------------------------------

if (!exists("data_list")) source(here::here("code", "00_setup.R"))

# packages
pkg_needed = c("lme4", "comorbidity", "mclogit")
for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(lme4)

# GLMER SMOKE TEST (2026.10.16, see header). If the environment cannot fit a
# trivial mixed model, every real fit below would silently take the fixed
# effects fallback and the site would distribute a bundle with no variation
# outputs -- which is exactly what happened at UMN on 261015. Fail loudly
# here instead, before any data are touched.
glmer_smoke = tryCatch({
  set.seed(1)
  toy = data.frame(y = rbinom(400L, 1L, 0.3),
                   g = rep(sprintf("h%02d", 1:8), each = 50L),
                   x = rnorm(400L))
  invisible(suppressWarnings(
    glmer(y ~ x + (1 | g), data = toy, family = binomial(link = "logit"),
          nAGQ = 1)))
  TRUE
}, error = function(e) conditionMessage(e))
if (!isTRUE(glmer_smoke)) {
  stop("lme4::glmer cannot fit a toy mixed model in this R session:\n    ",
       glmer_smoke,
       "\n  Every model in this script would silently fall back to fixed ",
       "effects.\n  This is usually a Matrix / lme4 binary mismatch after a ",
       "package update.\n  Reinstall BOTH from source in one session:\n",
       "    install.packages(c('Matrix', 'lme4'), type = 'source')\n",
       "  then RESTART R and re-run the pipeline.", call. = FALSE)
}
message("  glmer smoke test passed (mixed models estimable in this session).")
suppressWarnings(rm(glmer_smoke, toy))
library(ggplot2)
library(patchwork)


## Runtime configuration (added 2026.08.24 changes, folded in 2026.08.27) -------

# Change 1: reference group control subsampling. ENABLED 2026.08.28 per the
# runtime triage: model cuts alone cannot reach the 12 h target at large sites
# because the cost is the reference group itself (2.81M routine discharges at
# UMN riding in every fit). When enabled, fits whose reference group exceeds
# the threshold keep ALL events and sample controls WITHIN hospital at one
# global fraction f (about control_ratio controls per event), so hospital
# effects and sigma^2 are unchanged in expectation. The intercept shifts by
# exactly log(f); risk standardized outputs are corrected back to population
# scale via lp_offset in hospital_variation_profile. ORs need no correction;
# efficiency loss < 5% at 20 controls per event. THRESHOLD GATING means any
# site whose reference group is under 1,000,000 rows runs the untouched
# protocol, switch on or off. Threshold raised from 500,000 on 2026.10.07
# after the first subsampled OHSU run: OHSU (reference 602k) barely crossed
# the old threshold and paid the Prentice Pyke approximation for about one
# hour of saving it did not need, while at 1,000,000 OHSU runs the exact
# protocol and genuinely large sites (UMN 2.81M) still sample. First large
# site run should be compared against the UMN run 2 benchmarks (the compile
# script does this automatically).
control_subsample_enabled   = TRUE
control_subsample_threshold = 1000000L     # reference rows before sampling kicks in
control_ratio               = 20L          # target controls per event
control_subsample_seed      = 20260713L

# Package A reduction (2026.08.28, runtime triage): the multinomial and the
# era refits are DISABLED as switches, not deletions, so they can be re
# enabled later without redistributing a different script. The clinical
# course models are KEPT (Package A retains them; subsampling carries the
# runtime saving instead). Every run writes aim2a3_run_config_{site}.csv
# recording these states, so a SKIPPED model family is distinguishable from
# a FAILED one at the coordinating center.
# run_multinomial RE ENABLED 2026.10.07 at OHSU to exercise the multinomial
# on the v261007 configuration (fixed fallback below the 3 hospital floor;
# at OHSU the routine group, 602k, sits under the 1,000,000 threshold, so
# the multinomial fits UNSAMPLED — expect roughly the ~2 h cost observed on
# the 260828 run). Era refits remain off.
# run_era_models switched ON 2026.10.14: the era refits now run inside the
# end of life cohort (a few percent of admissions), so their cost is minutes.
run_clinical_course = TRUE
run_era_models      = TRUE
run_multinomial     = TRUE

# CHANGED 2026.10.14: end of life cohort definition in the RAW 01 categories
# (the two same facility categories are collapsed to one in section 4), the
# opportunity sensitivity window, and the terminal cohort window used at
# death linked sites. See section 4 for the cohorts themselves.
eol_categories_raw = c("in_hospital_death_acute", "same_facility_hospice_death",
                       "same_facility_hospice_alive", "external_hospice_discharge")
early_death_window_hours   = 48   # drop in hospital deaths this soon after admission (sensitivity)
terminal_death_window_days = 30   # terminal cohort: death within this many days of discharge

# Change 4 (2026.08.26): EXACT row compression, DEFAULT OFF. Rows identical on
# every model variable (all.vars of the formula + y, so any random effect
# structure stays valid automatically) collapse to one row with
# row_weight = count, and every binary fit runs with weights = row_weight
# (all 1s when disabled). The weighted binomial likelihood is identical up to
# a constant: fixed effects, SEs, sigma^2, and ranefs are unchanged; ONLY AIC
# shifts by a constant, so AIC is comparable only within the same compression
# setting. No approximation, no offset, no protocol amendment. Composes with
# subsampling (sample first, then compress). NOT applied to the multinomial or
# the era refits. A payoff diagnostic always prints before the model loop;
# when the diagnostic ratio is favorable, prefer compression over subsampling
# (exact beats approximate).
row_compress_enabled = FALSE

# Change 5 (2026.08.26): LOCKED network centering constants (always on).
# Site specific centers make fixed intercepts non comparable and non poolable
# across sites and put site data into the model spec. Hospital random
# intercepts, slopes, sigma^2, MOR, ICC, and RS probabilities are INVARIANT to
# centering (a constant covariate shift is absorbed by the fixed intercept);
# only the fixed intercept parameterization moves. The APACHE imputation
# constant is absorbed by the apache2_missing indicator and is locked purely
# for comparability of that coefficient. With these constants the intercept
# reads: log odds for a 65 year old, VW 0, APACHE II 10, reference level of
# every factor.
center_age_years   = 65
center_elix_vw     = 0
center_apache2     = 10
apache2_impute_val = 10   # was the site median

# Change 3: window (hours from admission) for the EARLY clinical course
# exposures used by the secondary model. 24 matches the APACHE II window.
course_window_hours = 24

# Change 2: per model timing collector (always on). Every completed fit appends
# a row; aim2a3_timing_{site}.csv is written in the export section.
fit_timings = list()
record_timing = function(stage, outcome, subgroup, rows, events, t0) {
  fit_timings[[length(fit_timings) + 1L]] <<- data.table(
    stage    = stage,
    outcome  = outcome,
    subgroup = subgroup,
    rows     = rows,
    events   = events,
    minutes  = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
  invisible(NULL)
}

# Change 1 helper: sample the reference (control) group within hospital at one
# global fraction f. Returns the rows to use plus the sampling metadata every
# caller must carry into its outputs.
subsample_controls = function(ref_dt, n_events) {
  n_ref_all = nrow(ref_dt)
  if (!isTRUE(control_subsample_enabled) ||
      n_ref_all <= control_subsample_threshold ||
      n_ref_all <= control_ratio * n_events) {
    return(list(ref = ref_dt, sampled = FALSE, fraction = 1, n_used = n_ref_all))
  }
  f = (control_ratio * n_events) / n_ref_all
  set.seed(control_subsample_seed)
  idx = ref_dt[, .I[sample.int(.N, max(1L, round(f * .N)))], by = hospital_unit]$V1
  out = ref_dt[idx]
  f_realized = nrow(out) / n_ref_all   # use the realized fraction for the offset
  list(ref = out, sampled = TRUE, fraction = f_realized, n_used = nrow(out))
}

# ADDED 2026.08.28: the run configuration travels with the shareable bundle.
# One row per switch, so the coordinating center can read a site's settings
# from the results themselves rather than from that site's console log.
# FIXED 2026.10.07: values now serialize per element. The old
# as.character(c(...)) coerced the mixed vector to numeric first, so logicals
# printed as 1/0 and the threshold as 5e+05; format() per element keeps
# TRUE/FALSE and plain integers.
run_config = data.table(
  setting = c("control_subsample_enabled", "control_subsample_threshold",
              "control_ratio", "control_subsample_seed",
              "row_compress_enabled",
              "run_clinical_course", "run_era_models", "run_multinomial",
              "course_window_hours",
              "center_age_years", "center_elix_vw", "center_apache2",
              "apache2_impute_val"),
  value   = vapply(list(control_subsample_enabled, control_subsample_threshold,
                        control_ratio, control_subsample_seed,
                        row_compress_enabled,
                        run_clinical_course, run_era_models, run_multinomial,
                        course_window_hours,
                        center_age_years, center_elix_vw, center_apache2,
                        apache2_impute_val),
                   function(x) format(x, scientific = FALSE), character(1))
)
# ENVIRONMENT MANIFEST (2026.10.16): the UMN 261015 failure (every glmer
# erroring, all models silently fixed effects) could not be diagnosed from the
# bundle because nothing recorded which package stack produced it. Record the
# versions that matter for the mixed models so the coordinating center can see
# an environment change between runs without asking the site for its console.
run_config = rbind(run_config, data.table(
  setting = c("r_version", "lme4_version", "matrix_version",
              "mclogit_version", "data_table_version"),
  value   = c(as.character(getRversion()),
              as.character(utils::packageVersion("lme4")),
              as.character(utils::packageVersion("Matrix")),
              tryCatch(as.character(utils::packageVersion("mclogit")),
                       error = function(e) "absent"),
              as.character(utils::packageVersion("data.table")))))

run_config[, `:=`(site = site_lowercase, script_version = "04_v261020")]
fwrite(run_config, shareable_path(
  paste0("aim2a3_run_config_", site_lowercase, "_", today, ".csv")))
message(sprintf("Run configuration saved (subsampling %s | clinical course %s | era %s | multinomial %s).",
                ifelse(control_subsample_enabled, "ON", "OFF"),
                ifelse(run_clinical_course, "ON", "OFF"),
                ifelse(run_era_models, "ON", "OFF"),
                ifelse(run_multinomial, "ON", "OFF")))


# 1. load cohort ---------------------------------------------------------------

message("\n===== 04_aim2a3: Risk-Adjusted Variation =====\n")

cohort_all = fread(latest_output("cohort_all_hospitalizations"),
                   colClasses = list(character = c("patient_id", "hospitalization_id")))
cohort_all[, hospitalization_id := as.character(hospitalization_id)]
cohort_all[, patient_id := as.character(patient_id)]

cohort = cohort_all[in_cohort == 1L]
n_cohort = nrow(cohort)

# patient demographics
patient = data_list[["patient"]] |>
  dplyr::select(patient_id, race_category, ethnicity_category, sex_category) |>
  dplyr::collect() |>
  as.data.table()
patient[, patient_id := as.character(patient_id)]

cohort = merge(cohort, patient, by = "patient_id", all.x = TRUE)
rm(patient)

# defensive LST flag check
for (col in c("imv_flag", "vaso_flag", "lst_any_flag")) {
  if (!col %in% names(cohort)) set(cohort, j = col, value = 0L)
}

# hospital grouping: use CCN if available, else hospital_id
if (!is.null(hospital_medicare_map)) {
  cohort = merge(cohort,
                 hospital_medicare_map[, .(hospital_id, medicare_id)],
                 by = "hospital_id", all.x = TRUE)
  cohort[is.na(medicare_id), medicare_id := hospital_id]
  cohort[, hospital_unit := medicare_id]
  message(sprintf("Using CCN as hospital unit: %d unique units.",
                  uniqueN(cohort$hospital_unit)))
} else {
  cohort[, hospital_unit := as.character(hospital_id)]
  message(sprintf("No CCN map — using hospital_id: %d unique units.",
                  uniqueN(cohort$hospital_unit)))
}

n_hospitals = uniqueN(cohort$hospital_unit)
message(sprintf("Primary cohort: %s encounters across %d hospital units.",
                format(n_cohort, big.mark = ","), n_hospitals))


# 2. Build covariates ----------------------------------------------------------


## 2a. Elixhauser comorbidity index --------------------------------------------

message("\n--- Building Elixhauser comorbidity index ---")

# detect a principal/priority/rank column once (varies across the consortium);
# we'll carry it through so 2b can pick the true principal dx, not a random row.
dx_schema = data_list[["hospital_diagnosis"]] |> head(0) |>
  dplyr::collect() |> as.data.table()
dx_flag_cols = grep("primary|principal|priority|rank|sequence|seq|position",
                    names(dx_schema), ignore.case = TRUE, value = TRUE)
flag_col = if (length(dx_flag_cols)) dx_flag_cols[1] else NA_character_

# CHANGED 2026.10.14: detect the present on admission column the same way
# (CLIF 2.1 names it poa_present; older extracts vary). It restricts the
# COMORBIDITY index below; the principal diagnosis in 2b is untouched.
poa_cols = grep("poa|present_on_admission", names(dx_schema), ignore.case = TRUE, value = TRUE)
poa_col  = if (length(poa_cols)) poa_cols[1] else NA_character_

dx = data_list[["hospital_diagnosis"]] |>
  dplyr::select(dplyr::all_of(unique(c(
    "hospitalization_id", "diagnosis_code", "diagnosis_code_format", dx_flag_cols, poa_cols)))) |>
  dplyr::collect() |>
  as.data.table()

dx[, hospitalization_id := as.character(hospitalization_id)]

# keep only ICD-10-CM codes. Normalize the format string first: across ~40 sites
# the value varies (ICD10CM, ICD-10-CM, icd10, 10, ...), so strip non-alphanumerics
# and match rather than enumerating spellings (a bare list silently matched zero
# rows at sites coding the dashless "ICD10CM").
dx[, fmt_norm := gsub("[^a-z0-9]", "", tolower(trimws(diagnosis_code_format)))]
dx_icd10 = dx[fmt_norm %in% c("icd10cm", "icd10")]

# clean codes: remove dots, upper-case, trim
dx_icd10[, code_clean := toupper(gsub("\\.", "", trimws(diagnosis_code)))]

# keep only encounters in the cohort
dx_icd10 = dx_icd10[hospitalization_id %in% cohort$hospitalization_id]

n_dx = nrow(dx_icd10)
n_enc_with_dx = uniqueN(dx_icd10$hospitalization_id)
message(sprintf("  ICD-10 codes: %s across %s encounters.",
                format(n_dx, big.mark = ","),
                format(n_enc_with_dx, big.mark = ",")))

# CHANGED 2026.10.14: comorbidity from PRESENT ON ADMISSION diagnoses only.
# Discharge diagnoses include conditions acquired during the stay (renal
# failure, sepsis, delirium), so without this restriction part of the harm
# being measured enters the expected risk, and the hospital's coding depth
# becomes part of its case mix. The principal diagnosis (2b) is not
# restricted: it is the reason for admission by definition and is chosen by
# the primary flag, not by POA.
#
# Value handling follows the AHRQ convention for the CMS indicator: Y / W
# (clinically undetermined) count as present, N / U (insufficient
# documentation) as not present, and exempt or blank values are not stated.
# CLIF ETLs usually deliver a 0/1 poa_present, which maps the same way. A code
# is DROPPED only when the site explicitly flags it as not present; a code
# whose POA status is not stated is kept, so a sparsely populated column
# weakens the restriction rather than discarding comorbidities. The shares
# are reported here and in aim2a3_run_config so the coordinating center can
# see how complete each site's restriction was. No column, or a column that
# is never populated, falls back to all diagnoses with a data quality flag.
.poa_norm = function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(fifelse(is.na(x), NA_integer_, as.integer(x != 0)))
  v = tolower(trimws(as.character(x)))
  fcase(v %in% c("1", "y", "yes", "true", "t", "w"), 1L,
        v %in% c("0", "n", "no", "false", "f", "u"), 0L,
        default = NA_integer_)
}
elix_poa_restricted = FALSE
poa_share = c(present = NA_real_, not_present = NA_real_, not_stated = NA_real_)
if (!is.na(poa_col) && n_dx > 0) {
  dx_icd10[, poa_norm := .poa_norm(get(poa_col))]
  poa_share = c(present     = dx_icd10[poa_norm %in% 1L, .N] / n_dx,
                not_present = dx_icd10[poa_norm %in% 0L, .N] / n_dx,
                not_stated  = dx_icd10[is.na(poa_norm), .N] / n_dx)
  if (poa_share[["not_stated"]] < 1) {
    elix_poa_restricted = TRUE
    dx_elix = dx_icd10[is.na(poa_norm) | poa_norm == 1L]
    message(sprintf("  Elixhauser restricted to present-on-admission codes ('%s'): %.1f%% present, %.1f%% not present (dropped), %.1f%% not stated (kept).",
                    poa_col, 100 * poa_share[["present"]],
                    100 * poa_share[["not_present"]], 100 * poa_share[["not_stated"]]))
    if (poa_share[["not_stated"]] > 0.5)
      qi_flag("Present-on-admission indicator sparsely populated",
              sprintf(paste0(
                "%.1f%% of ICD-10 codes have no POA status ('%s'). Those codes are kept in the\n",
                "Elixhauser index, so the POA restriction is only partial at this site and the\n",
                "comorbidity score may still carry in-hospital diagnoses. Review with the data team."),
                100 * poa_share[["not_stated"]], poa_col))
  } else {
    dx_elix = dx_icd10
    qi_flag("Present-on-admission indicator never populated",
            sprintf(paste0(
              "Column '%s' exists but holds no usable value, so the Elixhauser index uses ALL\n",
              "discharge diagnoses, including conditions acquired during the stay. The site's\n",
              "comorbidity adjustment is not comparable with POA-restricted sites."), poa_col))
  }
} else {
  dx_elix = dx_icd10
  qi_flag("No present-on-admission indicator",
          paste0("hospital_diagnosis has no POA column, so the Elixhauser index uses ALL\n",
                 "discharge diagnoses, including conditions acquired during the stay. The site's\n",
                 "comorbidity adjustment is not comparable with POA-restricted sites."))
}

# compute Elixhauser using the comorbidity package
elix = comorbidity::comorbidity(
  x   = dx_elix[, .(id = hospitalization_id, code = code_clean)],
  id  = "id",
  code = "code",
  map = "elixhauser_icd10_quan",
  assign0 = TRUE
)

# van Walraven weighted composite score
elix_scores = comorbidity::score(elix, weights = "vw", assign0 = TRUE)

elix_dt = data.table(
  hospitalization_id = elix$id,
  elixhauser_vw      = as.integer(elix_scores)
)

# count of comorbidities (unweighted)
flag_cols_elix = setdiff(names(elix), "id")
elix_dt[, elixhauser_count := as.integer(rowSums(
  as.matrix(elix[, flag_cols_elix, drop = FALSE])
))]

cohort = merge(cohort, elix_dt, by = "hospitalization_id", all.x = TRUE)
cohort[is.na(elixhauser_vw), elixhauser_vw := 0L]
cohort[is.na(elixhauser_count), elixhauser_count := 0L]

message(sprintf("  Elixhauser (%s): median VW score = %d, mean = %.1f",
                if (elix_poa_restricted) "POA codes" else "ALL codes",
                median(cohort$elixhauser_vw),
                mean(cohort$elixhauser_vw)))

rm(elix, elix_scores, elix_dt, flag_cols_elix, dx_elix)


## 2b. Admission diagnosis categories — AHRQ CCSR body system ------------------

message("\n--- Building diagnosis categories (AHRQ CCSR body system) ---")

# Classify the PRINCIPAL diagnosis using AHRQ CCSR: the inpatient-default CCSR
# category (what AHRQ assigns when a code is the principal/first-listed inpatient
# dx), rolled up to 22 body systems via the 3-letter category prefix. Standard,
# reproducible, and handles the long tail — preferred over hand-rolled prefix
# rules for a general inpatient cohort. The detailed 530-level CCSR category is
# retained (ccsr_category) so finer grouping is a column-swap away later.

ccsr_fn = list.files(
  unique(c(tryCatch(here::here(),                error = function(e) NULL),
           tryCatch(here::here("code"),          error = function(e) NULL),
           tryCatch(here::here("code", "utils"), error = function(e) NULL),
           if (exists("project_root")) c(project_root,
                                         file.path(project_root, "code"),
                                         file.path(project_root, "code", "utils")),
           getwd())) |> (\(x) x[!is.na(x) & dir.exists(x)])(),
  pattern = "^function_ccsr_bodysystem\\.R$", full.names = TRUE, recursive = TRUE)
if (!length(ccsr_fn))
  stop("function_ccsr_bodysystem.R not found in the project tree — ",
       "place it in code/utils/.")
if (length(ccsr_fn) > 1) message("  WARNING: multiple function_ccsr_bodysystem.R found; using ", ccsr_fn[1])
source(ccsr_fn[1])
ccsr_lk = load_ccsr_lookup()

message(sprintf("  Principal-dx selection: %s",
                if (!is.na(flag_col)) paste0("flag '", flag_col, "'")
                else "first-listed in data order (no primary flag found)"))

principal = select_principal_dx(dx_icd10, flag_col = flag_col)

ccsr = assign_ccsr_bodysystem(principal$principal_code, lookup = ccsr_lk)
principal[, `:=`(dx_category   = ccsr$body_system_short,   # 22-level body system
                 ccsr_category = ccsr$ccsr_category)]      # 530-level, for later

# coverage diagnostics — multi-year cohorts will have some retired/invalid codes
# (Unmapped) and some manifestation/sequela codes flagged unacceptable-principal
# (XXX). Both are surfaced rather than hidden so you can judge data quality.
n_unm = principal[dx_category == "Unmapped", .N]
n_xxx = principal[dx_category == "Unacceptable principal dx", .N]
message(sprintf("  CCSR mapped %s principal dx | unmapped: %s (%.2f%%) | XXX: %s (%.2f%%)",
                format(nrow(principal), big.mark = ","),
                format(n_unm, big.mark = ","), 100 * n_unm / max(nrow(principal), 1L),
                format(n_xxx, big.mark = ","), 100 * n_xxx / max(nrow(principal), 1L)))

cohort = merge(cohort,
               principal[, .(hospitalization_id, dx_category, ccsr_category)],
               by = "hospitalization_id", all.x = TRUE)
cohort[is.na(dx_category), dx_category := "Unmapped"]
cohort[is.na(ccsr_category), ccsr_category := "UNMAPPED"]

# Swap dx_category to the ACUITY grouping if the crosswalk is supplied -- the
# "column-swap" the note above anticipated. Acuity (cancer, dementia, frailty,
# chronic organ failure, ...) is more prognostically relevant for end-of-life
# disposition than body system. Body system is retained as dx_body_system; the
# model covariate dx_category becomes the acuity grouping. Unmapped CCSR codes
# fall to "Other/residual".
acuity_cw = tryCatch(load_ccsr_acuity_crosswalk(), error = function(e) NULL)
if (!is.null(acuity_cw) && all(c("ccsr_category", "acuity_category") %in% names(acuity_cw))) {
  acuity_cw = unique(as.data.table(acuity_cw)[, .(ccsr_category = as.character(ccsr_category),
                                                  acuity_category = as.character(acuity_category))])
  cohort[, dx_body_system := dx_category]
  # match() lookup (not an update-join): yields NA for any CCSR not in the
  # crosswalk, so unmapped codes correctly fall through to Other/residual below.
  cohort[, dx_category := acuity_cw$acuity_category[match(ccsr_category, acuity_cw$ccsr_category)]]
  n_resid = cohort[is.na(dx_category), .N]
  cohort[is.na(dx_category), dx_category := "Other/residual"]
  message(sprintf("  dx_category set to ACUITY grouping (%d levels); %s encounters with an\n  unmapped/uncrosswalked CCSR fell to Other/residual.",
                  uniqueN(cohort$dx_category), format(n_resid, big.mark = ",")))
  message("  Acuity dx_category distribution:")
} else {
  message("  CCSR body-system distribution (acuity crosswalk not supplied):")
}
dx_tab = cohort[, .N, by = dx_category][order(-N)]
for (i in seq_len(nrow(dx_tab))) {
  message(sprintf("    %-28s %s (%4.1f%%)",
                  dx_tab$dx_category[i],
                  format(dx_tab$N[i], big.mark = ","),
                  100 * dx_tab$N[i] / n_cohort))
}

rm(dx, dx_icd10, principal, ccsr, ccsr_lk, n_dx, n_enc_with_dx)
suppressWarnings(rm(acuity_cw))
suppressWarnings(rm(dx_schema, dx_flag_cols))


## 2c. APACHE II severity score ------------------------------------------------

message("\n--- Computing APACHE II scores ---")

# source the function (it's not loaded by default).
# The function lives alongside these scripts (in code/). The original fixed-path
# list with recursive = FALSE missed it whenever project_root points ABOVE code/
# (e.g. project_root = .../025_hospice_flip and the file is in .../code/).
# Search recursively from several anchors, dedupe, and take the first match.
apache_search_roots = unique(c(
  tryCatch(here::here(),                error = function(e) NULL),
  tryCatch(here::here("code"),          error = function(e) NULL),
  tryCatch(here::here("code", "utils"), error = function(e) NULL),
  tryCatch(here::here("R"),             error = function(e) NULL),
  project_root,
  file.path(project_root, "code"),
  file.path(project_root, "code", "utils"),
  file.path(project_root, "R"),
  file.path(project_root, "scripts"),
  getwd()
))
apache_search_roots = apache_search_roots[!is.na(apache_search_roots) &
                                            dir.exists(apache_search_roots)]

apache_fn_path = unique(list.files(
  apache_search_roots,
  pattern = "^(function_apache_II|fn_apache2)\\.R$",  # anchored: no backups/tests
  full.names = TRUE, recursive = TRUE
))

if (length(apache_fn_path) > 0) {
  if (length(apache_fn_path) > 1) message("  WARNING: multiple APACHE function files found; using ", apache_fn_path[1])
  source(apache_fn_path[1])
  message(sprintf("  Sourced APACHE II function from: %s", apache_fn_path[1]))

  # Protocol 2A.3: baseline severity must come from a prespecified EARLY window
  # (first 6 or 24 hours). 24h is also APACHE II's canonical definition (worst
  # values in the first 24 hours). Previously 48h.
  apache = calculate_apache2(data_list, cohort, window_hours = 24)

  cohort = merge(cohort,
                 apache[, .(hospitalization_id, apache2_score, apache2_complete)],
                 by = "hospitalization_id", all.x = TRUE)
  cohort[is.na(apache2_score), apache2_score := NA_real_]

  n_apache = sum(!is.na(cohort$apache2_score))
  message(sprintf("  APACHE II scores available for %s of %s encounters (%4.1f%%).",
                  format(n_apache, big.mark = ","),
                  format(n_cohort, big.mark = ","),
                  100 * n_apache / n_cohort))
  rm(apache)

} else {
  warning("APACHE II function file not found (function_apache_II.R / fn_apache2.R). ",
          "Searched recursively under:\n    ",
          paste(apache_search_roots, collapse = "\n    "),
          "\n  Proceeding without APACHE II (it will be dropped from the model).",
          call. = FALSE)
  cohort[, apache2_score := NA_real_]
  cohort[, apache2_complete := FALSE]
}


## 2d. Race/ethnicity and sex (shared recodes from 00) --------------------------

# race_eth: 5-level combined race/ethnicity with Other + Unknown merged.
# sex_female: binary Female vs not-female (Unknown grouped with not-female).
# These match the Table-1 reporting definitions so modeled and reported groups
# are identical, and the smallest cells are collapsed for adequate cell sizes.
recode_race_eth(cohort)
recode_sex_female(cohort)


## 2e. Centered and scaled continuous covariates --------------------------------

# Change 5: centering uses the LOCKED network constants from section 0, not
# site means, so fixed intercepts are directly comparable and poolable across
# sites and no site data enters the model spec.
cohort[, age_centered := age_at_admission - center_age_years]
cohort[, elix_centered := elixhauser_vw - center_elix_vw]

# APACHE II: center at the locked constant; impute the locked constant for
# missing (change 5; was the site median).
# Guard the all-missing case: at sites where the APACHE function isn't found
# (apache_fn_path empty) or required inputs are absent, apache2_score is wholly
# NA. The site median is still computed to DETECT that case (and for the
# standardization margins), but it no longer enters the model.
apache_med = median(cohort$apache2_score, na.rm = TRUE)   # detection + margins only
apache2_available = is.finite(apache_med)

if (apache2_available) {
  cohort[, apache2_imputed  := fifelse(is.na(apache2_score), apache2_impute_val, apache2_score)]
  cohort[, apache2_centered := apache2_imputed - center_apache2]
  cohort[, apache2_missing  := fifelse(is.na(apache2_score), 1L, 0L)]
} else {
  warning("No non-missing APACHE II scores at this site — ",
          "dropping APACHE from the covariate set for all models.", call. = FALSE)
  cohort[, apache2_imputed  := NA_real_]
  cohort[, apache2_centered := 0]   # constant; excluded from formula below
  cohort[, apache2_missing  := 1L]
}

message(sprintf("  APACHE II in model: %s", if (apache2_available) "yes" else "NO (all-missing at site)"))


# 2b. Baseline case-mix adjusters added in protocol v7.0 (§2A.3) ---------------
# Both are fixed at (or before) admission, so they belong in the PRIMARY model.

## Preexisting code status: the code status documented at or before admission +24h,
## using the same vocabulary as 02. A patient who arrives already DNR/AND/comfort
## is a materially different baseline risk. Standing pre-admission orders count.
pcs = tryCatch({
  cs = data_list[["code_status"]] |>
    dplyr::select(patient_id, start_dttm, code_status_category) |>
    dplyr::collect() |> as.data.table()
  cs[, patient_id := as.character(patient_id)]
  cs = merge(cs, cohort[, .(hospitalization_id, patient_id, admission_dttm)],
             by = "patient_id", allow.cartesian = TRUE)
  cs = cs[start_dttm <= admission_dttm + lubridate::hours(24)]
  setorder(cs, hospitalization_id, start_dttm)
  cs = cs[, .SD[.N], by = hospitalization_id]            # last order at/near admission
  cs[, cat_lower := tolower(trimws(code_status_category))]
  full_cats0    = c("full", "presume full")
  limited_cats0 = c("dnr", "dnar", "udnr", "dnr/dni", "dnar/dni", "special/partial")
  comfort_cats0 = c("and", "cmo", "comfort care", "comfort measures only", "comfort")
  cs[, preexisting_code_status := fcase(
    cat_lower %in% full_cats0,    "Full",
    cat_lower %in% comfort_cats0, "Comfort",
    cat_lower %in% limited_cats0, "Limited",
    default = "Other")]
  cs[, .(hospitalization_id, preexisting_code_status)]
}, error = function(e) { message("  Preexisting code status unavailable: ", e$message); NULL })

if (!is.null(pcs)) {
  cohort = merge(cohort, pcs, by = "hospitalization_id", all.x = TRUE)
  cohort[is.na(preexisting_code_status), preexisting_code_status := "None documented"]
  ref_pcs = names(which.max(table(cohort$preexisting_code_status)))  # most common = ref
  cohort[, preexisting_code_status := relevel(factor(preexisting_code_status), ref = ref_pcs)]
  message(sprintf("  Preexisting code status (ref = %s): %s", ref_pcs,
                  paste(sprintf("%s=%s", levels(cohort$preexisting_code_status),
                                table(cohort$preexisting_code_status)[levels(cohort$preexisting_code_status)]),
                        collapse = "  ")))
} else {
  cohort[, preexisting_code_status := factor("None documented")]  # constant -> excluded
}

## Admission source: CLIF field naming varies, so detect a source-like column and
## fall back to admission_type_category; if nothing usable, omit with a message.
adm_src_col = tryCatch({
  sch  = names(data_list[["hospitalization"]] |> head(0) |> dplyr::collect())
  cand = intersect(c("admission_source_category", "admission_source_name",
                     "admission_source", "admission_type_category"), sch)
  if (length(cand)) cand[1] else NA_character_
}, error = function(e) NA_character_)

if (!is.na(adm_src_col)) {
  asrc = data_list[["hospitalization"]] |>
    dplyr::select(hospitalization_id, dplyr::all_of(adm_src_col)) |>
    dplyr::collect() |> as.data.table()
  setnames(asrc, adm_src_col, "admission_source")
  asrc[, hospitalization_id := as.character(hospitalization_id)]
  asrc[, admission_source   := as.character(admission_source)]
  asrc[is.na(admission_source) | !nzchar(admission_source), admission_source := "Unknown"]
  cohort = merge(cohort, asrc, by = "hospitalization_id", all.x = TRUE)
  cohort[is.na(admission_source), admission_source := "Unknown"]
  # CHANGED 2026.10.14: the rare level floor is applied on the END OF LIFE rows
  # (the primary modeling sample, a few percent of admissions), not on all
  # admissions. A level with 50 admissions overall may have one or two among
  # deaths and hospice discharges, which invites separation in the primary
  # fits. Levels kept here are kept for every cohort, so the factor is coded
  # the same way in all fits and in both margins files.
  keep = cohort[outcome_category %in% eol_categories_raw, .N, by = admission_source
                ][N >= 50, admission_source]   # collapse rare
  cohort[!admission_source %in% keep, admission_source := "Other"]
  ref_src = names(which.max(table(cohort$admission_source)))
  cohort[, admission_source := relevel(factor(admission_source), ref = ref_src)]
  message(sprintf("  Admission source from '%s' (ref = %s): %d levels.",
                  adm_src_col, ref_src, nlevels(cohort$admission_source)))
} else {
  cohort[, admission_source := factor("Unknown")]   # constant -> excluded
  message("  Admission source: no source-like field in the hospitalization table; omitted.")
}


## 2g. Windowed clinical course exposures (change 3, protocol 2A.3) -------------
# The secondary clinical course model conditions on EARLY organ support, not
# whole stay exposure (a whole stay flag can postdate, and be caused by, the
# hospital's own practice). Definitions mirror script 01 exactly plus the
# window (course_window_hours from admission, set in section 0):
#   ICU:  adt location_category icu with 0 <= hours_to_icu <= window
#         (01 uses >= 0 without the cap)
#   IMV:  respiratory_support device_category imv,
#         recorded_dttm <= admission + window (no lower bound, so support
#         already running at admission counts)
#   Vaso: medication_admin_continuous med_category in the list below,
#         admin_dttm <= admission + window
# RRT stays excluded (as in 01; crrt_therapy coverage too sparse across the
# consortium). Early LST = IMV union vasopressor; partial unions are never
# used (if either component fails, BOTH fall back to the whole stay flag).
# If a source table is missing the affected flag falls back to the whole stay
# flag with a loud message, and the clinical course outputs carry
# course_window_hours + course_flag_provenance so the fallback is visible at
# the coordinating center.

message(sprintf("\n--- Building windowed clinical course exposures (first %d h) ---",
                course_window_hours))

# KEEP IN SYNC WITH 01 (vasopressor med_category list)
course_vaso_categories = c("norepinephrine", "vasopressin", "epinephrine",
                           "phenylephrine", "dopamine")

.window_ids = function(dt_expr) tryCatch(dt_expr(), error = function(e) {
  message("    source unavailable: ", conditionMessage(e)); NULL })

icu_prov = sprintf("windowed_%dh", course_window_hours)
lst_prov = sprintf("windowed_%dh", course_window_hours)

# ICU within the window
icu_early_ids = .window_ids(function() {
  adt = data_list[["adt"]] |>
    dplyr::select(hospitalization_id, in_dttm, location_category) |>
    dplyr::collect() |> as.data.table()
  adt[, hospitalization_id := as.character(hospitalization_id)]
  adt = adt[hospitalization_id %chin% cohort$hospitalization_id &
              tolower(trimws(location_category)) == "icu"]
  adt = merge(adt, cohort[, .(hospitalization_id, admission_dttm)],
              by = "hospitalization_id")
  adt[, hours_to_icu := as.numeric(difftime(in_dttm, admission_dttm, units = "hours"))]
  unique(adt[hours_to_icu >= 0 & hours_to_icu <= course_window_hours,
             hospitalization_id])
})

if (!is.null(icu_early_ids)) {
  cohort[, icu_early_flag := fifelse(hospitalization_id %chin% icu_early_ids, 1L, 0L)]
} else {
  message(sprintf("    FALLBACK: adt unavailable — icu_early_flag set to the WHOLE STAY icu_flag. Clinical course outputs will carry course_flag_provenance = icu=whole_stay_fallback."))
  cohort[, icu_early_flag := icu_flag]
  icu_prov = "whole_stay_fallback"
}

# IMV within the window (no lower bound)
imv_early_ids = .window_ids(function() {
  rs = data_list[["respiratory_support"]] |>
    dplyr::select(hospitalization_id, recorded_dttm, device_category) |>
    dplyr::collect() |> as.data.table()
  rs[, hospitalization_id := as.character(hospitalization_id)]
  rs = rs[hospitalization_id %chin% cohort$hospitalization_id &
            tolower(trimws(device_category)) == "imv"]
  rs = merge(rs, cohort[, .(hospitalization_id, admission_dttm)],
             by = "hospitalization_id")
  unique(rs[recorded_dttm <= admission_dttm + lubridate::hours(course_window_hours),
            hospitalization_id])
})

# vasopressors within the window
vaso_early_ids = .window_ids(function() {
  mac = data_list[["medication_admin_continuous"]] |>
    dplyr::select(hospitalization_id, admin_dttm, med_category) |>
    dplyr::collect() |> as.data.table()
  mac[, hospitalization_id := as.character(hospitalization_id)]
  mac = mac[hospitalization_id %chin% cohort$hospitalization_id &
              tolower(trimws(med_category)) %in% course_vaso_categories]
  mac = merge(mac, cohort[, .(hospitalization_id, admission_dttm)],
              by = "hospitalization_id")
  unique(mac[admin_dttm <= admission_dttm + lubridate::hours(course_window_hours),
             hospitalization_id])
})

if (!is.null(imv_early_ids) && !is.null(vaso_early_ids)) {
  cohort[, imv_early_flag  := fifelse(hospitalization_id %chin% imv_early_ids,  1L, 0L)]
  cohort[, vaso_early_flag := fifelse(hospitalization_id %chin% vaso_early_ids, 1L, 0L)]
  cohort[, lst_early_flag  := pmax(imv_early_flag, vaso_early_flag)]
} else {
  # never a partial union: if either component fails, all three fall back
  message("    FALLBACK: respiratory_support or medication_admin_continuous unavailable — lst_early_flag set to the WHOLE STAY lst_any_flag (never a partial union). Clinical course outputs will carry course_flag_provenance = lst=whole_stay_fallback.")
  cohort[, imv_early_flag  := imv_flag]
  cohort[, vaso_early_flag := vaso_flag]
  cohort[, lst_early_flag  := lst_any_flag]
  lst_prov = "whole_stay_fallback"
}

course_flag_provenance = sprintf("icu=%s;lst=%s", icu_prov, lst_prov)

message(sprintf("  Early flags (first %d h): ICU %.1f%% (whole stay %.1f%%) | LST %.1f%% (whole stay %.1f%%) | provenance: %s",
                course_window_hours,
                100 * mean(cohort$icu_early_flag, na.rm = TRUE),
                100 * mean(cohort$icu_flag, na.rm = TRUE),
                100 * mean(cohort$lst_early_flag, na.rm = TRUE),
                100 * mean(cohort$lst_any_flag, na.rm = TRUE),
                course_flag_provenance))

suppressWarnings(rm(icu_early_ids, imv_early_ids, vaso_early_ids))


# 2f. Federation standardization frame (protocol 2A.3, Aim 2A federation) -------
# For hospital estimates to be comparable ACROSS sites, every site must risk-
# standardize onto the SAME patient population. The coordinating center
# distributes a network standard population of RAW covariate profiles; each site
# reapplies ITS OWN centering/level mapping (captured here) and predicts on it.
# Absent a distributed standard, std_frame stays NULL and the profiling falls
# back to site-local standardization -- so single-site OHSU behavior is unchanged.
# Change 5: std_transform carries the LOCKED constants (field names unchanged
# so apply_std_transform works as before). The standardization margins CSV
# below still reports TRUE site means for pooled standard assembly.
std_transform = list(
  apache2_available = apache2_available,
  age_mean        = center_age_years,
  elix_mean       = center_elix_vw,
  apache_med      = apache2_impute_val,
  apache_imp_mean = if (apache2_available) center_apache2 else NA_real_,
  dx_levels       = sort(unique(as.character(cohort$dx_category))),
  race_levels     = levels(factor(cohort$race_eth)),
  sex_levels      = sort(unique(as.character(cohort$sex_female))),
  era_levels      = levels(factor(cohort$pandemic_period)),
  pcs_levels      = if ("preexisting_code_status" %in% names(cohort)) levels(cohort$preexisting_code_status) else NULL,
  src_levels      = if ("admission_source" %in% names(cohort)) levels(cohort$admission_source) else NULL
)

apply_std_transform = function(raw) {
  d = as.data.table(copy(raw))
  if (!"elixhauser_vw" %in% names(d)) stop("standard population missing elixhauser_vw")
  d[, age_centered  := age_at_admission - std_transform$age_mean]
  d[, elix_centered := elixhauser_vw     - std_transform$elix_mean]
  if (std_transform$apache2_available) {
    d[, apache2_imp      := fifelse(is.na(apache2_score), std_transform$apache_med, apache2_score)]
    d[, apache2_centered := apache2_imp - std_transform$apache_imp_mean]
    d[, apache2_missing  := fifelse(is.na(apache2_score), 1L, 0L)]
  } else {
    d[, apache2_centered := 0]; d[, apache2_missing := 1L]
  }
  # coerce categoricals to the SITE's factor levels; out-of-vocab -> NA -> dropped
  d[, dx_category     := factor(as.character(dx_category),     levels = std_transform$dx_levels)]
  d[, race_eth        := factor(as.character(race_eth),        levels = std_transform$race_levels)]
  if (!"sex_female" %in% names(d) && "sex_category" %in% names(d))
    d[, sex_female := fifelse(tolower(trimws(sex_category)) == "female", "Female", "Not female")]
  d[, sex_female      := factor(as.character(sex_female),      levels = std_transform$sex_levels)]
  # 2026.10.20: at a single era site the era term is not in the model, so the
  # standard population must not be filtered to that era; give every standard
  # row the site's one level (it does not enter any prediction)
  if (length(std_transform$era_levels) == 1L) d[, pandemic_period := std_transform$era_levels]
  d[, pandemic_period := factor(as.character(pandemic_period), levels = std_transform$era_levels)]
  if (!is.null(std_transform$pcs_levels) && "preexisting_code_status" %in% names(d))
    d[, preexisting_code_status := factor(as.character(preexisting_code_status), levels = std_transform$pcs_levels)]
  if (!is.null(std_transform$src_levels) && "admission_source" %in% names(d))
    d[, admission_source := factor(as.character(admission_source), levels = std_transform$src_levels)]
  keyvars = intersect(c("dx_category", "race_eth", "sex_female", "pandemic_period",
                        "preexisting_code_status", "admission_source"), names(d))
  d[complete.cases(d[, ..keyvars])]
}

standard_pop_raw = tryCatch(load_standard_population(), error = function(e) NULL)
if (!is.null(standard_pop_raw)) {
  std_frame = tryCatch(apply_std_transform(standard_pop_raw),
                       error = function(e) { message("  Standard-pop mapping failed: ",
                                                     conditionMessage(e)); NULL })
  if (!is.null(std_frame) && nrow(std_frame) > 0) {
    standardization_source = paste0("network:", standard_pop_version)
    message(sprintf("  Standardizing on NETWORK standard population (%s rows, version '%s').",
                    format(nrow(std_frame), big.mark = ","), standard_pop_version))
  } else {
    message("  Standard population unusable after level mapping -> site-local standardization.")
    std_frame = NULL; standardization_source = "site_local"
  }
} else {
  std_frame = NULL; standardization_source = "site_local"
}

# CHANGED 2026.10.14: the END OF LIFE fits standardize onto an end of life
# standard population, not the all admissions one above. The all admissions
# population describes a case mix in which almost nobody dies; predicting
# P(hospice | death or hospice) for it is meaningless and would distort the
# risk standardized probabilities. A second file, which the coordinating
# center assembles from the aim2a3_standardization_margins_eol family, is
# looked for next to the first (config key network_standard_population_eol_file,
# version key network_standard_population_eol_version). Absent, the end of
# life and terminal fits standardize site locally on their own analytic
# sample. The all admissions population still serves the incidence fits.
.cfg = function(key, default) if (!is.null(config[[key]])) config[[key]] else default
standard_pop_eol_path    = .cfg("network_standard_population_eol_file",
                                file.path(project_root, "config", "network_standard_population_eol.csv"))
standard_pop_eol_version = .cfg("network_standard_population_eol_version", "site_local")
std_frame_eol = NULL; standardization_source_eol = "site_local"
if (file.exists(standard_pop_eol_path)) {
  std_frame_eol = tryCatch(apply_std_transform(fread(standard_pop_eol_path)),
                           error = function(e) { message("  End of life standard-pop mapping failed: ",
                                                         conditionMessage(e)); NULL })
  if (!is.null(std_frame_eol) && nrow(std_frame_eol) > 0) {
    standardization_source_eol = paste0("network_eol:", standard_pop_eol_version)
    message(sprintf("  End of life fits standardize on the NETWORK end of life population (%s rows, version '%s').",
                    format(nrow(std_frame_eol), big.mark = ","), standard_pop_eol_version))
  } else {
    std_frame_eol = NULL
    message("  End of life standard population unusable after level mapping -> site-local standardization.")
  }
} else {
  message("  End of life standard population: not supplied -> end of life fits standardize site-locally.")
}

# Site covariate margins (shareable, non-identifying) so the coordinating center
# can assemble the pooled network standard population (design decision 1).
# CHANGED 2026.10.14: built twice, for all admissions and for the end of life
# rows, so the coordinating center can assemble both standard populations.
# An apache2 sd row is added (cc_03 previously had to assume one).
.marg_cat = function(col, dt) {
  if (!col %in% names(dt)) return(NULL)
  # group by the level only (a constant cannot go in `by`, which must be full
  # length); attach the constant covariate name afterward.
  out = dt[, .(value = round(.N / nrow(dt), 5)),
           by = .(level = as.character(get(col)))]
  out[, covariate := col]
  out[, .(covariate, level, value)]
}
# change 5: margins report TRUE site means (not the locked model centers) so
# the coordinating center can assemble the pooled network standard population
.build_margins = function(dt) {
  out = rbindlist(list(
    data.table(covariate = "age_at_admission", level = "mean", value = round(mean(dt$age_at_admission, na.rm = TRUE), 3)),
    data.table(covariate = "age_at_admission", level = "sd",   value = round(sd(dt$age_at_admission, na.rm = TRUE), 3)),
    data.table(covariate = "elixhauser_vw",    level = "mean", value = round(mean(dt$elixhauser_vw, na.rm = TRUE), 3)),
    data.table(covariate = "elixhauser_vw",    level = "sd",   value = round(sd(dt$elixhauser_vw, na.rm = TRUE), 3)),
    if (apache2_available) data.table(covariate = "apache2", level = "mean",   value = round(mean(dt$apache2_score, na.rm = TRUE), 3)),
    if (apache2_available) data.table(covariate = "apache2", level = "sd",     value = round(sd(dt$apache2_score, na.rm = TRUE), 3)),
    if (apache2_available) data.table(covariate = "apache2", level = "median", value = round(median(dt$apache2_score, na.rm = TRUE), 3)),
    .marg_cat("sex_female", dt), .marg_cat("race_eth", dt), .marg_cat("dx_category", dt),
    .marg_cat("pandemic_period", dt), .marg_cat("preexisting_code_status", dt), .marg_cat("admission_source", dt)
  ), fill = TRUE)
  out[, `:=`(site = site_lowercase, n_cohort = nrow(dt))]
  out
}
std_margins = .build_margins(cohort)
fwrite(std_margins, shareable_path(
  paste0("aim2a3_standardization_margins_", site_lowercase, "_", today, ".csv")))
message("  Saved aim2a3_standardization_margins (all admissions, for pooled-standard assembly).")

std_margins_eol = .build_margins(cohort[outcome_category %in% eol_categories_raw])
fwrite(std_margins_eol, shareable_path(
  paste0("aim2a3_standardization_margins_eol_", site_lowercase, "_", today, ".csv")))
message("  Saved aim2a3_standardization_margins_eol (end of life rows, for the eol standard population).")
rm(std_margins_eol)


# 3. Analytic dataset summary (shareable) --------------------------------------


message("\n--- Building shareable analytic summary ---")

# Covariate distributions by hospital unit — no patient-level data leaves the site.

analytic_summary = cohort[, .(
  n                     = .N,
  age_mean              = round(mean(age_at_admission, na.rm = TRUE), 1),
  age_sd                = round(sd(age_at_admission, na.rm = TRUE), 1),
  pct_female            = round(100 * mean(tolower(sex_category) == "female", na.rm = TRUE), 1),
  pct_icu               = round(100 * mean(icu_flag, na.rm = TRUE), 1),
  pct_lst_any           = round(100 * mean(lst_any_flag, na.rm = TRUE), 1),
  pct_icu_early         = round(100 * mean(icu_early_flag, na.rm = TRUE), 1),
  pct_lst_early         = round(100 * mean(lst_early_flag, na.rm = TRUE), 1),
  elix_mean             = round(mean(elixhauser_vw, na.rm = TRUE), 1),
  elix_sd               = round(sd(elixhauser_vw, na.rm = TRUE), 1),
  apache2_mean          = if (apache2_available) round(mean(apache2_imputed, na.rm = TRUE), 1) else NA_real_,
  apache2_sd            = if (apache2_available) round(sd(apache2_imputed,   na.rm = TRUE), 1) else NA_real_,
  pct_apache2_missing   = round(100 * mean(apache2_missing, na.rm = TRUE), 1),
  # dx_category now holds the ACUITY grouping (post-crosswalk), so these report
  # the acuity categories most relevant to end-of-life case mix, not the old
  # 22-level body-system names.
  pct_cancer            = round(100 * mean(dx_category == "Cancer/malignancy", na.rm = TRUE), 1),
  pct_chronic_organ_fail = round(100 * mean(dx_category == "Chronic organ failure", na.rm = TRUE), 1),
  pct_dementia_neuro    = round(100 * mean(dx_category == "Dementia/neurodegenerative/neurodevelopmental", na.rm = TRUE), 1),
  pct_infection_sepsis  = round(100 * mean(dx_category == "Acute infection/sepsis", na.rm = TRUE), 1),
  pct_hospice           = round(100 * mean(outcome_category %in% c(
    "same_facility_hospice_death", "same_facility_hospice_alive",
    "external_hospice_discharge"), na.rm = TRUE), 1),
  pct_death_acute       = round(100 * mean(outcome_category == "in_hospital_death_acute",
                                           na.rm = TRUE), 1)
), by = hospital_unit]

analytic_summary[, site := site_lowercase]

fwrite(
  analytic_summary,
  shareable_path(
    paste0("aim2a3_analytic_summary_", site_lowercase, "_", today, ".csv"))
)
message("  Saved aim2a3_analytic_summary.")


# =============================================================================
# 3b. Analysis 2A.1: Unadjusted hospital-level variation
# =============================================================================

message("\n--- Analysis 2A.1: Unadjusted variation ---")

hosp_rates = cohort[, .(
  n = .N,
  mortality_rate    = round(100 * mean(outcome_category == "in_hospital_death_acute"), 2),
  hospice_sf_rate   = round(100 * mean(outcome_category %in% c(
    "same_facility_hospice_death", "same_facility_hospice_alive")), 2),
  hospice_ext_rate  = round(100 * mean(outcome_category == "external_hospice_discharge"), 2),
  hospice_all_rate  = round(100 * mean(outcome_category %in% c(
    "same_facility_hospice_death", "same_facility_hospice_alive",
    "external_hospice_discharge")), 2),
  transfer_rate     = round(100 * mean(outcome_category == "transfer_acute_facility"), 2),
  has_hospice_unit  = as.integer(any(has_hospice_unit == 1L))
), by = hospital_unit]
hosp_rates[, site := site_lowercase]

n_hospice_unit_hosp = sum(hosp_rates$has_hospice_unit == 1L)
message(sprintf("  Hospitals with inpatient hospice unit: %d of %d",
                n_hospice_unit_hosp, nrow(hosp_rates)))

# variation statistics — two versions per §2.1:
#   (a) all hospitals (includes structural zeros at sites without hospice unit)
#   (b) restricted to hospitals with an inpatient hospice unit
# FIXED 2026.10.19 (UChicago, first single hospital site): this used to
# return NULL below 2 hospitals. At a site with ONE hospital_id both calls
# returned NULL, rbindlist() produced a zero column table, and the summary
# print four lines down crashed the whole stage on `version` (resolving to
# base R's version object). Return a TYPED empty table instead: the print
# loop runs zero times, fwrite writes a header only CSV the compile layer
# can read, and between hospital variation at a one hospital site correctly
# comes from the pooled coordinating center model, not from this file.
compute_var_stats = function(dt, version_label) {
  if (nrow(dt) < 2) {
    message(sprintf("  %s: %d hospital unit(s) — between hospital variation not defined; header only output.",
                    version_label, nrow(dt)))
    return(data.table(site = character(), version = character(),
                      n_hospitals = integer(), measure = character(),
                      iqr = numeric(), cv = numeric()))
  }
  data.table(
    site = site_lowercase, version = version_label, n_hospitals = nrow(dt),
    measure = c("mortality_rate", "hospice_sf_rate", "hospice_all_rate", "transfer_rate"),
    iqr = c(round(IQR(dt$mortality_rate), 2),
            round(IQR(dt$hospice_sf_rate), 2),
            round(IQR(dt$hospice_all_rate), 2),
            round(IQR(dt$transfer_rate), 2)),
    cv  = c(round(100 * sd(dt$mortality_rate) / max(mean(dt$mortality_rate), 0.01), 1),
            round(100 * sd(dt$hospice_sf_rate) / max(mean(dt$hospice_sf_rate), 0.01), 1),
            round(100 * sd(dt$hospice_all_rate) / max(mean(dt$hospice_all_rate), 0.01), 1),
            round(100 * sd(dt$transfer_rate) / max(mean(dt$transfer_rate), 0.01), 1))
  )
}

var_stats = rbindlist(list(
  compute_var_stats(hosp_rates, "all_hospitals"),
  compute_var_stats(hosp_rates[has_hospice_unit == 1L], "hospice_unit_only")
))

message("  Variation summary (all hospitals):")
for (i in seq_len(nrow(var_stats[version == "all_hospitals"]))) {
  row = var_stats[version == "all_hospitals"][i]
  message(sprintf("    %-25s IQR = %5.2f pp | CV = %5.1f%%",
                  row$measure, row$iqr, row$cv))
}

fwrite(hosp_rates, shareable_path(
  paste0("aim2a_hospital_rates_", site_lowercase, "_", today, ".csv")))
fwrite(var_stats, shareable_path(
  paste0("aim2a_variation_stats_", site_lowercase, "_", today, ".csv")))

# caterpillar + bivariate scatter plots
if (nrow(hosp_rates) >= 5) {
  hosp_rates[, rank := rank(hospice_sf_rate)]

  p_cat = ggplot(hosp_rates, aes(x = rank, y = hospice_sf_rate)) +
    geom_point(size = 2.5, color = pal_outcome[["Same-facility hospice"]]) +
    geom_hline(yintercept = median(hosp_rates$hospice_sf_rate),
               linetype = "dashed", color = "grey50") +
    labs(x = "Hospital rank", y = "Same-facility hospice rate (%)",
         title = "Unadjusted variation: same-facility hospice (all hospitals)",
         subtitle = paste0(toupper(site_lowercase), " | ", nrow(hosp_rates), " hospitals")) +
    theme_clif()

  # hospice-unit-only version (isolates practice from infrastructure variation)
  hu = hosp_rates[has_hospice_unit == 1L]
  if (nrow(hu) >= 3) {
    hu[, rank := rank(hospice_sf_rate)]
    p_cat_hu = ggplot(hu, aes(x = rank, y = hospice_sf_rate)) +
      geom_point(size = 2.5, color = pal_outcome[["Same-facility hospice"]]) +
      geom_hline(yintercept = median(hu$hospice_sf_rate),
                 linetype = "dashed", color = "grey50") +
      labs(x = "Hospital rank", y = "Same-facility hospice rate (%)",
           title = "Same-facility hospice (hospice-unit hospitals only)",
           subtitle = paste0(toupper(site_lowercase), " | ", nrow(hu), " hospitals")) +
      theme_clif()
  } else {
    p_cat_hu = NULL
  }

  p_scatter = ggplot(hosp_rates, aes(x = hospice_sf_rate, y = mortality_rate)) +
    geom_point(size = 2.5, color = pal_outcome[["Same-facility hospice"]]) +
    geom_smooth(method = "lm", se = TRUE, color = pal_outcome[["In-hospital death"]], linewidth = 0.6) +
    labs(x = "Same-facility hospice rate (%)", y = "In-hospital mortality rate (%)",
         title = "Hospital-level: hospice rate vs mortality rate",
         subtitle = toupper(site_lowercase)) +
    theme_clif()

  p_scatter2 = ggplot(hosp_rates, aes(x = transfer_rate, y = mortality_rate)) +
    geom_point(size = 2.5, color = pal_outcome[["Routine discharge"]]) +
    geom_smooth(method = "lm", se = TRUE, color = pal_outcome[["In-hospital death"]], linewidth = 0.6) +
    labs(x = "Transfer-out rate (%)", y = "In-hospital mortality rate (%)",
         title = "Hospital-level: transfer rate vs mortality rate",
         subtitle = toupper(site_lowercase)) +
    theme_clif()

  # scatter: external vs same-facility hospice rate (§2.1)
  p_scatter3 = ggplot(hosp_rates, aes(x = hospice_sf_rate, y = hospice_ext_rate)) +
    geom_point(size = 2.5, color = pal_outcome[["External hospice"]]) +
    geom_smooth(method = "lm", se = TRUE, color = pal_outcome[["In-hospital death"]], linewidth = 0.6) +
    labs(x = "Same-facility hospice rate (%)", y = "External hospice rate (%)",
         title = "Hospital-level: same-facility vs external hospice",
         subtitle = toupper(site_lowercase)) +
    theme_clif()

  top_row = if (!is.null(p_cat_hu)) (p_cat | p_cat_hu) else p_cat
  p_variation = top_row / (p_scatter | p_scatter2) / p_scatter3
  ggsave(shareable_path(
    paste0("aim2a_unadjusted_variation_", site_lowercase, "_", today, ".pdf")),
    plot = p_variation, width = 12, height = 14, dpi = 300)
  message("  Saved aim2a_unadjusted_variation plots.")
  suppressWarnings(rm(p_cat, p_cat_hu, p_scatter, p_scatter2, p_scatter3, p_variation, hu))
}


# =============================================================================
# 3c. Analysis 2A.2: Hospice-to-mortality ratios
# =============================================================================

message("\n--- Analysis 2A.2: Hospice-to-mortality ratios ---")

# hospital-level: (same-facility hospice) / (acute in-hospital deaths)
hosp_rates[, hospice_to_mort_ratio := fifelse(
  mortality_rate > 0,
  round(hospice_sf_rate / mortality_rate, 3),
  NA_real_
)]

# also compute on count basis (more intuitive)
hosp_counts = cohort[, .(
  n_death_acute = sum(outcome_category == "in_hospital_death_acute"),
  n_hospice_sf  = sum(outcome_category %in% c("same_facility_hospice_death",
                                              "same_facility_hospice_alive"))
), by = hospital_unit]
hosp_counts[, hospice_to_mort_count := round(n_hospice_sf / pmax(n_death_acute, 1), 3)]

hosp_rates = merge(hosp_rates, hosp_counts[, .(hospital_unit, n_death_acute,
                                               n_hospice_sf, hospice_to_mort_count)],
                   by = "hospital_unit", all.x = TRUE)

message(sprintf("  Hospice-to-mortality ratio: median = %.3f, IQR = %.3f–%.3f",
                median(hosp_rates$hospice_to_mort_count, na.rm = TRUE),
                quantile(hosp_rates$hospice_to_mort_count, 0.25, na.rm = TRUE),
                quantile(hosp_rates$hospice_to_mort_count, 0.75, na.rm = TRUE)))

# 30-day mortality version at death-linked sites
is_death_linked = "death_linked_site" %in% names(cohort) &&
  any(cohort$death_linked_site == 1L, na.rm = TRUE)

if (is_death_linked) {
  message("  Computing hospice-to-30-day-mortality ratio (death-linked site)...")

  hosp_30d = cohort[, .(
    n_died_30d = sum(!is.na(days_discharge_to_death) & days_discharge_to_death <= 30 |
                       outcome_category == "in_hospital_death_acute" |
                       outcome_category == "same_facility_hospice_death"),
    n_hospice_all = sum(outcome_category %in% c("same_facility_hospice_death",
                                                "same_facility_hospice_alive",
                                                "external_hospice_discharge"))
  ), by = hospital_unit]

  hosp_30d[, hospice_to_30d_mort := round(n_hospice_all / pmax(n_died_30d, 1), 3)]

  hosp_rates = merge(hosp_rates, hosp_30d[, .(hospital_unit, hospice_to_30d_mort)],
                     by = "hospital_unit", all.x = TRUE)

  message(sprintf("  30-day ratio: median = %.3f (national benchmark = 0.27)",
                  median(hosp_rates$hospice_to_30d_mort, na.rm = TRUE)))
}

# save
fwrite(hosp_rates, shareable_path(
  paste0("aim2a_hospital_ratios_", site_lowercase, "_", today, ".csv")))
message("  Saved aim2a_hospital_ratios.")

rm(hosp_rates, hosp_counts, var_stats)
suppressWarnings(rm(hosp_30d))



# 4. Model fitting: binary GLMMs per outcome category --------------------------


message("\n--- Fitting binary GLMMs ---")

# Outcome categories to model (each vs. routine_discharge as reference).
# Encounters with outcome = routine_discharge are always in the reference group.
# Each model is: logit(P(outcome_k)) = Xβ + (1|hospital_unit)

# CHANGED 2026.10.14 (research in progress feedback). Three analytic cohorts,
# all cut from model_base (one admission per patient) further down:
#
#   eol (PRIMARY)  admissions ending in death or hospice. Outcomes are modeled
#                  against in_hospital_death_acute, so each hospital effect is
#                  the adjusted hospice to death ratio on the log odds scale,
#                  the model version of the 2A.2 ratio.
#   full           all dispositions, one admission per patient. INCIDENCE fits
#                  only: any hospice and in hospital death, each against
#                  routine_discharge, so the paper can still say how often
#                  admissions end in hospice or death, not only how deaths are
#                  routed. These are the only fits that still subsample.
#   terminal       death linked sites only; membership by vital status (died
#                  in hospital or within terminal_death_window_days of
#                  discharge, any disposition) rather than by disposition.
#                  Same fits as eol plus a fourth category,
#                  died_30d_no_hospice. The concordance of its hospital
#                  effects with the eol fits is the check on whether
#                  disposition defined membership distorts the primary
#                  results.
#
# Same facility hospice: the death / alive split that 01 assigns is collapsed
# in model_base below, so same_facility_hospice is a single level here.
eol_ref            = "in_hospital_death_acute"
eol_outcomes       = c("same_facility_hospice", "external_hospice_discharge",
                       "any_hospice_disposition")
incidence_ref      = "routine_discharge"
incidence_outcomes = c("any_hospice_disposition", "in_hospital_death_acute")
terminal_outcomes  = c(eol_outcomes, "died_30d_no_hospice")

# protocol §3.4.4: same-facility hospice is only meaningful where the pathway
# exists and is captured, so that outcome is analysed CONDITIONAL on class (a)
# hospitals. The network-wide outcome, defined for every hospital, is any
# documented hospice disposition (same-facility OR external).
sf_outcomes = "same_facility_hospice"

# any_hospice_disposition is the network wide outcome of every cohort here, so
# it cannot be optional any more; derive it if 01 did not supply it
if (!"any_hospice_disposition" %in% names(cohort)) {
  message("  NOTE: any_hospice_disposition not found in the cohort -- derived here as any same-facility or external hospice disposition (re-run 01 to build it, protocol 3.4.4).")
  cohort[, any_hospice_disposition := as.integer(outcome_category %in% c(
    "same_facility_hospice_death", "same_facility_hospice_alive", "external_hospice_discharge"))]
}
if (!"pathway_class" %in% names(cohort)) {
  message("  NOTE: pathway_class not found -- the same-facility outcome will NOT be restricted")
  message("        to class (a) hospitals. Re-run 01 to build the §3.4.4 classification.")
}

# analysis dataset: routine_discharge + one outcome at a time
# Exclude encounters missing a genuinely-required covariate. NOTE: filter only
# on raw/required covariates — NOT on the imputed-and-centered columns
# (elix_centered, apache2_centered). Those are never NA under normal partial
# missingness (they're median/zero-imputed), so conditioning on them is a no-op
# in the normal case and silently drops the ENTIRE cohort when a whole source
# column is missing (e.g. APACHE absent at a site). That was the 0-row bug.
model_base = cohort[
  !is.na(age_centered) &
    !is.na(sex_female) &
    !is.na(race_eth) &
    !is.na(dx_category)
]

# CHANGED 2026.10.14: collapse the two same facility hospice categories into
# one, in the MODELING sample only. cohort keeps the split from 01, so 2A.1,
# 2A.2, and the analytic summary above still report death and alive
# separately. Every model below (binary GLMMs, multinomial, era refits) reads
# outcome_category from model_base, so all of them see the same five
# dispositions against routine discharge.
n_sf_split = model_base[outcome_category %in% c("same_facility_hospice_death",
                                                "same_facility_hospice_alive"), .N]
model_base[, outcome_category := fifelse(
  outcome_category %in% c("same_facility_hospice_death", "same_facility_hospice_alive"),
  "same_facility_hospice", as.character(outcome_category))]
message(sprintf("  Same-facility hospice: death and alive collapsed to one outcome (%s encounters).",
                format(n_sf_split, big.mark = ",")))
rm(n_sf_split)

# ensure factor levels
model_base[, sex_female := factor(sex_female)]
model_base[, race_eth := relevel(factor(race_eth), ref = "White non-Hispanic")]
# diagnosis-category reference = most frequent body system present in the model
# data. Robust to the CCSR scheme (and any future regrouping) — avoids hardcoding
# a level like "Other" that may not exist and would crash relevel().
dx_ref = model_base[, .N, by = dx_category][order(-N)][1L, dx_category]
model_base[, dx_category := relevel(factor(dx_category), ref = dx_ref)]
message(sprintf("  dx_category reference level: %s", dx_ref))

# pandemic period fixed effect (§3.3, §3.6.2)
if (!"pandemic_period" %in% names(model_base)) {
  if (exists("assign_pandemic_period")) {
    model_base[, pandemic_period := assign_pandemic_period(admission_dttm)]
  } else {
    model_base[, pandemic_period := "Pre-pandemic"]
  }
}
# FIXED 2026.10.20 (Emory): levels come from the eras OBSERVED here, in
# calendar order, and the reference is the earliest observed era. The old
# relevel(ref = "Pre-pandemic") stopped the stage at any site whose data
# begin after 2019. Unrecognized labels (none expected) sort after the three
# protocol eras.
era_order    = c("Pre-pandemic", "Acute pandemic", "Post-pandemic")
era_seen     = unique(as.character(model_base$pandemic_period))
era_seen     = era_seen[!is.na(era_seen)]
era_observed = c(intersect(era_order, era_seen), sort(setdiff(era_seen, era_order)))
model_base[, pandemic_period := factor(as.character(pandemic_period), levels = era_observed)]
era_reference = era_observed[1]
era_in_model  = length(era_observed) >= 2L
message(sprintf("  Pandemic eras observed: %s | reference: %s%s",
                paste(era_observed, collapse = ", "), era_reference,
                if (era_in_model) "" else " -> single era site: era term NOT in the model; era refits skipped"))
run_config = rbind(run_config, data.table(
  setting = c("era_levels_observed", "era_reference", "era_term_in_model"),
  value   = c(paste(era_observed, collapse = "; "), era_reference, toupper(format(era_in_model))),
  site           = site_lowercase,
  script_version = run_config$script_version[1]), use.names = TRUE, fill = TRUE)

n_model = nrow(model_base)
message(sprintf("  Analytic sample: %s encounters (%s excluded for missing covariates).",
                format(n_model, big.mark = ","),
                format(n_cohort - n_model, big.mark = ",")))

# health_system must be added BEFORE the reference group is split off, so
# ref_group carries the SAME columns as the outcome groups it is later rbind()ed
# with inside the GLMM (a data.table := after the split would not reach the
# already-separated ref_group). Default: one site == one health system; sites
# spanning several systems should populate health_system upstream (in 01).
if (!"health_system" %in% names(model_base))
  model_base[, health_system := site_lowercase]


## analytic unit: one admission per patient (protocol 2A.3) --------------------
# CHANGED 2026.10.13. The analytic unit is ONE admission per patient: the
# LAST admission in the study period. Each patient then appears once, so
# within patient dependence is removed by construction and no model carries a
# patient term. (The patient level random effect, formerly attempted where
# estimable and then kept as a sensitivity, was removed outright on
# 2026.10.14: almost no CLIF site is small enough for lme4 to fit it.)
# Everything fit from model_base below inherits the sample: the full cohort,
# LST and ICU fits, the clinical course models, the multinomial, and the era
# refits.
#
# "Last" is the latest admission_dttm among the patient's admissions in the
# analytic sample, with ties broken by hospitalization_id so every site
# resolves them the same way. The patient is attributed to the hospital of
# that admission. An encounter with no patient_id cannot be linked to any
# other, so it is kept as its own patient and the count is reported.
#
# Scope: this restricts the MODEL sample only. The unadjusted 2A.1 and 2A.2
# tables, the analytic summary, and the standardization margins above are
# built from cohort and still describe ALL admissions.
if (!all(c("patient_id", "admission_dttm") %in% names(model_base)))
  stop("patient_id and admission_dttm are both required to select one ",
       "admission per patient.")

setorder(model_base, patient_id, admission_dttm, hospitalization_id)
pid_ok     = model_base[, !is.na(patient_id) & nzchar(patient_id)]
model_base = model_base[!duplicated(patient_id, fromLast = TRUE) | !pid_ok]

n_patients = nrow(model_base)
pct_repeat = 100 * (n_model - n_patients) / n_model

message(sprintf("  Analytic sample: last admission per patient, %s patients from %s encounters (%.1f%% of encounters were earlier admissions of a repeat patient).",
                format(n_patients, big.mark = ","),
                format(n_model, big.mark = ","), pct_repeat))
if (any(!pid_ok))
  message(sprintf("  NOTE: %s encounters have no patient_id; each is kept as its own patient.",
                  format(sum(!pid_ok), big.mark = ",")))
rm(pid_ok)


## analytic cohorts: end of life (primary), full (incidence), terminal ----------
# CHANGED 2026.10.14. All three are cut from model_base (one admission per
# patient, same facility collapsed), so they share the covariate coding.

# end of life cohort: the admission ended in death or hospice
model_eol = model_base[outcome_category %in% c(eol_ref, "same_facility_hospice",
                                               "external_hospice_discharge")]
n_eol = nrow(model_eol)
message(sprintf("  End of life cohort (PRIMARY): %s patients, %.1f%% of the analytic sample (%s in hospital deaths, %s same facility hospice, %s external hospice).",
                format(n_eol, big.mark = ","), 100 * n_eol / n_patients,
                format(model_eol[outcome_category == eol_ref, .N], big.mark = ","),
                format(model_eol[outcome_category == "same_facility_hospice", .N], big.mark = ","),
                format(model_eol[outcome_category == "external_hospice_discharge", .N], big.mark = ",")))

# opportunity sensitivity: in hospital deaths within early_death_window_hours
# of admission had no realistic hospice option (arrests, trauma, deaths soon
# after arrival), and hospitals with more of them would look worse on hospice
# vs death whatever their practice. Hospice rows are untouched.
if (!"los_days" %in% names(model_eol))
  model_eol[, los_days := as.numeric(difftime(discharge_dttm, admission_dttm, units = "days"))]
model_eol[, early_death := as.integer(outcome_category == eol_ref &
                                        !is.na(los_days) &
                                        los_days < early_death_window_hours / 24)]
message(sprintf("  Early deaths (in hospital death within %d h of admission): %s, %.1f%% of in hospital deaths; dropped in the eol_excl_death_%dh sensitivity.",
                early_death_window_hours,
                format(sum(model_eol$early_death), big.mark = ","),
                100 * sum(model_eol$early_death) / max(model_eol[outcome_category == eol_ref, .N], 1L),
                early_death_window_hours))

# terminal cohort (death linked sites only): membership by vital status, so
# the hospital's disposition decision does not decide who is in the cohort.
# is_death_linked was set in 2A.2. Hospice rows keep their category; every
# other disposition that ended in death within the window (routine, transfer,
# AMA) becomes died_30d_no_hospice. Transfers are lumped in because the
# receiving hospital's care is not observed here; their count is reported so
# the lump is visible.
if (is_death_linked && "days_discharge_to_death" %in% names(model_base)) {
  died_in_window = model_base[, !is.na(days_discharge_to_death) &
                                days_discharge_to_death >= 0 &
                                days_discharge_to_death <= terminal_death_window_days]
  model_term = model_base[outcome_category == eol_ref | died_in_window]
  n_term_by_disp = model_term[, .N, by = outcome_category]
  model_term[!outcome_category %in% c(eol_ref, "same_facility_hospice",
                                      "external_hospice_discharge"),
             outcome_category := "died_30d_no_hospice"]
  message(sprintf("  Terminal cohort (death linked site, death within %d days of discharge): %s patients.",
                  terminal_death_window_days, format(nrow(model_term), big.mark = ",")))
  for (i in seq_len(nrow(n_term_by_disp)))
    message(sprintf("    %-30s %s%s", n_term_by_disp$outcome_category[i],
                    format(n_term_by_disp$N[i], big.mark = ","),
                    if (!n_term_by_disp$outcome_category[i] %in% c(eol_ref, "same_facility_hospice", "external_hospice_discharge"))
                      "  -> died_30d_no_hospice" else ""))
  rm(died_in_window, n_term_by_disp)
} else {
  model_term = NULL
  message("  Terminal cohort: not a death linked site (or no days_discharge_to_death) -> skipped.")
}

# reference group for the INCIDENCE fits (routine discharge). The end of life
# and terminal fits use in hospital death; fit_family() sets that per family.
ref_group = model_base[outcome_category == incidence_ref]
n_ref = nrow(ref_group)
message(sprintf("  Reference group, incidence fits (%s): %s encounters.",
                incidence_ref, format(n_ref, big.mark = ",")))

# build model formulas (calendar-era term included).
# base_rhs is the single source of truth: formula_fixed, formula_glmm, the
# multinomial mn_formula, and the per-era refit all derive from it, so dropping
# APACHE here drops it everywhere. When APACHE is unavailable, apache2_centered
# is constant 0 and apache2_missing is constant 1 (collinear with the intercept);
# feeding either to glm/glmer/mblogit would error or drop rank, so exclude both.
#
# PRIMARY model (protocol 2A.3): case-mix adjustment uses ONLY covariates fixed
# at or shortly after admission. Post-admission variables (ICU admission, organ
# support, later severity, length of stay) are EXCLUDED here because they may be
# consequences of hospital practice or disease evolution -- conditioning on them
# would adjust away part of the very variation Aim 2A is trying to measure.
# They move to the secondary clinical-course model below.
# NOTE: icu_flag was previously in this vector. Removing it also resolves a
# latent rank deficiency in the "icu_only" subgroup, where icu_flag is constant.
# the two v7.0 baseline adjusters enter only if they vary in the analytic sample
# (a constant factor would be collinear with the intercept and drop rank)
extra_base = c(
  if ("preexisting_code_status" %in% names(model_base) &&
      uniqueN(model_base$preexisting_code_status) >= 2) "preexisting_code_status",
  if ("admission_source" %in% names(model_base) &&
      uniqueN(model_base$admission_source) >= 2) "admission_source"
)

base_rhs = paste(
  c("age_centered", "sex_female", "race_eth", "dx_category",
    "elix_centered",
    if (apache2_available) c("apache2_centered", "apache2_missing"),
    extra_base,
    if (era_in_model) "pandemic_period"),   # 2026.10.20: constant at one era sites
  collapse = " + "
)
message(sprintf("  Primary model RHS:         %s", base_rhs))
if (length(extra_base))
  message(sprintf("  v7.0 baseline adjusters included: %s", paste(extra_base, collapse = ", ")))


## random-effect structure (protocol 2A.3) -------------------------------------
# The hospital random intercept is always included. One further term is
# protocol-specified but only conditionally estimable:
#
#   (1 | health_system) -- hospitals nested in health systems. A single CLIF site
#     is normally ONE health system, so this variance is not identifiable from
#     site-local data; it is recovered in the pooled coordinating-center model.
#     Included here only if a site genuinely spans several systems.
#
# There is no patient term. The analytic sample has one admission per
# patient, so (1 | patient_id) would be degenerate; the patient level random
# effect and its estimability test were removed 2026.10.14.

# health_system was populated above, before the reference-group split.
n_systems = model_base[, uniqueN(health_system)]

use_system_re = n_systems >= 3L

message(sprintf("  Health systems in local data: %d%s", n_systems,
                if (use_system_re) "" else " -> system random effect NOT estimable locally (pooled analysis only)"))

# ADDED 2026.10.13: the analytic unit travels with the run configuration, so
# the coordinating center can see which admission rule a site's rows were fit
# under. Same file as section 0, rewritten with more rows. EXTENDED
# 2026.10.14 with the cohort definitions and sizes, so a missing terminal_30d
# family reads as "not a death linked site" rather than as a failure.
run_config = rbind(run_config, data.table(
  setting = c("primary_admission_rule", "n_encounters_all_admissions",
              "n_patients", "pct_repeat_admissions",
              "primary_cohort", "eol_reference_category", "n_eol",
              "early_death_window_hours", "n_early_deaths",
              "terminal_cohort", "terminal_death_window_days", "n_terminal",
              "elixhauser_poa_restricted", "poa_column",
              "pct_dx_poa_present", "pct_dx_poa_not_present", "pct_dx_poa_not_stated"),
  value   = c("last_per_patient",
              format(n_model, scientific = FALSE),
              format(n_patients, scientific = FALSE),
              format(round(pct_repeat, 2), nsmall = 2),
              "end_of_life", eol_ref, format(n_eol, scientific = FALSE),
              format(early_death_window_hours, scientific = FALSE),
              format(sum(model_eol$early_death), scientific = FALSE),
              if (is.null(model_term)) "not_death_linked" else "run",
              format(terminal_death_window_days, scientific = FALSE),
              if (is.null(model_term)) "0" else format(nrow(model_term), scientific = FALSE),
              # ADDED 2026.10.14: how complete the POA restriction of the
              # comorbidity index was at this site (section 2a)
              format(elix_poa_restricted),
              if (is.na(poa_col)) "none" else poa_col,
              format(round(100 * poa_share[["present"]], 2), nsmall = 2),
              format(round(100 * poa_share[["not_present"]], 2), nsmall = 2),
              format(round(100 * poa_share[["not_stated"]], 2), nsmall = 2)),
  site           = site_lowercase,
  script_version = run_config$script_version[1]), use.names = TRUE)
fwrite(run_config, shareable_path(
  paste0("aim2a3_run_config_", site_lowercase, "_", today, ".csv")))

re_rhs = paste(c("(1 | hospital_unit)",
                 if (use_system_re) "(1 | health_system)"), collapse = " + ")

message(sprintf("  Random effects:            %s", re_rhs))

formula_fixed = as.formula(paste("y ~", base_rhs))
formula_glmm  = as.formula(paste("y ~", base_rhs, "+", re_rhs))

# SECONDARY clinical-course model (protocol 2A.3): adds post-admission clinical
# course back onto the admission-fixed base. The protocol interprets this as a
# DIFFERENT ESTIMAND and requires it be reported separately, so it is fit on the
# full cohort only and written to its own output files -- it never mixes into the
# primary coefficient/ICC tables.
# Change 3 (2026.08.24, folded 2026.08.27): the exposures are now the EARLY
# windowed flags built in section 2g (first course_window_hours from
# admission), closing the protocol caveat about whole stay exposures. The
# whole stay flags remain the SUBGROUP definitions (lst_any / icu_only),
# unchanged.
course_rhs = paste(c(base_rhs, "icu_early_flag", "lst_early_flag"), collapse = " + ")
message(sprintf("  Clinical-course model RHS: %s", course_rhs))

formula_fixed_course = as.formula(paste("y ~", course_rhs))
formula_glmm_course  = as.formula(paste("y ~", course_rhs, "+", re_rhs))

# ---- Hospital variation profiling (protocol 2A.3) ----------------------------
# The protocol replaced the ICC with a profiling suite. From the fitted GLMM we
# report, per disposition: between-hospital variance with uncertainty; risk-
# standardized hospital probabilities (shrunken) with uncertainty; the contrast
# in standardized probability between the 10th and 90th percentiles of the
# hospital random-effect distribution; the median odds ratio; the reliability of
# each hospital's estimate; and rank uncertainty rather than definitive ranks.
# ICC is retained as a secondary summary only (a monotone function of sigma^2).
#
# Efficiency: the standardized probability is monotone increasing in the hospital
# effect u, so an interval for u maps directly onto an interval for the
# probability; no nested simulation over patients is needed. Ranks depend only on
# u, so rank uncertainty is simulated on the hospital effects alone.
# ---- Interval for sigma^2 -----------------------------------------------------
# confint(fit, method = "profile") holds theta at each of up to 100 grid points
# and optimizes the ENTIRE fixed effect vector again at every one. At a site with
# a few million encounters that is 40 to 100 full refits of the model for a
# single number, and it is what stalled the first large-site test run.
#
# The interval below profiles the SAME likelihood in one dimension, using the
# nAGQ = 0 deviance function in which the fixed effects are profiled out inside
# PIRLS. Each grid point is then one PIRLS solve rather than one refit, and the
# bounds come from root finding on the chi-square(1) cutoff. Cost drops from tens
# of thousands of PIRLS solves to a few dozen.
compute_sigma_ci = TRUE   # pre-existing switch, unchanged. FALSE writes NA
# bounds and skips the interval entirely.

# one-dimensional profile. The nAGQ = 0 deviance function takes theta alone, so
# no inner optimization over the fixed effects is needed.
sigma_ci_fast = function(fit, model_dt, level = 0.95) {
  # change 4: the deviance surface must use the same weights as the fit
  # (the standalone validation script copy stays intentionally unweighted)
  dd = if ("row_weight" %in% names(model_dt)) {
    glmer(formula(fit), data = model_dt, family = binomial(link = "logit"),
          weights = row_weight, nAGQ = 0, devFunOnly = TRUE)
  } else {
    glmer(formula(fit), data = model_dt, family = binomial(link = "logit"),
          nAGQ = 0, devFunOnly = TRUE)
  }

  th = getME(fit, "theta")
  k  = which(names(th) == "hospital_unit.(Intercept)")
  if (length(k) != 1L) k = 1L

  # deviance with the hospital SD held at s. When the site also fits a system
  # random effect, that variance component is optimized out (cheap: one
  # nuisance parameter).
  pdev = if (length(th) == 1L) {
    function(s) dd(s)
  } else {
    function(s) {
      o = optim(as.numeric(th[-k]),
                function(v) { pp = numeric(length(th)); pp[k] = s; pp[-k] = v; dd(pp) },
                method = "L-BFGS-B", lower = 0,
                control = list(maxit = 100, factr = 1e9))
      o$value
    }
  }

  # the nAGQ = 0 surface has its own minimum, so re-locate it before cutting
  s_hat = optimize(pdev, interval = c(1e-4, max(10 * as.numeric(th[k]), 1)))$minimum
  cut   = pdev(s_hat) + qchisq(level, 1)
  f     = function(s) pdev(s) - cut

  lo = if (f(1e-6) <= 0) 0 else
    tryCatch(uniroot(f, c(1e-6, s_hat), tol = 1e-4)$root, error = \(e) NA_real_)

  # walk the upper bound out until the deviance clears the cutoff
  s_up = NA_real_
  for (m in c(2, 4, 8, 16, 32)) if (f(m * s_hat) > 0) { s_up = m * s_hat; break }
  hi = if (is.na(s_up)) NA_real_ else
    tryCatch(uniroot(f, c(s_hat, s_up), tol = 1e-4)$root, error = \(e) NA_real_)

  c(lo^2, hi^2)
}

# lp_offset (change 1): when the fit used control subsampling at fraction f,
# the fitted intercept is on the sampled scale (sampled = true - log(f), with
# log(f) negative). Passing lp_offset = log(f) restores the population scale
# linear predictor inside the standardization, so RS probabilities and the
# p10/p90 contrast are population scale. sigma^2, MOR, ICC, reliability, and
# ranks depend only on the random effects and need no correction.
# CHANGED 2026.10.14: the standardization frame, its provenance, and the
# reference category are ARGUMENTS (std_frame, std_source, std_version,
# ref_label) rather than globals, because the end of life fits standardize on
# a different population from the incidence fits. std_frame = NULL means
# site local standardization on model_dt.
hospital_variation_profile = function(fit, model_dt, outcome_label, subgroup_label,
                                      n_outcome, n_hospitals,
                                      B = 2000, std_n = 10000, lp_offset = 0,
                                      std_frame = NULL, std_source = "site_local",
                                      std_version = "site_local",
                                      ref_label = NA_character_) {

  vc     = as.data.table(VarCorr(fit))
  sigma2 = vc[grp == "hospital_unit", vcov][1]

  # a singular / zero-variance fit has no between-hospital variation to profile
  if (!is.finite(sigma2) || sigma2 <= 1e-8) {
    message("      Between-hospital variance is ~0 (singular fit); skipping profiling.")
    return(NULL)
  }
  sigma = sqrt(sigma2)

  set.seed(20260713)   # reproducible standardization sample and rank draws

  # between-hospital variance, with a profile-likelihood interval (SD scale)
  sig_ci = c(NA_real_, NA_real_)
  if (isTRUE(compute_sigma_ci)) {
    t_ci   = Sys.time()
    sig_ci = tryCatch(sigma_ci_fast(fit, model_dt), error = function(e) {
      message("      sigma^2 CI unavailable: ", conditionMessage(e))
      c(NA_real_, NA_real_)
    })
    message(sprintf("      sigma^2 CI in %.1f min.",
                    as.numeric(difftime(Sys.time(), t_ci, units = "mins"))))
  }

  # median odds ratio (Larsen & Merlo)
  mor = exp(sqrt(2 * sigma2) * qnorm(0.75))

  # empirical-Bayes hospital effects and posterior SDs
  re   = ranef(fit, condVar = TRUE)$hospital_unit
  u    = re[, 1]
  pv   = as.numeric(attr(re, "postVar")[1, 1, ])
  sd_u = sqrt(pv)

  # reliability: share of a hospital's apparent variation that is signal.
  # For the normal-normal model, postVar = sigma2 * (1 - reliability).
  reliability = pmin(pmax(1 - pv / sigma2, 0), 1)

  # standardized patient population: linear predictor with the hospital effect
  # removed, so every hospital is evaluated on the SAME case mix
  # standardize on the NETWORK standard population if one was distributed
  # (comparable across sites); otherwise fall back to a site-local sample.
  # change 4: when rows are compressed patterns, draw patterns with prob
  # proportional to row_weight (with replacement) to recover the encounter
  # level case mix in the site local standardization sample
  std_dt = if (!is.null(std_frame)) std_frame
  else if ("row_weight" %in% names(model_dt) && any(model_dt$row_weight > 1L))
    model_dt[sample(.N, std_n, replace = TRUE, prob = row_weight)]
  else if (nrow(model_dt) > std_n) model_dt[sample(.N, std_n)] else model_dt
  xb     = predict(fit, newdata = std_dt, re.form = NA, type = "link") + lp_offset
  rsp    = function(uu) mean(plogis(xb + uu))

  rs_prob = vapply(u,               rsp, numeric(1))
  rs_lo   = vapply(u - 1.96 * sd_u, rsp, numeric(1))
  rs_hi   = vapply(u + 1.96 * sd_u, rsp, numeric(1))

  # rank uncertainty: draw from each hospital's posterior, rank within each draw
  draws    = vapply(seq_along(u), function(j) rnorm(B, u[j], sd_u[j]), numeric(B))
  if (is.null(dim(draws))) draws = matrix(draws, nrow = B)
  rk       = t(apply(draws, 1, function(r) rank(-r, ties.method = "average")))
  rank_med = apply(rk, 2, median)
  rank_lo  = apply(rk, 2, quantile, probs = 0.025, names = FALSE)
  rank_hi  = apply(rk, 2, quantile, probs = 0.975, names = FALSE)

  hosp_tbl = data.table(
    hospital_unit   = rownames(re),
    outcome         = outcome_label,
    subgroup        = subgroup_label,
    reference_category = ref_label,
    site            = site_lowercase,
    ranef_intercept = round(u, 4),
    ranef_sd        = round(sd_u, 4),
    rs_probability  = round(rs_prob, 5),
    rs_prob_lower   = round(rs_lo, 5),
    rs_prob_upper   = round(rs_hi, 5),
    rs_prob_se      = round((rs_hi - rs_lo) / (2 * 1.96), 5),  # for CC meta-analysis
    standardization = std_source,
    standard_pop_version = std_version,
    reliability     = round(reliability, 3),
    rank_median     = round(rank_med, 1),
    rank_lower      = round(rank_lo, 1),
    rank_upper      = round(rank_hi, 1),
    n_outcome       = n_outcome
  )

  # 10th vs 90th percentile of the hospital random-effect distribution
  p10 = rsp(qnorm(0.10, 0, sigma))
  p90 = rsp(qnorm(0.90, 0, sigma))

  summ = data.table(
    outcome               = outcome_label,
    subgroup              = subgroup_label,
    reference_category    = ref_label,
    site                  = site_lowercase,
    model_type            = "glmm",
    sigma2_hospital       = round(sigma2, 4),
    sigma2_lower          = round(sig_ci[1], 4),
    sigma2_upper          = round(sig_ci[2], 4),
    median_or             = round(mor, 3),
    icc                   = round(sigma2 / (sigma2 + pi^2 / 3), 4),  # secondary
    rs_prob_p10           = round(p10, 5),
    rs_prob_p90           = round(p90, 5),
    rs_prob_p90_minus_p10 = round(p90 - p10, 5),
    rs_prob_p90_over_p10  = round(p90 / max(p10, 1e-12), 3),
    median_reliability    = round(median(reliability), 3),
    n_outcome             = n_outcome,
    n_hospitals           = n_hospitals,
    n_standardization     = nrow(std_dt)
  )

  message(sprintf("      MOR = %.2f | sigma2 = %.4f | RS prob p10->p90: %.4f -> %.4f (diff %.4f) | median reliability = %.2f",
                  mor, sigma2, p10, p90, p90 - p10, median(reliability)))

  list(summary = summ, hospitals = hosp_tbl)
}


# ---- Hospital fixed intercept for the fixed effects fallbacks ----------------
# ADDED 2026.10.16. Below the 3 hospital floor the hospital random intercept
# is not estimable: two hospitals give one contrast, not a variance, and the
# mblogit inner solve is singular. The fallbacks used to drop hospital from
# the model entirely, which made the fallback ORs population averaged rather
# than within hospital (not comparable with the GLMM's conditional ORs at
# larger sites) and left their SEs blind to the clustering. Hospital now
# enters each fallback as a FIXED intercept: one dummy per hospital beyond
# the reference, which is the hospital with the most encounters in the
# analytic sample. That restores the within hospital interpretation and
# absorbs the intercept level clustering, so no separate SE correction is
# needed.
#
# The hospital contrasts travel as ordinary coefficient rows (term =
# hospital_fe<unit>, log odds vs the reference hospital) with
# hospital_fe_reference naming the reference unit, so the coordinating
# center can report them as site level hospital contrasts and keep them out
# of the covariate pooling by the hospital_fe prefix. A hospital with no
# events in a fit separates on its dummy (OR near 0, very large SE), exactly
# as a rare covariate level does; the covariate ORs are unaffected and the
# compile step drops the separated row as it already does for those levels.
# An analytic sample with a single hospital gets no dummy (it would alias
# the intercept). sigma^2, the MOR, the shrunken hospital profiles and the
# risk standardized composition remain GLMM only.
add_hospital_fe = function(dt) {
  counts = if ("row_weight" %in% names(dt)) dt[, .(N = sum(row_weight)), by = hospital_unit]
           else dt[, .N, by = hospital_unit]
  ref_h = counts[order(-N, hospital_unit)][1L, hospital_unit]
  dt[, hospital_fe := relevel(factor(hospital_unit), ref = ref_h)]
  ref_h
}


# function to fit one binary GLMM
# CHANGED 2026.10.14: ref_label names the reference category the globals
# ref_group / n_ref currently hold (set per family by fit_family()); it is
# carried into every output table as reference_category so a vs death row can
# never be confused with a vs routine row. std_frame / std_source /
# std_version are passed through to hospital_variation_profile().
fit_binary_glmm = function(outcome_label, subgroup_label = "eol",
                           f_fixed = formula_fixed, f_glmm = formula_glmm,
                           stage_label = "primary",
                           ref_label = eol_ref,
                           std_frame = NULL, std_source = "site_local",
                           std_version = "site_local") {

  # protocol §3.4.4: the PRIMARY same-facility analysis is CONDITIONAL on class (a)
  # hospitals (pathway available and reliably captured). At other hospitals the
  # absence of the pathway is a structural / data-capture fact, not a clinical
  # choice, so leaving them in would rank billing and infrastructure rather than
  # practice. any_hospice_disposition is defined for EVERY hospital, so it is not
  # restricted.
  base_dt    = model_base
  ref_dt     = ref_group
  restricted = FALSE
  if (outcome_label %in% sf_outcomes && "pathway_class" %in% names(model_base)) {
    base_dt    = model_base[pathway_class == "a_available_captured"]
    ref_dt     = ref_group[pathway_class == "a_available_captured"]
    restricted = TRUE
  }

  # any_hospice_disposition is a derived flag spanning same-facility AND external
  # hospice, not a level of outcome_category
  outcome_group = if (outcome_label == "any_hospice_disposition") {
    base_dt[any_hospice_disposition == 1L]
  } else {
    base_dt[outcome_category == outcome_label]
  }
  n_outcome = nrow(outcome_group)

  if (n_outcome < 20) {
    message(sprintf("    %s: only %d events — skipping.", outcome_label, n_outcome))
    return(NULL)
  }

  # change 1: control subsampling (no-op unless enabled and the reference
  # group is above the threshold). All events kept; controls sampled within
  # hospital at one global fraction.
  ss = subsample_controls(ref_dt, n_outcome)
  if (ss$sampled)
    message(sprintf("    %s: reference subsampled %s -> %s controls (f = %.5f); intercept on sampled scale, RS outputs offset corrected.",
                    outcome_label,
                    format(nrow(ref_dt), big.mark = ","),
                    format(ss$n_used, big.mark = ","), ss$fraction))

  # combine reference + outcome
  model_dt = rbind(ss$ref, outcome_group)
  model_dt[, y := if (outcome_label == "any_hospice_disposition") {
    as.integer(any_hospice_disposition)
  } else {
    fifelse(outcome_category == outcome_label, 1L, 0L)
  }]

  # hospital counts must reflect the (possibly restricted) analytic sample
  n_hospitals_eff  = model_dt[, uniqueN(hospital_unit)]
  hosp_with_events = model_dt[y == 1L, uniqueN(hospital_unit)]

  # change 4: exact row compression (after subsampling: sample first, then
  # compress). Keys = all.vars of the fitted formula + y, so any random effect
  # structure stays valid automatically. Every fit below runs with
  # weights = row_weight (all 1s when compression is off); the weighted
  # binomial likelihood is identical up to a constant, so only AIC shifts
  # (comparable only within the same compression setting).
  model_dt[, row_weight := 1L]
  compression_ratio = 1
  if (isTRUE(row_compress_enabled)) {
    comp_keys = intersect(unique(c(all.vars(f_glmm), "y")), names(model_dt))
    n_before  = nrow(model_dt)
    model_dt  = model_dt[, .(row_weight = .N), by = comp_keys]
    compression_ratio = round(nrow(model_dt) / n_before, 4)
    message(sprintf("    %s: row compressed %s -> %s rows (ratio %.3f).",
                    outcome_label, format(n_before, big.mark = ","),
                    format(nrow(model_dt), big.mark = ","), compression_ratio))
  }

  if (restricted) {
    message(sprintf("    %s: restricted to class (a) hospitals — %d of %d, %s encounters.",
                    outcome_label, n_hospitals_eff, n_hospitals,
                    format(nrow(model_dt), big.mark = ",")))
  }

  if (n_hospitals_eff < 3 || hosp_with_events < 2) {
    # ADDED 2026.10.16: hospital enters as a FIXED intercept wherever the
    # analytic sample has 2 or more hospitals (add_hospital_fe above)
    hosp_fixed = n_hospitals_eff >= 2L
    hosp_ref   = if (hosp_fixed) add_hospital_fe(model_dt) else NA_character_
    f_fit      = if (hosp_fixed) update(f_fixed, . ~ . + hospital_fe) else f_fixed
    message(sprintf("    %s: insufficient hospital units for a random intercept (%d total, %d with events) -> fixed effects fit%s.",
                    outcome_label, n_hospitals_eff, hosp_with_events,
                    if (hosp_fixed) sprintf(" with hospital fixed intercept (reference %s)", hosp_ref)
                    else " (single hospital, no hospital term)"))

    # fixed-effects logistic regression (no random intercepts)
    tryCatch({
      t_fit = Sys.time()
      fit_fixed = glm(
        f_fit,
        data = model_dt, family = binomial(link = "logit"),
        weights = row_weight            # change 4: all 1s unless compressed
      )
      record_timing(paste0(stage_label, "_fixed"), outcome_label, subgroup_label,
                    nrow(model_dt), n_outcome, t_fit)

      coefs = as.data.table(summary(fit_fixed)$coefficients, keep.rownames = "term")
      setnames(coefs, c("term", "estimate", "se", "z_value", "p_value"))
      coefs[, `:=`(
        or         = round(exp(estimate), 3),
        ci_lower   = round(exp(estimate - 1.96 * se), 3),
        ci_upper   = round(exp(estimate + 1.96 * se), 3),
        outcome    = outcome_label,
        subgroup   = subgroup_label,
        site       = site_lowercase,
        model_type = "fixed_only",
        n_outcome  = n_outcome,
        n_reference = n_ref,
        reference_category = ref_label,
        n_hospitals_eff = n_hospitals_eff,
        aic        = round(AIC(fit_fixed), 1),
        control_sampled  = ss$sampled,
        control_fraction = round(ss$fraction, 6),
        n_reference_used = ss$n_used,
        row_compressed   = isTRUE(row_compress_enabled),
        compression_ratio = compression_ratio,
        hospital_fixed_intercept = hosp_fixed,
        hospital_fe_reference    = hosp_ref
      )]

      # no ICC for fixed-effects model
      icc_row = data.table(
        outcome    = outcome_label,
        subgroup   = subgroup_label,
        site       = site_lowercase,
        model_type = "fixed_only",
        sigma2_hospital = NA_real_,
        icc             = NA_real_,
        n_outcome       = n_outcome,
        n_hospitals_eff     = n_hospitals_eff,
        reference_category  = ref_label
      )

      ranef_dt = NULL

      return(list(coefficients = coefs, icc = icc_row, ranef = ranef_dt))

    }, error = \(e) {
      message(sprintf("    %s: fixed-effects model failed: %s", outcome_label, e$message))
      return(NULL)
    })

  } else {

    # mixed-effects logistic regression
    message(sprintf("    %s: fitting GLMM (n=%s events, %d hospitals)...",
                    outcome_label,
                    format(n_outcome, big.mark = ","),
                    n_hospitals_eff))

    tryCatch({
      t_fit = Sys.time()
      fit = glmer(
        f_glmm,
        data    = model_dt,
        family  = binomial(link = "logit"),
        weights = row_weight,              # change 4: all 1s unless compressed
        nAGQ    = 1,                       # Laplace for speed; 1 is standard
        control = glmerControl(
          optimizer  = "bobyqa",
          optCtrl    = list(maxfun = 50000),
          calc.derivs = FALSE              # skip Hessian for speed
        )
      )
      record_timing(paste0(stage_label, "_glmm"), outcome_label, subgroup_label,
                    nrow(model_dt), n_outcome, t_fit)

      # check convergence
      if (!is.null(fit@optinfo$conv$lme4)) {
        message(sprintf("      Convergence warning: %s",
                        paste(fit@optinfo$conv$lme4$messages, collapse = "; ")))
      }

      # a rank deficient design matrix makes lme4 drop columns silently. Two
      # sites can drop DIFFERENT columns, which misaligns the pooled coefficient
      # table, so name them here rather than leaving it to the console log.
      col_dropped = attr(getME(fit, "X"), "col.dropped")
      if (!is.null(col_dropped) && length(col_dropped) > 0) {
        message(sprintf("      Rank deficient: dropped %s. Usually an empty or collinear factor level.",
                        paste(names(col_dropped), collapse = ", ")))
      }

      # extract fixed effects
      coefs = as.data.table(summary(fit)$coefficients, keep.rownames = "term")
      setnames(coefs, c("term", "estimate", "se", "z_value", "p_value"))
      coefs[, `:=`(
        or         = round(exp(estimate), 3),
        ci_lower   = round(exp(estimate - 1.96 * se), 3),
        ci_upper   = round(exp(estimate + 1.96 * se), 3),
        outcome    = outcome_label,
        subgroup   = subgroup_label,
        site       = site_lowercase,
        model_type = "glmm",
        n_outcome  = n_outcome,
        n_reference = n_ref,
        reference_category = ref_label,
        n_hospitals_eff = n_hospitals_eff,
        aic        = round(AIC(fit), 1),
        control_sampled  = ss$sampled,
        control_fraction = round(ss$fraction, 6),
        n_reference_used = ss$n_used,
        row_compressed   = isTRUE(row_compress_enabled),
        compression_ratio = compression_ratio
      )]

      # ICC: σ²_hospital / (σ²_hospital + π²/3)
      vc = as.data.table(VarCorr(fit))
      sigma2_hosp = vc[grp == "hospital_unit", vcov]
      icc_val = sigma2_hosp / (sigma2_hosp + (pi^2 / 3))

      icc_row = data.table(
        outcome         = outcome_label,
        subgroup        = subgroup_label,
        site            = site_lowercase,
        model_type      = "glmm",
        sigma2_hospital = round(sigma2_hosp, 4),
        icc             = round(icc_val, 4),
        n_outcome       = n_outcome,
        n_hospitals_eff     = n_hospitals_eff,
        reference_category  = ref_label
      )

      message(sprintf("      ICC = %.4f (σ² = %.4f)", icc_val, sigma2_hosp))

      # hospital random intercepts (BLUPs)
      re = ranef(fit)$hospital_unit
      ranef_dt = data.table(
        hospital_unit = rownames(re),
        ranef_intercept = round(re[, 1], 4),
        outcome    = outcome_label,
        subgroup   = subgroup_label,
        reference_category = ref_label,
        site       = site_lowercase
      )

      # protocol 2A.3 profiling suite (supersedes ICC as the headline summary)
      # lp_offset restores population scale when controls were subsampled
      prof = tryCatch(
        hospital_variation_profile(fit, model_dt, outcome_label, subgroup_label,
                                   n_outcome, n_hospitals_eff,
                                   lp_offset = if (ss$sampled) log(ss$fraction) else 0,
                                   std_frame = std_frame, std_source = std_source,
                                   std_version = std_version, ref_label = ref_label),
        error = \(e) {
          message(sprintf("      Hospital profiling failed: %s", e$message))
          NULL
        })

      return(list(coefficients = coefs, icc = icc_row, ranef = ranef_dt,
                  var_summary  = if (!is.null(prof)) prof$summary,
                  hosp_profile = if (!is.null(prof)) prof$hospitals))

    }, error = \(e) {
      message(sprintf("    %s: GLMM failed: %s", outcome_label, e$message))

      # fallback to fixed effects
      # ADDED 2026.10.16: hospital enters as a FIXED intercept (add_hospital_fe
      # above); n_hospitals_eff is 3 or more on this branch, so always on
      hosp_fixed = n_hospitals_eff >= 2L
      hosp_ref   = if (hosp_fixed) add_hospital_fe(model_dt) else NA_character_
      f_fit      = if (hosp_fixed) update(f_fixed, . ~ . + hospital_fe) else f_fixed
      message(sprintf("    %s: falling back to fixed-effects model%s.", outcome_label,
                      if (hosp_fixed) sprintf(" with hospital fixed intercept (reference %s)", hosp_ref)
                      else ""))
      tryCatch({
        t_fit = Sys.time()
        fit_fixed = glm(
          f_fit,
          data = model_dt, family = binomial(link = "logit"),
          weights = row_weight          # change 4: all 1s unless compressed
        )
        record_timing(paste0(stage_label, "_fixed_fallback"), outcome_label,
                      subgroup_label, nrow(model_dt), n_outcome, t_fit)

        coefs = as.data.table(summary(fit_fixed)$coefficients, keep.rownames = "term")
        setnames(coefs, c("term", "estimate", "se", "z_value", "p_value"))
        coefs[, `:=`(
          or         = round(exp(estimate), 3),
          ci_lower   = round(exp(estimate - 1.96 * se), 3),
          ci_upper   = round(exp(estimate + 1.96 * se), 3),
          outcome    = outcome_label,
          subgroup   = subgroup_label,
          site       = site_lowercase,
          model_type = "fixed_fallback",
          n_outcome  = n_outcome,
          n_reference = n_ref,
          reference_category = ref_label,
          n_hospitals_eff = n_hospitals_eff,
          aic        = round(AIC(fit_fixed), 1),
          control_sampled  = ss$sampled,
          control_fraction = round(ss$fraction, 6),
          n_reference_used = ss$n_used,
          row_compressed   = isTRUE(row_compress_enabled),
          compression_ratio = compression_ratio,
          hospital_fixed_intercept = hosp_fixed,
          hospital_fe_reference    = hosp_ref
        )]

        icc_row = data.table(
          outcome = outcome_label, subgroup = subgroup_label,
          site = site_lowercase, model_type = "fixed_fallback",
          sigma2_hospital = NA_real_, icc = NA_real_,
          n_outcome = n_outcome, n_hospitals_eff = n_hospitals_eff,
          reference_category = ref_label
        )

        return(list(coefficients = coefs, icc = icc_row, ranef = NULL))

      }, error = \(e2) {
        message(sprintf("    %s: fixed-effects fallback also failed: %s",
                        outcome_label, e2$message))
        return(NULL)
      })
    })
  }
}


# --- Change 4 payoff diagnostic (always printed) ---
# How much would exact row compression shrink the design? Computed on the
# primary model variables over the END OF LIFE sample (y adds at most a factor
# of ~2 per fit). Decide row_compress_enabled from this line in the site log.
# With the primary fits now a few percent of admissions, compression is
# unlikely to be worth switching on; the incidence fits carry their own cost
# through subsampling.
diag_vars  = intersect(all.vars(formula_glmm), names(model_eol))
n_patterns = nrow(unique(model_eol[, ..diag_vars]))
message(sprintf("\n  Row compression diagnostic (eol): %s unique covariate patterns / %s encounters (shrink factor %.2f). row_compress_enabled = %s.",
                format(n_patterns, big.mark = ","),
                format(nrow(model_eol), big.mark = ","),
                nrow(model_eol) / max(n_patterns, 1L),
                row_compress_enabled))
rm(diag_vars, n_patterns)


# --- fit one model family ---------------------------------------------------
# ADDED 2026.10.14. fit_binary_glmm() reads model_base, ref_group, n_ref and
# n_hospitals from the global environment (so the §3.4.4 class (a) restriction
# and the reference group are applied consistently inside it). The subgroup
# blocks used to swap those globals by hand and restore them afterwards; with
# eight families that is now done here, once. The globals are restored on
# exit whatever happens inside the fits.
fit_family = function(dt, subgroup_label, ref_label, outcomes,
                      f_fixed = formula_fixed, f_glmm = formula_glmm,
                      stage_label = "primary",
                      std_frame = NULL, std_source = "site_local",
                      std_version = "site_local",
                      min_ref = 20L, min_rows = 100L) {
  saved = list(model_base = model_base, ref_group = ref_group,
               n_ref = n_ref, n_hospitals = n_hospitals)
  on.exit({
    model_base  <<- saved$model_base
    ref_group   <<- saved$ref_group
    n_ref       <<- saved$n_ref
    n_hospitals <<- saved$n_hospitals
  }, add = TRUE)

  model_base  <<- dt
  ref_group   <<- dt[outcome_category == ref_label]
  n_ref       <<- nrow(ref_group)
  n_hospitals <<- uniqueN(dt$hospital_unit)

  message(sprintf("\n  === %s | %s: %s encounters, %s reference (%s), %d hospitals ===",
                  stage_label, subgroup_label,
                  format(nrow(dt), big.mark = ","),
                  format(n_ref, big.mark = ","), ref_label, n_hospitals))

  if (n_ref < min_ref || nrow(dt) < min_rows) {
    message("  Too small for modeling -- skipping.")
    return(list())
  }
  res = lapply(outcomes, \(out) {
    fit_binary_glmm(out, subgroup_label, f_fixed, f_glmm, stage_label,
                    ref_label = ref_label, std_frame = std_frame,
                    std_source = std_source, std_version = std_version)
  })
  names(res) = outcomes
  Filter(Negate(is.null), res)
}


# --- PRIMARY: end of life cohort, outcomes vs in hospital death -------------
# subgroup = "eol". same_facility_hospice is fit at class (a) hospitals only
# (inside fit_binary_glmm); external_hospice_discharge and
# any_hospice_disposition at every hospital.

eol_results = fit_family(model_eol, "eol", eol_ref, eol_outcomes,
                         std_frame = std_frame_eol,
                         std_source = standardization_source_eol,
                         std_version = standard_pop_eol_version)


# --- LST and ICU subgroups within the end of life cohort ---------------------

eol_lst_results = fit_family(model_eol[lst_any_flag == 1L], "eol_lst_any",
                             eol_ref, eol_outcomes,
                             std_frame = std_frame_eol,
                             std_source = standardization_source_eol,
                             std_version = standard_pop_eol_version)

# §3.4.3 sensitivity
eol_icu_results = fit_family(model_eol[icu_flag == 1L], "eol_icu_only",
                             eol_ref, eol_outcomes,
                             std_frame = std_frame_eol,
                             std_source = standardization_source_eol,
                             std_version = standard_pop_eol_version)


# --- Opportunity sensitivity: drop early in hospital deaths ------------------
# In hospital deaths within early_death_window_hours of admission are removed
# from the reference group (early_death is 0 for every hospice row).

eol_early_results = fit_family(model_eol[early_death == 0L],
                               sprintf("eol_excl_death_%dh", early_death_window_hours),
                               eol_ref, eol_outcomes,
                               std_frame = std_frame_eol,
                               std_source = standardization_source_eol,
                               std_version = standard_pop_eol_version)


# --- Secondary clinical-course models (end of life cohort) -------------------
# Protocol 2A.3: this model conditions on post-admission clinical course, so it
# answers a DIFFERENT question than the primary model and is reported separately.
# It is never merged into all_coefs / all_icc. Fit on the end of life cohort,
# like the primary.

if (isTRUE(run_clinical_course)) {
  course_results = fit_family(model_eol, "eol", eol_ref, eol_outcomes,
                              formula_fixed_course, formula_glmm_course,
                              stage_label = "clinical_course",
                              std_frame = std_frame_eol,
                              std_source = standardization_source_eol,
                              std_version = standard_pop_eol_version)
} else {
  message("  DISABLED by run_clinical_course = FALSE (recorded in aim2a3_run_config).")
  course_results = list()
}


# --- INCIDENCE fits: full one admission sample, outcomes vs routine discharge -
# subgroup = "full_cohort_incidence". How often admissions end in hospice or
# in hospital death, which the end of life fits cannot say. These are the only
# fits with a routine discharge reference group, so the only ones that still
# subsample controls at large sites (change 1) and the only ones that use the
# all admissions standard population.

incidence_results = fit_family(model_base, "full_cohort_incidence",
                               incidence_ref, incidence_outcomes,
                               std_frame = std_frame,
                               std_source = standardization_source,
                               std_version = standard_pop_version)


# --- Terminal cohort (death linked sites): membership by vital status -------
# subgroup = "terminal_30d". The same outcomes as eol plus died_30d_no_hospice,
# all vs in hospital death. Standardized on the end of life population so its
# profile is comparable with the eol profile; the concordance table in the
# export section compares the two sets of hospital effects.

if (!is.null(model_term)) {
  terminal_results = fit_family(model_term,
                                sprintf("terminal_%dd", terminal_death_window_days),
                                eol_ref, terminal_outcomes,
                                std_frame = std_frame_eol,
                                std_source = standardization_source_eol,
                                std_version = standard_pop_eol_version)
} else {
  terminal_results = list()
}


# --- Hospital characteristics (conditional on reference file) ---

# Bed size, teaching status, and other hospital-level covariates for examining
# whether hospice-to-mortality ratios vary by hospital type (Analysis 2A.2).
# Expected file: reference_files/hospital_characteristics.csv
#   Columns: medicare_id, bed_size, teaching_status, [others]

hosp_char_path = file.path(project_root, "reference_files", "hospital_characteristics.csv")

if (file.exists(hosp_char_path)) {
  message("\n  Loading hospital characteristics...")

  hosp_char = fread(hosp_char_path, colClasses = c(medicare_id = "character"))
  hosp_char[, medicare_id := stringr::str_pad(medicare_id, 6, "left", "0")]

  # merge onto analytic summary (already saved above, reload for enrichment)
  hosp_char_cols = intersect(names(hosp_char), c("medicare_id", "bed_size",
                                                 "teaching_status", "icu_beds", "has_hospice_unit"))
  hosp_char_slim = hosp_char[, ..hosp_char_cols]

  # save enriched hospital file for cross-site compilation
  fwrite(hosp_char_slim,
         shareable_path(
           paste0("aim2a_hospital_characteristics_", site_lowercase, "_", today, ".csv")))
  message(sprintf("  Hospital characteristics loaded for %d hospitals.", nrow(hosp_char_slim)))
  rm(hosp_char, hosp_char_slim)
} else {
  message("\n  Hospital characteristics file not found at: ", hosp_char_path)
  message("  To examine variation by hospital type, create reference_files/hospital_characteristics.csv")
  message("  with columns: medicare_id, bed_size, teaching_status")
  message("  Sources: AHA Annual Survey, CMS Provider of Services, or manual compilation.")
}


# =============================================================================
# Analysis 2.2 (primary): Multinomial mixed model
# =============================================================================

# §3.6.2: multinomial is the PRIMARY model (respects mutual exclusivity of
# disposition categories); binary GLMMs above are the complementary 2.3 analysis.

message("\n--- Analysis 2.2: Multinomial mixed model (primary) ---")

# BUGFIX 2026.08.27: the Option C consistency checks below were guarded by
# exists("all_var_summary"), but all_var_summary was only assembled in the
# export section AFTER this block, so on a fresh session the guard was always
# FALSE and aim2a3_model_consistency was never written at any site (dead code).
# Assemble the binary GLMM profiling results here, before the multinomial fit;
# the export section rebuilds the same objects later, unchanged.
all_results_pre = c(eol_results, eol_lst_results, eol_icu_results, eol_early_results,
                    incidence_results, terminal_results)
all_var_summary = rbindlist(
  Filter(Negate(is.null), lapply(all_results_pre, \(r) r$var_summary)), fill = TRUE)
all_hosp_profile = rbindlist(
  Filter(Negate(is.null), lapply(all_results_pre, \(r) r$hosp_profile)), fill = TRUE)
rm(all_results_pre)

mn_fitted          = FALSE
mn_error           = NA_character_
mn_model_type      = "none"
mn_fallback_reason = NA_character_
mn_hosp_fixed      = FALSE            # ADDED 2026.10.16: fallback hospital fixed intercept
mn_hosp_ref        = NA_character_
mn_ok = requireNamespace("mclogit", quietly = TRUE)

# CHANGED 2026.10.14: the multinomial is fit on the END OF LIFE cohort at
# CLASS (a) hospitals only: three levels (in hospital death, same facility
# hospice, external hospice) against in hospital death. At class (b) and (c)
# hospitals same facility hospice is a structural zero, which a three level
# model cannot carry, so those hospitals are left out here; the network wide
# composition that exists for every hospital is the any_hospice_disposition
# binary vs death in the eol family. mn_dt is built before the switch so that
# a site with NO class (a) hospitals (UMN until its ADT mapping is fixed)
# records a SKIP with a reason rather than a failed fit.
if ("pathway_class" %in% names(model_eol)) {
  mn_dt = model_eol[pathway_class == "a_available_captured"]
} else {
  message("  NOTE: pathway_class not found -- multinomial fit on ALL hospitals; same facility hospice may be a structural zero at some.")
  mn_dt = model_eol
}

if (!isTRUE(run_multinomial)) {
  # Package A: skipped by switch, not failed — the status file says so
  mn_model_type      = "disabled"
  mn_fallback_reason = "disabled by run_multinomial switch (Package A reduction)"
  message("  DISABLED by run_multinomial = FALSE (recorded in aim2a3_run_config and aim2_multinomial_status).")
} else if (mn_ok && nrow(mn_dt) == 0L) {
  mn_model_type      = "skipped"
  mn_fallback_reason = "no class (a) hospitals in the end of life cohort; multinomial not attempted"
  message("  SKIPPED: no class (a) hospitals in the end of life cohort (recorded in aim2_multinomial_status).")
} else if (mn_ok) {
  library(mclogit)

  # change 1: the reference category is sampled the same way as the binary
  # reference groups. With in hospital death as the reference this is a no op
  # at every site (deaths never reach the threshold); the machinery is kept
  # so the composition rescale below stays correct if the threshold changes.
  mn_ref_rows = mn_dt[outcome_category == eol_ref]
  mn_events   = mn_dt[outcome_category != eol_ref]
  ss_mn = subsample_controls(mn_ref_rows, nrow(mn_events))
  if (ss_mn$sampled) {
    message(sprintf("  Multinomial: reference (%s) subsampled %s -> %s (f = %.5f).",
                    eol_ref, format(nrow(mn_ref_rows), big.mark = ","),
                    format(ss_mn$n_used, big.mark = ","), ss_mn$fraction))
    mn_dt = rbind(ss_mn$ref, mn_events)
  }
  rm(mn_ref_rows, mn_events)

  mn_dt[, outcome_mn := relevel(factor(outcome_category), ref = eol_ref)]

  message(sprintf("  Multinomial: %s encounters, %d hospitals, %d outcome levels.",
                  format(nrow(mn_dt), big.mark = ","),
                  uniqueN(mn_dt$hospital_unit), nlevels(mn_dt$outcome_mn)))

  # BUGFIX 2026.08.27 (fifth): a CONSTANT design column is fatal to mblogit —
  # the random effects fit dies with "'a' is computationally singular,
  # rcond(a)=0" and the fixed fit with "no valid set of coefficients has been
  # found: please supply starting values" (both reproduced empirically).
  # glm/glmer silently DROP aliased columns, so the binary fits never noticed;
  # mblogit does not. The known offender is apache2_missing, constant 0
  # wherever the site has no missing APACHE scores (OHSU: 0% missing at both
  # hospitals; the binary fits show "dropped apache2_missing"), which is what
  # killed the multinomial at OHSU and on the same evidence most likely at
  # UMN. Keep only terms that vary in the multinomial analytic data, then
  # verify the remaining design is full rank on a sample.
  mn_terms = strsplit(base_rhs, "\\s*\\+\\s*")[[1]]
  mn_keep = vapply(mn_terms, \(tm) {
    if (!tm %in% names(mn_dt)) return(TRUE)
    uniqueN(mn_dt[[tm]][!is.na(mn_dt[[tm]])]) >= 2L
  }, logical(1))
  if (any(!mn_keep))
    message(sprintf("  Multinomial: dropping constant term(s) from the design: %s",
                    paste(mn_terms[!mn_keep], collapse = ", ")))
  mn_rhs = paste(mn_terms[mn_keep], collapse = " + ")
  mn_formula = as.formula(paste("outcome_mn ~", mn_rhs))

  # rank diagnostic on a sample (a rare level absent from the sample can flag
  # falsely, so this warns rather than stops)
  mn_rank_ok = tryCatch({
    Xs = model.matrix(as.formula(paste("~", mn_rhs)),
                      mn_dt[sample(.N, min(.N, 200000L))])
    qr(Xs)$rank == ncol(Xs)
  }, error = function(e) TRUE)
  if (!isTRUE(mn_rank_ok))
    message("  WARNING: multinomial design may still be rank deficient (aliased factor level?); mblogit can fail on this.")

  # Random structure mirrors the binary GLMMs (system term where a site spans
  # multiple systems; no patient term, since model_base holds one admission
  # per patient). mblogit takes a LIST of formulas.
  mn_random = c(list(~ 1 | hospital_unit),
                if (exists("use_system_re") && isTRUE(use_system_re)) list(~ 1 | health_system))
  if (length(mn_random) == 1L) mn_random = mn_random[[1]]
  message(sprintf("  Multinomial random effects: %s",
                  paste(vapply(if (is.list(mn_random)) mn_random else list(mn_random),
                               \(f) deparse(f[[2]]), character(1)), collapse = " + ")))

  # ADDED 2026.08.27: fixed effects fallback, mirroring the binary GLMM
  # convention. Below the 3 hospital floor the random effect inner solve is
  # singular (observed at OHSU, 2 units: "'a' is computationally singular"),
  # so the random effects fit is not attempted; and any random effects error
  # falls back the same way. The fallback still delivers the protocol 3.6.2
  # coefficients (model_type = multinomial_fixed); the risk standardized
  # composition and the consistency checks are skipped because both need
  # hospital random effects. CHANGED 2026.10.16: as in the binary fallbacks,
  # hospital now enters the fallback as a FIXED intercept wherever the
  # multinomial sample has 2 or more hospitals (add_hospital_fe, section 4),
  # so the fallback ORs are within hospital and the intercept level
  # clustering is absorbed by the hospital dummies.
  mn_hosp_n      = mn_dt[, uniqueN(hospital_unit)]
  mn_hosp_events = mn_dt[outcome_category != eol_ref, uniqueN(hospital_unit)]
  mn_use_re      = mn_hosp_n >= 3 && mn_hosp_events >= 2

  t_fit_mn = Sys.time()
  mn_fit = NULL
  if (mn_use_re) {
    mn_fit = tryCatch(
      mblogit(mn_formula, data = mn_dt, random = mn_random,
              control = mmclogit.control(maxit = 100)),
      error = function(e) {
        mn_error <<- conditionMessage(e)
        message("  Multinomial (random effects) failed: ", e$message)
        NULL
      }
    )
    if (is.null(mn_fit))
      mn_fallback_reason = paste("random effects fit failed:", mn_error)
  } else {
    mn_fallback_reason = sprintf(
      "below hospital floor (%d units, %d with events); random effect not estimable",
      mn_hosp_n, mn_hosp_events)
    message(sprintf("  Multinomial: %s -> fixed effects fallback.", mn_fallback_reason))
  }

  if (!is.null(mn_fit)) {
    mn_fitted = TRUE
    mn_model_type = "multinomial_glmm"
    record_timing("multinomial", "all_dispositions", "eol",
                  nrow(mn_dt), nrow(mn_dt[outcome_category != eol_ref]),
                  t_fit_mn)
    mn_coefs = as.data.table(summary(mn_fit)$coefficients, keep.rownames = "term")
    setnames(mn_coefs, 1:5, c("term", "estimate", "se", "z_value", "p_value"))
    mn_coefs[, `:=`(
      or       = round(exp(estimate), 3),
      ci_lower = round(exp(estimate - 1.96 * se), 3),
      ci_upper = round(exp(estimate + 1.96 * se), 3),
      model_type = "multinomial_mblogit", site = site_lowercase,
      subgroup = "eol", reference_category = eol_ref,
      # change 1: intercept rows are on the sampled scale when subsampled
      # (population intercept = sampled + log(control_fraction))
      control_sampled  = ss_mn$sampled,
      control_fraction = round(ss_mn$fraction, 6),
      n_reference_used = ss_mn$n_used
    )]

    # BUGFIX 2026.08.27: current mclogit (0.9.x) no longer exports a VarCov()
    # function; the random effect covariance lives on the fit object as
    # mn_fit$VarCov (a list by grouping factor). The old call errored with
    # "could not find function", the tryCatch turned that into NULL, and
    # aim2_multinomial_icc plus the consistency checks silently skipped even
    # when the multinomial fitted. Use the component directly, falling back to
    # the function for older mclogit versions.
    mn_vc = mn_fit$VarCov
    if (is.null(mn_vc))
      mn_vc = tryCatch(VarCov(mn_fit), error = function(e) NULL)
    if (!is.null(mn_vc)) {
      mn_icc = rbindlist(lapply(seq_along(diag(mn_vc[[1]])), \(i) {
        s2 = diag(mn_vc[[1]])[i]
        # dimnames carry a "~1" suffix (e.g. "in_hospital_death_acute~1")
        data.table(outcome = sub("~.*$", "", names(diag(mn_vc[[1]]))[i]),
                   sigma2_hospital = round(s2, 4),
                   icc = round(s2 / (s2 + pi^2/3), 4),
                   model_type = "multinomial_mblogit", site = site_lowercase,
                   subgroup = "eol", reference_category = eol_ref)
      }))
      fwrite(mn_icc, shareable_path(
        paste0("aim2_multinomial_icc_", site_lowercase, "_", today, ".csv")))
    }

    fwrite(mn_coefs, shareable_path(
      paste0("aim2_multinomial_coefficients_", site_lowercase, "_", today, ".csv")))
    message("  Saved multinomial model results.")

    # -- Risk-standardized disposition composition (protocol 2A.3, Option C) -----
    # For each hospital, predict the full disposition distribution for a COMMON
    # standardized patient population (the same rows for every hospital), then
    # average. Because the multinomial is one joint model, the predicted category
    # probabilities for a patient sum to 1 by construction, so each hospital's
    # standardized composition also sums to 1 -- the property separate
    # binary models cannot provide, and the reason the composition comes from
    # this model rather than from the GLMMs.
    mn_composition = tryCatch({
      std_n  = 10000L
      set.seed(20260713)
      # end of life standard population (section 2f), never the all
      # admissions one: this composition is conditional on death or hospice
      std_dt = if (!is.null(std_frame_eol)) copy(std_frame_eol)
      else if (nrow(mn_dt) > std_n) mn_dt[sample(.N, std_n)] else copy(mn_dt)
      hosp_ids = mn_dt[, sort(unique(hospital_unit))]
      lvls = levels(mn_dt$outcome_mn)

      comp = rbindlist(lapply(hosp_ids, \(h) {
        nd = copy(std_dt); nd[, hospital_unit := h]
        # mblogit predict returns an N x (#categories) matrix of probabilities
        pr = predict(mn_fit, newdata = nd, type = "response")
        pr = matrix(pr, ncol = length(lvls), dimnames = list(NULL, lvls))
        # change 1: predictions are on the sampled scale when the reference
        # category was subsampled; scale it back up by 1/f and renormalize
        # so the composition is population scale and still sums to 1.
        if (ss_mn$sampled) {
          pr[, eol_ref] = pr[, eol_ref] / ss_mn$fraction
          pr = pr / rowSums(pr)
        }
        data.table(hospital_unit = h, disposition = lvls,
                   rs_probability = round(colMeans(pr), 5))
      }))
      comp[, `:=`(site = site_lowercase, model = "multinomial",
                  subgroup = "eol", reference_category = eol_ref,
                  standardization = standardization_source_eol,
                  standard_pop_version = standard_pop_eol_version,
                  n_standardization = nrow(std_dt),
                  control_sampled = ss_mn$sampled,
                  control_fraction = round(ss_mn$fraction, 6))]
      comp
    }, error = \(e) { message("  Multinomial standardization failed: ", e$message); NULL })

    if (!is.null(mn_composition)) {
      # QC: every hospital's standardized composition must sum to ~1
      sums = mn_composition[, .(tot = sum(rs_probability)), by = hospital_unit]
      if (any(abs(sums$tot - 1) > 1e-3)) {
        message(sprintf("  WARNING: standardized composition does not sum to 1 (range %.4f-%.4f).",
                        min(sums$tot), max(sums$tot)))
      } else {
        message(sprintf("  Standardized composition sums to 1 for all %d hospitals (max dev %.1e).",
                        nrow(sums), max(abs(sums$tot - 1))))
      }
      fwrite(mn_composition, shareable_path(
        paste0("aim2a3_rs_composition_", site_lowercase, "_", today, ".csv")))
      message("  Saved aim2a3_rs_composition (risk-standardized disposition composition).")
    }

    # -- Option C consistency checks (protocol 2A.3 design) ---------------------
    # The multinomial and the binary GLMMs estimate the SAME contrasts (both use
    # in_hospital_death_acute as reference), so they should agree. These run
    # only where the binary GLMMs also fitted, i.e. multi-hospital sites; at a
    # single- or two-hospital site neither model fits and the block is silently
    # skipped.
    # BUGFIX 2026.08.27: guard formerly required exists("all_var_summary"),
    # which was never true here on a fresh run (assembled only in the export
    # section); it is now built above, before this block. The binary
    # comparators are restricted to the primary eol spec, which is what the
    # multinomial is fit on (the old unfiltered merge would have duplicated
    # rows across subgroups). NOTE 2026.10.14: the external_hospice_discharge
    # binary runs at every hospital while the multinomial runs at class (a)
    # hospitals only, so that comparison is like for like only where every
    # hospital is class (a); the same_facility_hospice comparison always is.
    if (nrow(all_var_summary) > 0 && !is.null(mn_vc)) {

      # Check 1: sigma^2 ordering. mblogit is PQL (biased downward); lme4 is
      # Laplace. So for each disposition the multinomial sigma^2 should be <= the
      # binary sigma^2. A multinomial value LARGER than the binary one contradicts
      # the expected PQL bias and points to a convergence problem, not a finding.
      mn_s2 = data.table(outcome = sub("~.*$", "", names(diag(mn_vc[[1]]))),
                         sigma2_multinomial = round(diag(mn_vc[[1]]), 4))
      chk = merge(all_var_summary[subgroup == "eol",
                                  .(outcome, sigma2_binary = sigma2_hospital)],
                  mn_s2, by = "outcome")
      if (nrow(chk)) {
        chk[, sigma2_ordering_ok := sigma2_multinomial <= sigma2_binary * 1.05]
        for (i in seq_len(nrow(chk)))
          message(sprintf("    [consistency] %-28s sigma2 binary=%.4f multinomial=%.4f  %s",
                          chk$outcome[i], chk$sigma2_binary[i], chk$sigma2_multinomial[i],
                          if (chk$sigma2_ordering_ok[i]) "OK"
                          else "** multinomial > binary: check convergence **"))
      }

      # Check 2: hospital random effects should correlate across the two models
      # for each disposition. Best-effort -- mblogit's random-effect storage is
      # version-dependent, so a failure here degrades to a message, not an error.
      chk[, ranef_correlation := NA_real_]
      tryCatch({
        re_mn = mn_fit$random.effects
        re_mn = if (is.list(re_mn)) as.matrix(re_mn[[1]]) else as.matrix(re_mn)
        # hospital_unit is character in this pipeline; levels() on a character
        # vector is NULL, which made this check bail with "unexpected RE shape".
        # Use sorted unique values, matching mblogit's internal factor coding.
        hosp_lv = if (is.factor(mn_dt$hospital_unit)) levels(mn_dt$hospital_unit)
                  else sort(unique(mn_dt$hospital_unit))
        cats    = setdiff(levels(mn_dt$outcome_mn), eol_ref)
        if (nrow(re_mn) == length(hosp_lv) * length(cats)) {
          # rows are stacked hospital-within-category
          bl = data.table(hospital_unit = rep(hosp_lv, times = length(cats)),
                          outcome        = rep(cats, each = length(hosp_lv)),
                          blup_mn        = as.numeric(re_mn))
          bin = all_hosp_profile[subgroup == "eol",
                                 .(hospital_unit, outcome, blup_bin = ranef_intercept)]
          j = merge(bl, bin, by = c("hospital_unit", "outcome"))
          cors = j[, .(r = if (.N >= 3) cor(blup_mn, blup_bin) else NA_real_), by = outcome]
          chk[cors, on = "outcome", ranef_correlation := round(i.r, 3)]
          for (i in seq_len(nrow(cors)))
            if (!is.na(cors$r[i]))
              message(sprintf("    [consistency] %-28s ranef correlation (mn vs binary) = %.3f%s",
                              cors$outcome[i], cors$r[i],
                              if (cors$r[i] < 0.5) "  ** low: possible misspecification **" else ""))
        } else {
          message("    [consistency] ranef correlation skipped (unexpected mblogit RE shape).")
        }
      }, error = \(e) message("    [consistency] ranef correlation unavailable: ", conditionMessage(e)))

      chk[, `:=`(site = site_lowercase, model_spec = "primary")]
      fwrite(chk, shareable_path(
        paste0("aim2a3_model_consistency_", site_lowercase, "_", today, ".csv")))
      message("  Saved aim2a3_model_consistency (multinomial vs binary GLMM agreement).")
    }

    rm(mn_dt, mn_fit, mn_coefs)
    suppressWarnings(rm(mn_vc, mn_icc, mn_composition))
  }

  # fixed effects multinomial fallback (see the fallback design note above)
  if (mn_model_type == "none") {
    # ADDED 2026.10.16: hospital enters the fallback as a FIXED intercept
    # (add_hospital_fe, section 4), mirroring the binary fallbacks. The
    # random effects fit above used mn_formula unchanged.
    mn_hosp_fixed    = mn_hosp_n >= 2L
    mn_hosp_ref      = if (mn_hosp_fixed) add_hospital_fe(mn_dt) else NA_character_
    mn_formula_fixed = if (mn_hosp_fixed) update(mn_formula, . ~ . + hospital_fe) else mn_formula
    message(sprintf("  Fitting fixed effects multinomial fallback (%s)...",
                    if (mn_hosp_fixed) sprintf("hospital fixed intercept, %d units, reference %s",
                                               mn_hosp_n, mn_hosp_ref)
                    else "single hospital, no hospital term"))
    t_fit_mnf = Sys.time()
    mn_fit_fixed = tryCatch(
      mblogit(mn_formula_fixed, data = mn_dt, control = mclogit.control(maxit = 100)),
      error = function(e) {
        mn_error <<- conditionMessage(e)
        message("  Fixed effects multinomial also failed: ", e$message)
        NULL
      }
    )
    if (!is.null(mn_fit_fixed)) {
      mn_fitted     = TRUE
      mn_model_type = "multinomial_fixed"
      record_timing("multinomial_fixed", "all_dispositions", "eol",
                    nrow(mn_dt),
                    nrow(mn_dt[outcome_category != eol_ref]),
                    t_fit_mnf)
      mn_coefs = as.data.table(summary(mn_fit_fixed)$coefficients, keep.rownames = "term")
      setnames(mn_coefs, 1:5, c("term", "estimate", "se", "z_value", "p_value"))
      mn_coefs[, `:=`(
        or       = round(exp(estimate), 3),
        ci_lower = round(exp(estimate - 1.96 * se), 3),
        ci_upper = round(exp(estimate + 1.96 * se), 3),
        model_type = "multinomial_fixed", site = site_lowercase,
        subgroup = "eol", reference_category = eol_ref,
        control_sampled  = ss_mn$sampled,
        control_fraction = round(ss_mn$fraction, 6),
        n_reference_used = ss_mn$n_used,
        hospital_fixed_intercept = mn_hosp_fixed,
        hospital_fe_reference    = mn_hosp_ref
      )]
      fwrite(mn_coefs, shareable_path(
        paste0("aim2_multinomial_coefficients_", site_lowercase, "_", today, ".csv")))
      message("  Saved multinomial coefficients (fixed effects fallback).")
      message("  RS composition and consistency checks skipped: both need hospital random effects.")
      rm(mn_coefs, mn_fit_fixed)
    }
    suppressWarnings(rm(mn_dt, mn_formula_fixed))
  }
} else {
  message("  mclogit not available — multinomial model skipped.")
}

# ADDED 2026.08.27: always record whether the multinomial fitted, so a missing
# coefficient file in a site's shareable bundle is diagnosable at the
# coordinating center rather than only from that site's console log.
fwrite(data.table(site               = site_lowercase,
                  mclogit_available  = mn_ok,
                  multinomial_fitted = isTRUE(mn_fitted),
                  model_type         = mn_model_type,
                  fallback_reason    = fifelse(is.na(mn_fallback_reason), "", mn_fallback_reason),
                  error              = fifelse(is.na(mn_error), "", mn_error),
                  # ADDED 2026.10.16: whether the fallback carried hospital dummies
                  hospital_fixed_intercept = isTRUE(mn_hosp_fixed),
                  hospital_fe_reference    = fifelse(is.na(mn_hosp_ref), "", mn_hosp_ref)),
       shareable_path(paste0("aim2_multinomial_status_", site_lowercase, "_", today, ".csv")))
message(sprintf("  Saved aim2_multinomial_status (fitted = %s, model_type = %s).",
                isTRUE(mn_fitted), mn_model_type))


# =============================================================================
# Analysis 2.5: Temporal variation — ICC per pandemic era
# =============================================================================

# §3.6.2: re-fit binary GLMMs within each pandemic era; report σ²_hospital and
# ICC per era to assess convergence/divergence. Plus Spearman correlation of
# hospital random intercepts across eras.

message("\n--- Analysis 2.5: ICC over pandemic eras ---")

# CHANGED 2026.10.14: refits run inside the END OF LIFE cohort, each outcome
# vs in hospital death, mirroring the primary eol fits (same facility hospice
# at class (a) hospitals only).
era_levels = levels(model_eol$pandemic_period)
era_outcomes = c("same_facility_hospice", "external_hospice_discharge")

# BUGFIX 2026.08.27 (era refits silently produced no output):
#   (a) formula_glmm carries pandemic_period as a fixed effect. Within a single
#       era it has one observed level; model.frame drops unused factor levels,
#       so glmer stopped with "contrasts can be applied only to factors with 2
#       or more levels" on EVERY era fit, the tryCatch swallowed it, and both
#       era files were skipped without a trace. The RHS is now rebuilt per era,
#       keeping only terms that still vary in that era's analytic data (this
#       also guards any other covariate that collapses within an era).
#   (b) the old collection pattern returned list(summary=, ranef=) through
#       nested rbindlist(), which errors on mixed column counts the moment a
#       fit succeeds ("Column 2 of item 1 is length 4 inconsistent..."), and
#       would then have crashed the whole script. Summary and ranef rows are
#       now collected in separate lists.
#   Era fit failures are also now reported loudly instead of returning NULL
#   silently.

era_summary_list = list()
era_ranef_list   = list()

if (!isTRUE(run_era_models)) {
  message("  DISABLED by run_era_models = FALSE (recorded in aim2a3_run_config).")
} else if (length(era_levels) < 2L) {
  # 2026.10.20: one observed era (e.g. data beginning after 2021) -> an era
  # refit would just repeat the primary fit, and nothing exists to compare
  message(sprintf("  SKIPPED: only one pandemic era observed (%s); era refits need 2 or more.",
                  paste(era_levels, collapse = ", ")))
} else for (era in era_levels) {
  era_dt  = model_eol[pandemic_period == era]
  era_ref = era_dt[outcome_category == eol_ref]
  if (nrow(era_dt) < 100 || nrow(era_ref) < 50 || uniqueN(era_dt$hospital_unit) < 3) {
    message(sprintf("    %s: too small — skipping.", era))
    next
  }
  for (out in era_outcomes) {
    # protocol §3.4.4: same facility hospice conditional on class (a) hospitals,
    # as in the primary fits; the reference group is restricted the same way
    if (out %in% sf_outcomes && "pathway_class" %in% names(era_dt)) {
      od = era_dt[pathway_class == "a_available_captured" & outcome_category == out]
      er = era_ref[pathway_class == "a_available_captured"]
    } else {
      od = era_dt[outcome_category == out]
      er = era_ref
    }
    if (nrow(od) < 20 || nrow(er) < 20) {
      message(sprintf("    %s / %s: only %d events, %d reference -- skipping.", era, out, nrow(od), nrow(er)))
      next
    }
    # change 1: era refits use the same control subsampling (sigma^2 only is
    # reported from these fits, so no offset correction is needed)
    ss_era = subsample_controls(er, nrow(od))
    if (ss_era$sampled)
      message(sprintf("    %s / %s: era reference subsampled to %s controls (f = %.5f).",
                      era, out, format(ss_era$n_used, big.mark = ","), ss_era$fraction))
    md = rbind(ss_era$ref, od)
    md[, y := fifelse(outcome_category == out, 1L, 0L)]

    # keep only fixed effect terms that still vary within this era's data
    era_terms = strsplit(base_rhs, "\\s*\\+\\s*")[[1]]
    keep = vapply(era_terms, \(tm) {
      if (!tm %in% names(md)) return(TRUE)
      uniqueN(md[[tm]][!is.na(md[[tm]])]) >= 2L
    }, logical(1))
    if (any(!keep))
      message(sprintf("    %s / %s: dropping constant term(s): %s",
                      era, out, paste(era_terms[!keep], collapse = ", ")))
    f_era = as.formula(paste("y ~", paste(era_terms[keep], collapse = " + "),
                             "+", re_rhs))

    res = tryCatch({
      t_fit = Sys.time()
      fit = glmer(f_era, data = md, family = binomial, nAGQ = 1,
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 30000),
                                         calc.derivs = FALSE))
      record_timing("era_glmm", out, era, nrow(md), nrow(od), t_fit)
      s2 = as.data.table(VarCorr(fit))[grp == "hospital_unit", vcov]
      re = ranef(fit)$hospital_unit
      list(summary = data.table(pandemic_period = era, outcome = out,
                                subgroup = "eol", reference_category = eol_ref,
                                sigma2_hospital = round(s2, 4), icc = round(s2/(s2 + pi^2/3), 4),
                                n_hospitals = uniqueN(md$hospital_unit), site = site_lowercase),
           ranef = data.table(hospital_unit = rownames(re),
                              ranef_intercept = re[,1], era = era, outcome = out))
    }, error = function(e) {
      message(sprintf("    %s / %s: era GLMM failed: %s", era, out, conditionMessage(e)))
      NULL
    })
    if (!is.null(res)) {
      key = paste(era, out, sep = "|")
      era_summary_list[[key]] = res$summary
      era_ranef_list[[key]]   = res$ranef
    }
  }
}

icc_summary = rbindlist(era_summary_list, fill = TRUE)
icc_ranef   = rbindlist(era_ranef_list,   fill = TRUE)

if (nrow(icc_summary) > 0) {

  if (!is.null(icc_summary) && nrow(icc_summary) > 0) {
    message("  ICC by era:")
    for (i in seq_len(nrow(icc_summary)))
      message(sprintf("    %-16s %-30s ICC = %6.4f",
                      icc_summary$pandemic_period[i], icc_summary$outcome[i],
                      icc_summary$icc[i]))
    fwrite(icc_summary, shareable_path(
      paste0("aim2_icc_over_time_", site_lowercase, "_", today, ".csv")))
  }

  # Spearman correlation of hospital random intercepts across eras
  if (!is.null(icc_ranef) && nrow(icc_ranef) > 0) {
    spear_rows = rbindlist(lapply(era_outcomes, \(out) {
      sub = icc_ranef[outcome == out]
      eras_present = unique(sub$era)
      if (length(eras_present) < 2) return(NULL)
      wide = dcast(sub, hospital_unit ~ era, value.var = "ranef_intercept")
      pairs = combn(eras_present, 2, simplify = FALSE)
      rbindlist(lapply(pairs, \(pr) {
        x = wide[[pr[1]]]; y = wide[[pr[2]]]
        ok = !is.na(x) & !is.na(y)
        if (sum(ok) < 4) return(NULL)
        data.table(outcome = out, era_1 = pr[1], era_2 = pr[2],
                   spearman_rho = round(cor(x[ok], y[ok], method = "spearman"), 3),
                   n_hospitals = sum(ok), site = site_lowercase)
      }))
    }))
    if (nrow(spear_rows) > 0) {
      message("  Spearman correlation of random intercepts across eras:")
      for (i in seq_len(nrow(spear_rows)))
        message(sprintf("    %-30s %s vs %s: rho = %.3f",
                        spear_rows$outcome[i], spear_rows$era_1[i],
                        spear_rows$era_2[i], spear_rows$spearman_rho[i]))
      fwrite(spear_rows, shareable_path(
        paste0("aim2_ranef_spearman_", site_lowercase, "_", today, ".csv")))
    }
  }
  rm(icc_summary, icc_ranef)
  suppressWarnings(rm(spear_rows))
}


# 5. Assemble and export results -----------------------------------------------


message("\n--- Assembling results ---")

all_results = c(eol_results, eol_lst_results, eol_icu_results, eol_early_results,
                incidence_results, terminal_results)

# coefficients
all_coefs = rbindlist(
  lapply(all_results, \(r) r$coefficients),
  fill = TRUE
)

# ICC
all_icc = rbindlist(
  lapply(all_results, \(r) r$icc),
  fill = TRUE
)

# random effects
all_ranef = rbindlist(
  Filter(Negate(is.null), lapply(all_results, \(r) r$ranef)),
  fill = TRUE
)


# --- Print ICC summary ---

message("\n========== ICC SUMMARY ==========")
for (sg in unique(all_icc$subgroup)) {
  sub = all_icc[subgroup == sg & !is.na(icc)]
  if (nrow(sub) == 0) {
    message(sprintf("\n  %s: no ICC available (fixed-effects only).", sg))
    next
  }
  message(sprintf("\n  %s:", sg))
  for (i in seq_len(nrow(sub))) {
    message(sprintf("    %-40s ICC = %6.4f  (σ² = %6.4f, n_events = %s)",
                    sub$outcome[i], sub$icc[i], sub$sigma2_hospital[i],
                    format(sub$n_outcome[i], big.mark = ",")))
  }
}
message("==================================\n")


# --- Print key fixed effects ---

message("========== KEY FIXED EFFECTS (end of life cohort, ORs vs in hospital death) ==========")
key_terms = c("age_centered", "icu_flag", "elix_centered", "apache2_centered")
fc_coefs = all_coefs[subgroup == "eol" & term %in% key_terms]
if (nrow(fc_coefs) > 0) {
  for (out in unique(fc_coefs$outcome)) {
    message(sprintf("\n  Outcome: %s", out))
    sub = fc_coefs[outcome == out]
    for (i in seq_len(nrow(sub))) {
      message(sprintf("    %-25s OR = %5.3f [%5.3f, %5.3f]  p = %.4f",
                      sub$term[i], sub$or[i],
                      sub$ci_lower[i], sub$ci_upper[i],
                      sub$p_value[i]))
    }
  }
}
message("==========================================================\n")


# --- Save ---

# label the estimand so primary and clinical-course results can never be confused
if (nrow(all_coefs)) all_coefs[, model_spec := "primary"]
if (nrow(all_icc))   all_icc[,   model_spec := "primary"]

fwrite(
  all_coefs,
  shareable_path(
    paste0("aim2a3_model_coefficients_", site_lowercase, "_", today, ".csv"))
)
message("  Saved aim2a3_model_coefficients.")

# secondary clinical-course model: separate files, explicitly a different estimand
course_coefs = rbindlist(lapply(course_results, \(r) r$coefficients), fill = TRUE)
course_icc   = rbindlist(lapply(course_results, \(r) r$icc),          fill = TRUE)
if (nrow(course_coefs)) {
  # change 3: carry the window and any fallback so the coordinating center can
  # see which sites ran windowed vs whole stay exposures
  course_coefs[, `:=`(model_spec = "clinical_course",
                      course_window_hours = course_window_hours,
                      course_flag_provenance = course_flag_provenance)]
  fwrite(course_coefs, shareable_path(
    paste0("aim2a3_model_coefficients_clinical_course_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim2a3_model_coefficients_clinical_course (secondary estimand).")
}
if (nrow(course_icc)) {
  course_icc[, `:=`(model_spec = "clinical_course",
                    course_window_hours = course_window_hours,
                    course_flag_provenance = course_flag_provenance)]
  fwrite(course_icc, shareable_path(
    paste0("aim2a3_icc_clinical_course_", site_lowercase, "_", today, ".csv")))
}

fwrite(
  all_icc,
  shareable_path(
    paste0("aim2a3_icc_", site_lowercase, "_", today, ".csv"))
)
message("  Saved aim2a3_icc.")

# protocol 2A.3 variation profiling (present only where a GLMM actually fitted)
all_var_summary = rbindlist(
  Filter(Negate(is.null), lapply(all_results, \(r) r$var_summary)), fill = TRUE)
all_hosp_profile = rbindlist(
  Filter(Negate(is.null), lapply(all_results, \(r) r$hosp_profile)), fill = TRUE)

if (nrow(all_var_summary) > 0) {
  all_var_summary[, model_spec := "primary"]
  fwrite(all_var_summary, shareable_path(
    paste0("aim2a3_variation_summary_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim2a3_variation_summary (sigma2, MOR, p10-p90 contrast, reliability).")
} else {
  message("  No GLMM fitted -> no variation summary (expected at single-hospital sites).")
}

if (nrow(all_hosp_profile) > 0) {
  all_hosp_profile[, model_spec := "primary"]
  if (isTRUE(suppress_small_cells_enabled))
    all_hosp_profile = suppress_small_cells(all_hosp_profile,
                                            count_cols = "n_outcome", pct_cols = character(0))
  fwrite(all_hosp_profile, shareable_path(
    paste0("aim2a3_hospital_profile_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim2a3_hospital_profile (risk-standardized probs, reliability, rank CIs).")
}

if (nrow(all_ranef) > 0) {
  fwrite(
    all_ranef,
    shareable_path(
      paste0("aim2a3_hospital_ranef_", site_lowercase, "_", today, ".csv"))
  )
  message("  Saved aim2a3_hospital_ranef.")
}

# ADDED 2026.10.14: concordance of the end of life (disposition defined) and
# terminal (vital status defined) hospital effects at death linked sites.
# Per outcome: the two sigma^2 and the correlation of the hospital random
# intercepts across the two cohorts. High rank correlation means disposition
# defined membership is not distorting the primary results; low correlation
# means it is, and the terminal cohort becomes the one to read. Needs a GLMM
# in both cohorts, so it is absent below the 3 hospital floor.
terminal_label = sprintf("terminal_%dd", terminal_death_window_days)
if (!is.null(model_term) && nrow(all_ranef) > 0 &&
    all(c("eol", terminal_label) %in% all_ranef$subgroup)) {
  conc = merge(all_ranef[subgroup == "eol",
                         .(hospital_unit, outcome, u_eol = ranef_intercept)],
               all_ranef[subgroup == terminal_label,
                         .(hospital_unit, outcome, u_terminal = ranef_intercept)],
               by = c("hospital_unit", "outcome"))
  conc = conc[, .(n_hospitals  = .N,
                  spearman_rho = if (.N >= 4) round(cor(u_eol, u_terminal, method = "spearman"), 3) else NA_real_,
                  pearson_r    = if (.N >= 4) round(cor(u_eol, u_terminal), 3) else NA_real_),
              by = outcome]
  vs2 = all_var_summary[subgroup %in% c("eol", terminal_label),
                        .(outcome, subgroup, sigma2_hospital)]
  if (nrow(vs2)) {
    vs2 = dcast(vs2, outcome ~ subgroup, value.var = "sigma2_hospital")
    setnames(vs2, old = c("eol", terminal_label), new = c("sigma2_eol", "sigma2_terminal"),
             skip_absent = TRUE)
    conc = merge(conc, vs2, by = "outcome", all.x = TRUE)
  }
  conc[, `:=`(site = site_lowercase, reference_category = eol_ref,
              terminal_death_window_days = terminal_death_window_days)]
  fwrite(conc, shareable_path(
    paste0("aim2a3_eol_terminal_concordance_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim2a3_eol_terminal_concordance (end of life vs terminal hospital effects).")
  for (i in seq_len(nrow(conc)))
    message(sprintf("    %-28s Spearman rho = %s over %d hospitals",
                    conc$outcome[i], format(conc$spearman_rho[i]), conc$n_hospitals[i]))
  rm(conc, vs2)
} else if (!is.null(model_term)) {
  message("  Concordance table skipped: need a GLMM (3+ hospitals) in both the eol and terminal cohorts.")
}

# change 2: per model timing (always on) — one row per completed fit
if (length(fit_timings) > 0) {
  timing_dt = rbindlist(fit_timings, fill = TRUE)
  timing_dt[, site := site_lowercase]
  setcolorder(timing_dt, c("site", "stage", "outcome", "subgroup",
                           "rows", "events", "minutes"))
  fwrite(timing_dt, shareable_path(
    paste0("aim2a3_timing_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim2a3_timing.")

  top = timing_dt[order(-minutes)][seq_len(min(5L, .N))]
  message("\n  Five slowest fits:")
  for (i in seq_len(nrow(top)))
    message(sprintf("    %-25s %-30s %-25s %8.1f min (%s rows, %s events)",
                    top$stage[i], top$outcome[i], top$subgroup[i], top$minutes[i],
                    format(top$rows[i], big.mark = ","),
                    format(top$events[i], big.mark = ",")))
}


# 6. save intermediate covariate file ------------------------------------------

# save the covariates so 05_aim1 and others can use them without recomputing
covariate_cols = c(
  "hospitalization_id", "hospital_unit",
  "elixhauser_vw", "elixhauser_count",
  "dx_category", "ccsr_category",
  "apache2_score", "apache2_complete", "apache2_imputed", "apache2_missing",
  "race_eth",
  # change 3: the four early flag columns, so 05 and others can reuse them
  "icu_early_flag", "lst_early_flag", "imv_early_flag", "vaso_early_flag"
)
covariate_cols = intersect(covariate_cols, names(cohort))

save_output(
  cohort[, ..covariate_cols],
  name      = "covariates_aim2a3",
  shareable = FALSE
)
message("  Saved intermediate covariates (covariates_aim2a3).")



# cleanup ----------------------------------------------------------------------

rm(cohort, cohort_all, model_base, model_eol, ref_group,
   analytic_summary, dx_tab,
   all_results, eol_results, eol_lst_results, eol_icu_results, eol_early_results,
   incidence_results, terminal_results,
   all_coefs, all_icc, all_ranef,
   eol_outcomes, incidence_outcomes, terminal_outcomes,
   formula_fixed, formula_glmm, base_rhs)

suppressWarnings(rm(all_var_summary, all_hosp_profile, model_term, mn_dt,
                    std_frame_eol, std_margins, terminal_label,
                    course_results, course_coefs, course_icc,
                    formula_fixed_course, formula_glmm_course, course_rhs))

suppressWarnings(rm(era_summary_list, era_ranef_list, era_levels, era_outcomes,
                    hosp_rates, var_stats, n_hospice_unit_hosp,
                    mn_ok, mn_fitted, mn_error, mn_model_type, mn_fallback_reason,
                    mn_use_re, mn_hosp_n, mn_hosp_events, mn_hosp_fixed, mn_hosp_ref,
                    all_var_summary, all_hosp_profile))

gc()

message("\n===== 04_aim2a3_risk_adjusted_variation.R complete. =====\n")
