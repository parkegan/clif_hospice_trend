# 06_aim2a3_risk_adjusted_variation.R
# Aim 2A.3: Risk-adjusted hospital-level variation in discharge disposition
# Protocol v5.0 Section 3.6.2

# Approach:
#   The protocol specifies a mixed-effects multinomial logistic regression with
#   hospital random intercepts. We implement this as a set of separate binary
#   GLMMs using lme4::glmer().

# Covariates (protocol 3.6.2):
#     - Hospital (random)
#     - Age (continuous, centered)
#     - Sex
#     - Race/ethnicity (collapsed to 5-level)
#     - Admission diagnosis category (ICD-10 based)
#     - Comorbidity burden (Elixhauser van Walraven index)
#     - Admission severity (APACHE II score, first 48h)
#     - ICU utilization (binary)


# Outputs (shareable CSVs in final_output/ for cross-site compilation):
#   - aim2a3_analytic_summary_{site}.csv   — covariate distributions by hospital
#   - aim2a3_model_coefficients_{site}.csv — fixed effects from each binary GLMM
#   - aim2a3_icc_{site}.csv               — ICC per outcome category
#   - aim2a3_hospital_ranef_{site}.csv     — hospital random intercepts (BLUPs)

# Last updated bpe 16.07.2026


# 0. setup ---------------------------------------------------------------------

if (!exists("data_list")) source(here::here("code", "00_setup.R"))

# packages
pkg_needed = c("lme4", "comorbidity", "mclogit")
for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(lme4)
library(ggplot2)
library(patchwork)


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

dx = data_list[["hospital_diagnosis"]] |>
  dplyr::select(dplyr::all_of(unique(c(
    "hospitalization_id", "diagnosis_code", "diagnosis_code_format", dx_flag_cols)))) |>
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

# compute Elixhauser using the comorbidity package
elix = comorbidity::comorbidity(
  x   = dx_icd10[, .(id = hospitalization_id, code = code_clean)],
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

message(sprintf("  Elixhauser: median VW score = %d, mean = %.1f",
                median(cohort$elixhauser_vw),
                mean(cohort$elixhauser_vw)))

rm(elix, elix_scores, elix_dt, flag_cols_elix)


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


## 2d. Race/ethnicity collapsing -----------------------------------------------

# Collapse to 5 levels + unknown to ensure adequate cell sizes.
cohort[, race_eth := fcase(
  tolower(ethnicity_category) %in% c("hispanic", "hispanic or latino"),
  "Hispanic",
  tolower(race_category) %in% c("white", "caucasian"),
  "White non-Hispanic",
  tolower(race_category) %in% c("black", "african american",
                                "black or african american"),
  "Black non-Hispanic",
  tolower(race_category) %in% c("asian"),
  "Asian",
  tolower(race_category) %in% c("other", "multiracial", "multiple",
                                "american indian", "alaska native",
                                "native hawaiian", "pacific islander",
                                "american indian or alaska native",
                                "native hawaiian or other pacific islander"),
  "Other/Multiracial",
  default = "Unknown"
)]


## 2e. Centered and scaled continuous covariates --------------------------------

cohort[, age_centered := age_at_admission - mean(age_at_admission, na.rm = TRUE)]
cohort[, elix_centered := elixhauser_vw - mean(elixhauser_vw, na.rm = TRUE)]

# APACHE II: center among those with scores; impute median for missing.
# Guard the all-missing case: at sites where the APACHE function isn't found
# (apache_fn_path empty) or required inputs are absent, apache2_score is wholly
# NA. median(.) is then NA, imputation fills NA with NA, centering propagates
# NaN, and every row gets dropped by the model_base covariate filter below.
# When that happens, drop APACHE from the model rather than nuking the cohort.
apache_med = median(cohort$apache2_score, na.rm = TRUE)
apache2_available = is.finite(apache_med)

if (apache2_available) {
  cohort[, apache2_imputed  := fifelse(is.na(apache2_score), apache_med, apache2_score)]
  cohort[, apache2_centered := apache2_imputed - mean(apache2_imputed, na.rm = TRUE)]
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
  keep = cohort[, .N, by = admission_source][N >= 50, admission_source]   # collapse rare
  cohort[!admission_source %in% keep, admission_source := "Other"]
  ref_src = names(which.max(table(cohort$admission_source)))
  cohort[, admission_source := relevel(factor(admission_source), ref = ref_src)]
  message(sprintf("  Admission source from '%s' (ref = %s): %d levels.",
                  adm_src_col, ref_src, nlevels(cohort$admission_source)))
} else {
  cohort[, admission_source := factor("Unknown")]   # constant -> excluded
  message("  Admission source: no source-like field in the hospitalization table; omitted.")
}


# 2f. Federation standardization frame (protocol 2A.3, Aim 2A federation) -------
# For hospital estimates to be comparable ACROSS sites, every site must risk-
# standardize onto the SAME patient population. The coordinating center
# distributes a network standard population of RAW covariate profiles; each site
# reapplies ITS OWN centering/level mapping (captured here) and predicts on it.
# Absent a distributed standard, std_frame stays NULL and the profiling falls
# back to site-local standardization -- so single-site OHSU behavior is unchanged.
std_transform = list(
  apache2_available = apache2_available,
  age_mean        = mean(cohort$age_at_admission, na.rm = TRUE),
  elix_mean       = mean(cohort$elixhauser_vw,    na.rm = TRUE),
  apache_med      = apache_med,
  apache_imp_mean = if (apache2_available) mean(cohort$apache2_imputed, na.rm = TRUE) else NA_real_,
  dx_levels       = sort(unique(as.character(cohort$dx_category))),
  race_levels     = levels(factor(cohort$race_eth)),
  sex_levels      = sort(unique(as.character(cohort$sex_category))),
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
  d[, sex_category    := factor(as.character(sex_category),    levels = std_transform$sex_levels)]
  d[, pandemic_period := factor(as.character(pandemic_period), levels = std_transform$era_levels)]
  if (!is.null(std_transform$pcs_levels) && "preexisting_code_status" %in% names(d))
    d[, preexisting_code_status := factor(as.character(preexisting_code_status), levels = std_transform$pcs_levels)]
  if (!is.null(std_transform$src_levels) && "admission_source" %in% names(d))
    d[, admission_source := factor(as.character(admission_source), levels = std_transform$src_levels)]
  keyvars = intersect(c("dx_category", "race_eth", "sex_category", "pandemic_period",
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

# Site covariate margins (shareable, non-identifying) so the coordinating center
# can assemble the pooled network standard population (design decision 1).
.marg_cat = function(col) {
  if (!col %in% names(cohort)) return(NULL)
  # group by the level only (a constant cannot go in `by`, which must be full
  # length); attach the constant covariate name afterward.
  out = cohort[, .(value = round(.N / nrow(cohort), 5)),
               by = .(level = as.character(get(col)))]
  out[, covariate := col]
  out[, .(covariate, level, value)]
}
std_margins = rbindlist(list(
  data.table(covariate = "age_at_admission", level = "mean", value = round(std_transform$age_mean, 3)),
  data.table(covariate = "age_at_admission", level = "sd",   value = round(sd(cohort$age_at_admission, na.rm = TRUE), 3)),
  data.table(covariate = "elixhauser_vw",    level = "mean", value = round(std_transform$elix_mean, 3)),
  data.table(covariate = "elixhauser_vw",    level = "sd",   value = round(sd(cohort$elixhauser_vw, na.rm = TRUE), 3)),
  if (apache2_available) data.table(covariate = "apache2", level = "mean", value = round(std_transform$apache_imp_mean, 3)),
  .marg_cat("sex_category"), .marg_cat("race_eth"), .marg_cat("dx_category"),
  .marg_cat("pandemic_period"), .marg_cat("preexisting_code_status"), .marg_cat("admission_source")
), fill = TRUE)
std_margins[, `:=`(site = site_lowercase, n_cohort = nrow(cohort))]
fwrite(std_margins, shareable_path(
  paste0("aim2a3_standardization_margins_", site_lowercase, "_", today, ".csv")))
message("  Saved aim2a3_standardization_margins (for pooled-standard assembly).")


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
  elix_mean             = round(mean(elixhauser_vw, na.rm = TRUE), 1),
  elix_sd               = round(sd(elixhauser_vw, na.rm = TRUE), 1),
  apache2_mean          = if (apache2_available) round(mean(apache2_imputed, na.rm = TRUE), 1) else NA_real_,
  apache2_sd            = if (apache2_available) round(sd(apache2_imputed,   na.rm = TRUE), 1) else NA_real_,
  pct_apache2_missing   = round(100 * mean(apache2_missing, na.rm = TRUE), 1),
  pct_circulatory       = round(100 * mean(dx_category == "Circulatory", na.rm = TRUE), 1),
  pct_respiratory       = round(100 * mean(dx_category == "Respiratory", na.rm = TRUE), 1),
  pct_infectious        = round(100 * mean(dx_category == "Infectious/parasitic", na.rm = TRUE), 1),
  pct_neoplasm          = round(100 * mean(dx_category == "Neoplasms", na.rm = TRUE), 1),
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
compute_var_stats = function(dt, version_label) {
  if (nrow(dt) < 2) return(NULL)
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

outcomes_to_model = c(
  "in_hospital_death_acute",
  "same_facility_hospice_death",
  "same_facility_hospice_alive",
  "external_hospice_discharge",
  "transfer_acute_facility",
  "other_ama"
)

# protocol §3.4.4: same-facility hospice is only meaningful where the pathway
# exists and is captured, so those outcomes are analysed CONDITIONAL on class (a)
# hospitals. The network-wide outcome, defined for every hospital, is any
# documented hospice disposition (same-facility OR external).
sf_outcomes = c("same_facility_hospice_death", "same_facility_hospice_alive")

if ("any_hospice_disposition" %in% names(cohort)) {
  outcomes_to_model = c(outcomes_to_model, "any_hospice_disposition")
} else {
  message("  NOTE: any_hospice_disposition not found in the cohort -- re-run 01 to build it (§3.4.4).")
}
if (!"pathway_class" %in% names(cohort)) {
  message("  NOTE: pathway_class not found -- same-facility outcomes will NOT be restricted")
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
    !is.na(sex_category) &
    !is.na(race_eth) &
    !is.na(dx_category)
]

# ensure factor levels
model_base[, sex_category := factor(sex_category)]
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
model_base[, pandemic_period := relevel(factor(pandemic_period), ref = "Pre-pandemic")]

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

# reference group
ref_group = model_base[outcome_category == "routine_discharge"]
n_ref = nrow(ref_group)
message(sprintf("  Reference group (routine_discharge): %s encounters.",
                format(n_ref, big.mark = ",")))

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
  c("age_centered", "sex_category", "race_eth", "dx_category",
    "elix_centered",
    if (apache2_available) c("apache2_centered", "apache2_missing"),
    extra_base,
    "pandemic_period"),
  collapse = " + "
)
message(sprintf("  Primary model RHS:         %s", base_rhs))
if (length(extra_base))
  message(sprintf("  v7.0 baseline adjusters included: %s", paste(extra_base, collapse = ", ")))


## random-effect structure (protocol 2A.3) -------------------------------------
# The hospital random intercept is always included. Two further terms are
# protocol-specified but only conditionally estimable:
#
#   (1 | health_system) -- hospitals nested in health systems. A single CLIF site
#     is normally ONE health system, so this variance is not identifiable from
#     site-local data; it is recovered in the pooled coordinating-center model.
#     Included here only if a site genuinely spans several systems.
#
#   (1 | patient_id) -- within-patient dependence from repeat admissions. The
#     protocol asks for this "where estimable"; with hundreds of thousands of
#     patients lme4 cannot fit it in reasonable time, so it is attempted only
#     below patient_re_max_levels. Where it is not estimable, the prespecified
#     one-admission-per-patient sensitivity analysis carries the argument instead
#     (it removes within-patient dependence by construction).

# health_system was populated above, before the reference-group split.
n_systems = model_base[, uniqueN(health_system)]

patient_re_max_levels = 50000L   # beyond this, glmer is impractical
have_pid     = "patient_id" %in% names(model_base)
n_patients   = if (have_pid) model_base[, uniqueN(patient_id)] else NA_integer_
pct_repeat   = if (have_pid) 100 * model_base[, sum(duplicated(patient_id))] / nrow(model_base) else 0

use_system_re  = n_systems >= 3L
use_patient_re = have_pid && pct_repeat >= 2 && n_patients <= patient_re_max_levels

message(sprintf("  Health systems in local data: %d%s", n_systems,
                if (use_system_re) "" else " -> system random effect NOT estimable locally (pooled analysis only)"))
if (have_pid) {
  message(sprintf("  Patients: %s | encounters from repeat patients: %.1f%%",
                  format(n_patients, big.mark = ","), pct_repeat))
  if (use_patient_re) {
    message("  -> including a patient-level random intercept.")
  } else if (pct_repeat < 2) {
    message("  -> repeat admissions negligible; patient random intercept omitted.")
  } else {
    message(sprintf("  -> %s patients exceeds patient_re_max_levels (%s); patient random intercept NOT fitted.",
                    format(n_patients, big.mark = ","), format(patient_re_max_levels, big.mark = ",")))
    message("     Within-patient dependence is addressed by the one-admission-per-patient sensitivity analysis.")
  }
}

re_terms = c("(1 | hospital_unit)",
             if (use_system_re)  "(1 | health_system)",
             if (use_patient_re) "(1 | patient_id)")
re_rhs = paste(re_terms, collapse = " + ")

# random-effect structure WITHOUT the patient term, for the one-admission-per-
# patient sensitivity (where each patient appears once, so the term is degenerate)
re_rhs_nopat = paste(c("(1 | hospital_unit)",
                       if (use_system_re) "(1 | health_system)"), collapse = " + ")

message(sprintf("  Random effects:            %s", re_rhs))

formula_fixed = as.formula(paste("y ~", base_rhs))
formula_glmm  = as.formula(paste("y ~", base_rhs, "+", re_rhs))
formula_glmm_nopat = as.formula(paste("y ~", base_rhs, "+", re_rhs_nopat))

# SECONDARY clinical-course model (protocol 2A.3): adds post-admission clinical
# course back onto the admission-fixed base. The protocol interprets this as a
# DIFFERENT ESTIMAND and requires it be reported separately, so it is fit on the
# full cohort only and written to its own output files -- it never mixes into the
# primary coefficient/ICC tables.
# CAVEAT: icu_flag and lst_any_flag are whole-stay exposures. The protocol asks
# for EARLY organ support (or maximum early severity). Windowing these to the
# first 6/24 hours is a follow-up; until then this model is an approximation.
course_rhs = paste(c(base_rhs, "icu_flag", "lst_any_flag"), collapse = " + ")
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
compute_sigma_ci = TRUE   # profile-likelihood CI for sigma^2. This is the slow
# step; set FALSE at very large sites if runtime bites.

hospital_variation_profile = function(fit, model_dt, outcome_label, subgroup_label,
                                      n_outcome, n_hospitals,
                                      B = 2000, std_n = 10000) {
  
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
    sig_ci = tryCatch({
      pr = suppressWarnings(confint(fit, parm = "theta_", method = "profile",
                                    oldNames = FALSE))
      c(pr[1, 1]^2, pr[1, 2]^2)          # theta is on the SD scale -> square it
    }, error = function(e) {
      message("      sigma^2 profile CI unavailable: ", conditionMessage(e))
      c(NA_real_, NA_real_)
    })
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
  std_dt = if (exists("std_frame") && !is.null(std_frame)) std_frame
  else if (nrow(model_dt) > std_n) model_dt[sample(.N, std_n)] else model_dt
  xb     = predict(fit, newdata = std_dt, re.form = NA, type = "link")
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
    site            = site_lowercase,
    ranef_intercept = round(u, 4),
    ranef_sd        = round(sd_u, 4),
    rs_probability  = round(rs_prob, 5),
    rs_prob_lower   = round(rs_lo, 5),
    rs_prob_upper   = round(rs_hi, 5),
    rs_prob_se      = round((rs_hi - rs_lo) / (2 * 1.96), 5),  # for CC meta-analysis
    standardization = if (exists("standardization_source")) standardization_source else "site_local",
    standard_pop_version = standard_pop_version,
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


# function to fit one binary GLMM
fit_binary_glmm = function(outcome_label, subgroup_label = "full_cohort",
                           f_fixed = formula_fixed, f_glmm = formula_glmm) {
  
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
  
  # combine reference + outcome
  model_dt = rbind(ref_dt, outcome_group)
  model_dt[, y := if (outcome_label == "any_hospice_disposition") {
    as.integer(any_hospice_disposition)
  } else {
    fifelse(outcome_category == outcome_label, 1L, 0L)
  }]
  
  # hospital counts must reflect the (possibly restricted) analytic sample
  n_hospitals_eff  = model_dt[, uniqueN(hospital_unit)]
  hosp_with_events = model_dt[y == 1L, uniqueN(hospital_unit)]
  
  if (restricted) {
    message(sprintf("    %s: restricted to class (a) hospitals — %d of %d, %s encounters.",
                    outcome_label, n_hospitals_eff, n_hospitals,
                    format(nrow(model_dt), big.mark = ",")))
  }
  
  if (n_hospitals_eff < 3 || hosp_with_events < 2) {
    message(sprintf("    %s: insufficient hospital units (%d total, %d with events) — fitting fixed effects only.",
                    outcome_label, n_hospitals_eff, hosp_with_events))
    
    # fixed-effects logistic regression (no random intercepts)
    tryCatch({
      fit_fixed = glm(
        f_fixed,
        data = model_dt, family = binomial(link = "logit")
      )
      
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
        n_hospitals_eff = n_hospitals_eff,
        aic        = round(AIC(fit_fixed), 1)
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
        n_hospitals_eff     = n_hospitals_eff
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
      fit = glmer(
        f_glmm,
        data    = model_dt,
        family  = binomial(link = "logit"),
        nAGQ    = 1,                       # Laplace for speed; 1 is standard
        control = glmerControl(
          optimizer  = "bobyqa",
          optCtrl    = list(maxfun = 50000),
          calc.derivs = FALSE              # skip Hessian for speed
        )
      )
      
      # check convergence
      if (!is.null(fit@optinfo$conv$lme4)) {
        message(sprintf("      Convergence warning: %s",
                        paste(fit@optinfo$conv$lme4$messages, collapse = "; ")))
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
        n_hospitals_eff = n_hospitals_eff,
        aic        = round(AIC(fit), 1)
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
        n_hospitals_eff     = n_hospitals_eff
      )
      
      message(sprintf("      ICC = %.4f (σ² = %.4f)", icc_val, sigma2_hosp))
      
      # hospital random intercepts (BLUPs)
      re = ranef(fit)$hospital_unit
      ranef_dt = data.table(
        hospital_unit = rownames(re),
        ranef_intercept = round(re[, 1], 4),
        outcome    = outcome_label,
        subgroup   = subgroup_label,
        site       = site_lowercase
      )
      
      # protocol 2A.3 profiling suite (supersedes ICC as the headline summary)
      prof = tryCatch(
        hospital_variation_profile(fit, model_dt, outcome_label, subgroup_label,
                                   n_outcome, n_hospitals_eff),
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
      message(sprintf("    %s: falling back to fixed-effects model.", outcome_label))
      tryCatch({
        fit_fixed = glm(
          f_fixed,
          data = model_dt, family = binomial(link = "logit")
        )
        
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
          n_hospitals_eff = n_hospitals_eff,
          aic        = round(AIC(fit_fixed), 1)
        )]
        
        icc_row = data.table(
          outcome = outcome_label, subgroup = subgroup_label,
          site = site_lowercase, model_type = "fixed_fallback",
          sigma2_hospital = NA_real_, icc = NA_real_,
          n_outcome = n_outcome, n_hospitals_eff = n_hospitals_eff
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


# --- Run models for full cohort ---

message("\n  === Full cohort models ===")

full_results = lapply(outcomes_to_model, \(out) {
  fit_binary_glmm(out, "full_cohort")
})
names(full_results) = outcomes_to_model
full_results = Filter(Negate(is.null), full_results)


# --- Run models for LST subgroup ---

message("\n  === LST (any) subgroup models ===")

# temporarily replace model_base with LST subset
model_base_full = copy(model_base)
ref_group_full  = copy(ref_group)
n_ref_full      = n_ref
n_hospitals_full = n_hospitals

model_base = model_base_full[lst_any_flag == 1L]
ref_group  = model_base[outcome_category == "routine_discharge"]
n_ref      = nrow(ref_group)
n_hospitals = uniqueN(model_base$hospital_unit)

if (n_ref >= 20 && nrow(model_base) >= 100) {
  lst_results = lapply(outcomes_to_model, \(out) {
    fit_binary_glmm(out, "lst_any")
  })
  names(lst_results) = outcomes_to_model
  lst_results = Filter(Negate(is.null), lst_results)
} else {
  message("  LST subgroup too small for modeling — skipping.")
  lst_results = list()
}

# restore
model_base  = model_base_full
ref_group   = ref_group_full
n_ref       = n_ref_full
n_hospitals = n_hospitals_full
rm(model_base_full, ref_group_full)


# --- Run models for ICU-only subgroup (§3.4.3 sensitivity) ---

message("\n  === ICU-only subgroup models ===")

model_base_full = copy(model_base)
ref_group_full  = copy(ref_group)
n_ref_full      = n_ref
n_hospitals_full = n_hospitals

model_base = model_base_full[icu_flag == 1L]
ref_group  = model_base[outcome_category == "routine_discharge"]
n_ref      = nrow(ref_group)
n_hospitals = uniqueN(model_base$hospital_unit)

if (n_ref >= 20 && nrow(model_base) >= 100) {
  icu_results = lapply(outcomes_to_model, \(out) {
    fit_binary_glmm(out, "icu_only")
  })
  names(icu_results) = outcomes_to_model
  icu_results = Filter(Negate(is.null), icu_results)
} else {
  message("  ICU subgroup too small for modeling — skipping.")
  icu_results = list()
}

# restore
model_base  = model_base_full
ref_group   = ref_group_full
n_ref       = n_ref_full
n_hospitals = n_hospitals_full
rm(model_base_full, ref_group_full)


# --- Prespecified sensitivity: one admission per patient (protocol 2A.3) ------
# Removes within-patient dependence by construction (each patient contributes a
# single, first, admission). This is the protocol's stated recourse when the
# patient-level random intercept is not estimable. Fitted WITHOUT the patient
# random effect, which would be degenerate here. Results carry
# subgroup = "one_admission_per_patient" in the standard coefficient tables.

if ("patient_id" %in% names(model_base)) {
  
  model_base_full  = copy(model_base)
  ref_group_full   = ref_group
  n_ref_full       = n_ref
  n_hospitals_full = n_hospitals
  
  if ("admission_dttm" %in% names(model_base)) {
    setorder(model_base, patient_id, admission_dttm)
  } else {
    setorder(model_base, patient_id)
  }
  model_base  = model_base[!duplicated(patient_id)]
  ref_group   = model_base[outcome_category == "routine_discharge"]
  n_ref       = nrow(ref_group)
  n_hospitals = uniqueN(model_base$hospital_unit)
  
  message(sprintf("\n  === One-admission-per-patient sensitivity (%s of %s encounters) ===",
                  format(nrow(model_base), big.mark = ","),
                  format(nrow(model_base_full), big.mark = ",")))
  
  if (n_ref >= 20 && nrow(model_base) >= 100) {
    oneadm_results = lapply(outcomes_to_model, \(out) {
      fit_binary_glmm(out, "one_admission_per_patient",
                      formula_fixed, formula_glmm_nopat)
    })
    names(oneadm_results) = outcomes_to_model
    oneadm_results = Filter(Negate(is.null), oneadm_results)
  } else {
    message("  Too small for modeling -- skipping.")
    oneadm_results = list()
  }
  
  # restore the full analytic sample
  model_base  = model_base_full
  ref_group   = ref_group_full
  n_ref       = n_ref_full
  n_hospitals = n_hospitals_full
  rm(model_base_full, ref_group_full)
  
} else {
  oneadm_results = list()
}


# --- Secondary clinical-course models (full cohort) ---
# Protocol 2A.3: this model conditions on post-admission clinical course, so it
# answers a DIFFERENT question than the primary model and is reported separately.
# It is never merged into all_coefs / all_icc.

message("\n  === Secondary clinical-course models (full cohort) ===")

course_results = lapply(outcomes_to_model, \(out) {
  fit_binary_glmm(out, "full_cohort", formula_fixed_course, formula_glmm_course)
})
names(course_results) = outcomes_to_model
course_results = Filter(Negate(is.null), course_results)


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

mn_ok = requireNamespace("mclogit", quietly = TRUE)

if (mn_ok) {
  library(mclogit)
  
  mn_dt = model_base[outcome_category %in% c("routine_discharge", outcomes_to_model)]
  mn_dt[, outcome_mn := relevel(factor(outcome_category), ref = "routine_discharge")]
  
  message(sprintf("  Multinomial: %s encounters, %d hospitals, %d outcome levels.",
                  format(nrow(mn_dt), big.mark = ","),
                  uniqueN(mn_dt$hospital_unit), nlevels(mn_dt$outcome_mn)))
  
  mn_formula = as.formula(paste("outcome_mn ~", base_rhs))
  
  # Random structure mirrors the binary GLMMs (system term where a site spans
  # multiple systems; patient term is NOT added -- mblogit cannot carry hundreds
  # of thousands of patient levels, and within-patient dependence is handled by
  # the one-admission-per-patient sensitivity). mblogit takes a LIST of formulas.
  mn_random = c(list(~ 1 | hospital_unit),
                if (exists("use_system_re") && isTRUE(use_system_re)) list(~ 1 | health_system))
  if (length(mn_random) == 1L) mn_random = mn_random[[1]]
  message(sprintf("  Multinomial random effects: %s",
                  paste(vapply(if (is.list(mn_random)) mn_random else list(mn_random),
                               \(f) deparse(f[[2]]), character(1)), collapse = " + ")))
  
  mn_fit = tryCatch(
    mblogit(mn_formula, data = mn_dt, random = mn_random,
            control = mmclogit.control(maxit = 100)),
    error = function(e) { message("  Multinomial failed: ", e$message); NULL }
  )
  
  if (!is.null(mn_fit)) {
    mn_coefs = as.data.table(summary(mn_fit)$coefficients, keep.rownames = "term")
    setnames(mn_coefs, 1:5, c("term", "estimate", "se", "z_value", "p_value"))
    mn_coefs[, `:=`(
      or       = round(exp(estimate), 3),
      ci_lower = round(exp(estimate - 1.96 * se), 3),
      ci_upper = round(exp(estimate + 1.96 * se), 3),
      model_type = "multinomial_mblogit", site = site_lowercase
    )]
    
    mn_vc = tryCatch(VarCov(mn_fit), error = function(e) NULL)
    if (!is.null(mn_vc)) {
      mn_icc = rbindlist(lapply(seq_along(diag(mn_vc[[1]])), \(i) {
        s2 = diag(mn_vc[[1]])[i]
        data.table(outcome = names(diag(mn_vc[[1]]))[i],
                   sigma2_hospital = round(s2, 4),
                   icc = round(s2 / (s2 + pi^2/3), 4),
                   model_type = "multinomial_mblogit", site = site_lowercase)
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
    # standardized composition also sums to 1 -- the property seven separate
    # binary models cannot provide, and the reason the composition comes from
    # this model rather than from the GLMMs.
    mn_composition = tryCatch({
      std_n  = 10000L
      set.seed(20260713)
      std_dt = if (exists("std_frame") && !is.null(std_frame)) copy(std_frame)
      else if (nrow(mn_dt) > std_n) mn_dt[sample(.N, std_n)] else copy(mn_dt)
      hosp_ids = mn_dt[, sort(unique(hospital_unit))]
      lvls = levels(mn_dt$outcome_mn)
      
      comp = rbindlist(lapply(hosp_ids, \(h) {
        nd = copy(std_dt); nd[, hospital_unit := h]
        # mblogit predict returns an N x (#categories) matrix of probabilities
        pr = predict(mn_fit, newdata = nd, type = "response")
        pr = matrix(pr, ncol = length(lvls), dimnames = list(NULL, lvls))
        data.table(hospital_unit = h, disposition = lvls,
                   rs_probability = round(colMeans(pr), 5))
      }))
      comp[, `:=`(site = site_lowercase, model = "multinomial",
                  standardization = if (exists("standardization_source")) standardization_source else "site_local",
                  standard_pop_version = standard_pop_version,
                  n_standardization = nrow(std_dt))]
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
    # routine_discharge as reference), so they should agree. These run only where
    # the binary GLMMs also fitted, i.e. multi-hospital sites; at a single- or
    # two-hospital site neither model fits and the block is silently skipped.
    if (exists("all_var_summary") && nrow(all_var_summary) > 0 && !is.null(mn_vc)) {
      
      # Check 1: sigma^2 ordering. mblogit is PQL (biased downward); lme4 is
      # Laplace. So for each disposition the multinomial sigma^2 should be <= the
      # binary sigma^2. A multinomial value LARGER than the binary one contradicts
      # the expected PQL bias and points to a convergence problem, not a finding.
      mn_s2 = data.table(outcome = sub("~.*$", "", names(diag(mn_vc[[1]]))),
                         sigma2_multinomial = round(diag(mn_vc[[1]]), 4))
      chk = merge(all_var_summary[, .(outcome, sigma2_binary = sigma2_hospital)],
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
        hosp_lv = levels(mn_dt$hospital_unit)
        cats    = setdiff(levels(mn_dt$outcome_mn), "routine_discharge")
        if (nrow(re_mn) == length(hosp_lv) * length(cats)) {
          # rows are stacked hospital-within-category
          bl = data.table(hospital_unit = rep(hosp_lv, times = length(cats)),
                          outcome        = rep(cats, each = length(hosp_lv)),
                          blup_mn        = as.numeric(re_mn))
          bin = all_hosp_profile[, .(hospital_unit, outcome, blup_bin = ranef_intercept)]
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
} else {
  message("  mclogit not available — multinomial model skipped.")
}


# =============================================================================
# Analysis 2.5: Temporal variation — ICC per pandemic era
# =============================================================================

# §3.6.2: re-fit binary GLMMs within each pandemic era; report σ²_hospital and
# ICC per era to assess convergence/divergence. Plus Spearman correlation of
# hospital random intercepts across eras.

message("\n--- Analysis 2.5: ICC over pandemic eras ---")

era_levels = levels(model_base$pandemic_period)
era_outcomes = c("in_hospital_death_acute", "same_facility_hospice_death")

icc_over_time = rbindlist(lapply(era_levels, \(era) {
  era_dt = model_base[pandemic_period == era]
  era_ref = era_dt[outcome_category == "routine_discharge"]
  if (nrow(era_dt) < 100 || nrow(era_ref) < 50 || uniqueN(era_dt$hospital_unit) < 3) {
    message(sprintf("    %s: too small — skipping.", era))
    return(NULL)
  }
  rbindlist(lapply(era_outcomes, \(out) {
    od = era_dt[outcome_category == out]
    if (nrow(od) < 20) return(NULL)
    md = rbind(era_ref, od)
    md[, y := fifelse(outcome_category == out, 1L, 0L)]
    tryCatch({
      fit = glmer(formula_glmm, data = md, family = binomial, nAGQ = 1,
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 30000),
                                         calc.derivs = FALSE))
      s2 = as.data.table(VarCorr(fit))[grp == "hospital_unit", vcov]
      re = ranef(fit)$hospital_unit
      list(summary = data.table(pandemic_period = era, outcome = out,
                                sigma2_hospital = round(s2, 4), icc = round(s2/(s2 + pi^2/3), 4),
                                n_hospitals = uniqueN(md$hospital_unit), site = site_lowercase),
           ranef = data.table(hospital_unit = rownames(re),
                              ranef_intercept = re[,1], era = era, outcome = out))
    }, error = function(e) NULL)
  }), fill = TRUE)
}), fill = TRUE)

# icc_over_time currently holds nested lists; separate summary and ranef
if (length(icc_over_time) > 0 && nrow(icc_over_time) > 0) {
  icc_summary = rbindlist(lapply(icc_over_time$summary, identity), fill = TRUE)
  icc_ranef   = rbindlist(lapply(icc_over_time$ranef, identity), fill = TRUE)
  
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

all_results = c(full_results, lst_results, icu_results, oneadm_results)

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

message("========== KEY FIXED EFFECTS (full cohort, ORs) ==========")
key_terms = c("age_centered", "icu_flag", "elix_centered", "apache2_centered")
fc_coefs = all_coefs[subgroup == "full_cohort" & term %in% key_terms]
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
  course_coefs[, model_spec := "clinical_course"]
  fwrite(course_coefs, shareable_path(
    paste0("aim2a3_model_coefficients_clinical_course_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim2a3_model_coefficients_clinical_course (secondary estimand).")
}
if (nrow(course_icc)) {
  course_icc[, model_spec := "clinical_course"]
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


# 6. save intermediate covariate file ------------------------------------------

# save the covariates so 05_aim1 and others can use them without recomputing
covariate_cols = c(
  "hospitalization_id", "hospital_unit",
  "elixhauser_vw", "elixhauser_count",
  "dx_category", "ccsr_category",
  "apache2_score", "apache2_complete", "apache2_imputed", "apache2_missing",
  "race_eth"
)
covariate_cols = intersect(covariate_cols, names(cohort))

save_output(
  cohort[, ..covariate_cols],
  name      = "covariates_aim2a3",
  shareable = FALSE
)
message("  Saved intermediate covariates (covariates_aim2a3).")



# cleanup ----------------------------------------------------------------------

rm(cohort, cohort_all, model_base, ref_group,
   analytic_summary, dx_tab,
   all_results, full_results, lst_results, icu_results, oneadm_results,
   all_coefs, all_icc, all_ranef,
   outcomes_to_model, formula_fixed, formula_glmm, base_rhs)

suppressWarnings(rm(all_var_summary, all_hosp_profile,
                    course_results, course_coefs, course_icc,
                    formula_fixed_course, formula_glmm_course, course_rhs))

suppressWarnings(rm(icc_over_time, era_levels, era_outcomes,
                    hosp_rates, var_stats, n_hospice_unit_hosp, mn_ok))

gc()

message("\n===== 04_aim2a3_risk_adjusted_variation.R complete. =====\n")