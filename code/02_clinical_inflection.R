
# 02_clinical_inflection.R

# Clinical inflection point detection and post-inflection trajectory
# Protocol section 3.5.4
#
# Framing: Rather than "clinical inflection" (which implies clinical intent
# we cannot observe), we identify clinical inflection points — observable
# changes in care patterns that precede hospice transition or death.
#
# Stage 1: INFLECTION POINT — the earlier of:
#   (a) First full → non-full code status transition, OR
#   (b) First comfort-oriented LST cessation (vaso withdrawal, MV withdrawal)
#
# Stage 2: POST-INFLECTION EVENTS — what happens after:
#   - Additional code status changes
#   - Opioid initiation/escalation
#   - Lab cessation
#   - Time to hospice transition and death
#
# Components (flagged independently):
#   1. Code status change to comfort-focused within 72h of discharge/death
#      + first_comfort_dttm (Aim 3 timing)
#   2. Opioid infusion initiation/escalation in final 48h
#      + windowed opioid counts: opioid_{72_48h,48_24h,24_0h}
#      + opioid_trajectory classification
#   3. Diagnostic testing cessation: no routine labs in final 48h
#      + windowed lab counts: labs_{72_48h,48_24h,24_0h}
#   4. Withdrawal of life-sustaining treatment:
#      a. Vasopressor discontinuation + last_vaso_dttm
#      b. MV discontinued (not extubated to room air) + last_imv_dttm
#      (CRRT: excluded — limited availability across sites and high overlap
#       with IMV/vasopressor subcohort reduces marginal contribution)
#      + Derived: inflection_point_dttm, first_lst_withdrawal_dttm,
#        hours to hospice/death, code_status_sequence
#
# Exploratory flags:
#   5. Broad medication discontinuation
#   6. Absence of routine labs
#   7. Diminished vital signs frequency in final 24h vs prior 24h
#
# Output: one row per hospitalization_id with all component flags,
#         composites, inflection point timing, and Aim 3 variables

# Last updated bpe 22.07.2026


# 0. setup ---------------------------------------------------------------------

source(file.path(here::here("code"), "00_setup.R"))

# load cohort from 01_cohort.R output
cohort_all = fread(latest_output("cohort_all_hospitalizations"),
                   colClasses = list(character = c("patient_id", "hospitalization_id")))
cohort_all[, hospitalization_id := as.character(hospitalization_id)]
cohort_all[, patient_id := as.character(patient_id)]

# we need a reference endpoint for each encounter:
# use death_dttm if died, otherwise discharge_dttm
cohort_all[, endpoint_dttm := fifelse(
  !is.na(death_dttm) & death_dttm <= discharge_dttm + days(30),
  death_dttm,
  discharge_dttm
)]

# limit to primary cohort
cohort = cohort_all[in_cohort == 1L]

# ensure datetimes are POSIXct
dt_cols = c("admission_dttm", "discharge_dttm", "death_dttm", "endpoint_dttm",
            "hospice_transition_dttm")
for (col in dt_cols) {
  if (col %in% names(cohort) && !inherits(cohort[[col]], "POSIXct")) {
    set(cohort, j = col, value = as.POSIXct(cohort[[col]], tz = "UTC"))
  }
}

hids = unique(cohort$hospitalization_id)
n_cohort = length(hids)

message(sprintf("Building clinical inflection for %s encounters.", format(n_cohort, big.mark = ",")))


# reference: endpoint minus window boundaries
ref = cohort[, .(
  hospitalization_id, patient_id,
  admission_dttm, discharge_dttm, endpoint_dttm,
  hospice_transition_dttm,
  ep_minus_48h = endpoint_dttm - hours(48),
  ep_minus_72h = endpoint_dttm - hours(72),
  ep_minus_24h = endpoint_dttm - hours(24)
)]



# 1. code status change - 3-category mapping (§3.5.4) --------------------------

message("  Component 1: Code status change...")

code_status = data_list[["code_status"]] |>
  dplyr::select(patient_id, start_dttm, code_status_category) |>
  dplyr::collect() |>
  as.data.table()
code_status[, patient_id := as.character(patient_id)]

code_status[, cat_lower := tolower(trimws(code_status_category))]

# Protocol v6.0 §3.5.4 three-category mapping:
#   Full code:            Full, Presume Full
#   Limited resuscitation: DNR, DNAR, UDNR, DNR/DNI, DNAR/DNI
#   Comfort-focused:       AND (Allow Natural Death) = the network comfort-focused signal
#                          (a site may also record an explicit CMO / comfort-care value)
#
# NOTE: The CLIF schema has no dedicated CMO value. Network decision (2026-06):
# rather than asking each site to configure a comfort mapping, the "AND" (Allow
# Natural Death) category is treated as the comfort-focused anchor, fixed
# network-wide so sites have nothing to configure. At sites where "AND"
# is used as a DNR-equivalent, this mapping must be revisited (protocol §3.8,
# code-status mapping heterogeneity). Where comfort is not captured at all, the
# framework still falls back to the Full -> Limited transition (handled in §4e).

full_cats     = c("full", "presume full")
# "special/partial" covers partial-code / limited orders (at OHSU the "LIMITED"
# name maps to this category); note it also sweeps in a few device (VAD) rows
# the category cannot distinguish from genuine limited orders.
limited_cats  = c("dnr", "dnar", "udnr", "dnr/dni", "dnar/dni", "dni_only", "special/partial")
# "AND" (Allow Natural Death) is this network's comfort-focused signal; the CMO
# synonyms are a safety net for any site that records an explicit comfort value.
# Fixed network-wide -- there is nothing for sites to configure.
comfort_focused_cats = c("and", "cmo", "comfort care", "comfort measures only", "comfort")

code_status[, cs_category := fcase(
  cat_lower %in% full_cats,            "Full",
  cat_lower %in% comfort_focused_cats, "Comfort-focused",
  cat_lower %in% limited_cats,         "Limited resuscitation",
  default = "Other"
)]

# any non-Full status counts as a reduction in resuscitation intensity
nonfull_cats = c(limited_cats, comfort_focused_cats)

# merge for endpoint windows
cs = merge(
  code_status,
  ref[, .(hospitalization_id,
          patient_id,
          admission_dttm,
          discharge_dttm,
          ep_minus_72h,
          endpoint_dttm)],
  by = "patient_id",
  allow.cartesian = TRUE
)

# non-Full code status orders within 72h of endpoint
cs_72h = cs[
  start_dttm >= admission_dttm &
    start_dttm <= endpoint_dttm &
    start_dttm >= ep_minus_72h &
    cat_lower %in% nonfull_cats
]

# all code status during hospitalization (for prior-full and timing)
all_cs = merge(
  code_status,
  ref[, .(hospitalization_id, patient_id, admission_dttm, endpoint_dttm)],
  by = "patient_id",
  allow.cartesian = TRUE
)
all_cs = all_cs[start_dttm >= admission_dttm & start_dttm <= endpoint_dttm]

prior_full = all_cs[
  cat_lower %in% full_cats,
  .(had_prior_full = 1L),
  by = hospitalization_id
]

cs_flag = cs_72h[, .(cs_nonfull_72h = 1L), by = hospitalization_id]
cs_flag = merge(cs_flag, prior_full, by = "hospitalization_id", all.x = TRUE)
cs_flag[is.na(had_prior_full), had_prior_full := 0L]

# strict: Full → non-Full transition; lenient: any non-Full order in final 72h
cs_flag[, flag_code_status_strict := fifelse(cs_nonfull_72h == 1L & had_prior_full == 1L, 1L, 0L)]
cs_flag[, flag_code_status_lenient := cs_nonfull_72h]

cs_out = cs_flag[, .(hospitalization_id, flag_code_status_strict, flag_code_status_lenient)]

# Aim 3 timing — capture BOTH:
#   first_comfort_dttm: first comfort-focused order (if captured at this site)
#   first_limited_dttm: first Limited Resuscitation order (Full→Limited fallback anchor)
comfort_timing = all_cs[
  cs_category == "Comfort-focused",
  .(first_comfort_dttm = min(start_dttm)),
  by = hospitalization_id
]
limited_timing = all_cs[
  cs_category == "Limited resuscitation",
  .(first_limited_dttm = min(start_dttm)),
  by = hospitalization_id
]
# first non-Full of any kind
nonfull_timing = all_cs[
  cat_lower %in% nonfull_cats,
  .(first_nonfull_dttm = min(start_dttm)),
  by = hospitalization_id
]

# QI diagnostic: comfort-order capture FUNNEL. When first_comfort_dttm is NA for
# nearly everyone, code_status_sequence (Aim 3) collapses to a single "No comfort
# order" level. The funnel separates the possible causes, which have different
# remedies:
#   (i)   no comfort orders exist at all              -> DATA gap (documented elsewhere)
#   (ii)  they exist but belong to non-cohort patients -> cohort restriction, expected
#   (iii) they exist for cohort patients but sit outside [admission, endpoint]
#         -> WINDOWING gap; standing/pre-admission orders (POLST, outpatient AND)
#            are the usual cause, and are arguably baseline status rather than an
#            in-hospital transition
#   (iv)  they land in-window but not in the LST + hospice subset Aim 3 reports on
#         -> the counts are simply too small to intersect
cs_comfort           = code_status[cs_category == "Comfort-focused"]
n_comfort_classified = nrow(cs_comfort)
n_comfort_cohort     = cs_comfort[patient_id %in% ref$patient_id, .N]
n_enc_comfort_window = nrow(comfort_timing)

message("  Comfort-order capture funnel:")
message(sprintf("    %s order(s) classified Comfort-focused (all patients)",
                format(n_comfort_classified, big.mark = ",")))
message(sprintf("    %s of those belong to cohort patients",
                format(n_comfort_cohort, big.mark = ",")))
message(sprintf("    %s cohort encounter(s) with a comfort order inside [admission, endpoint]",
                format(n_enc_comfort_window, big.mark = ",")))

# where do the out-of-window orders sit relative to the patient's encounters?
if (n_comfort_cohort > 0L) {
  cf = merge(cs_comfort[, .(patient_id, start_dttm)],
             ref[, .(patient_id, admission_dttm, endpoint_dttm)],
             by = "patient_id", allow.cartesian = TRUE)
  cf[, in_win := start_dttm >= admission_dttm & start_dttm <= endpoint_dttm]
  ord = cf[, .(any_in     = any(in_win),
               all_before = all(start_dttm < admission_dttm),
               all_after  = all(start_dttm > endpoint_dttm),
               days_before_first_adm = as.numeric(
                 difftime(min(admission_dttm), start_dttm[1], units = "days"))),
           by = .(patient_id, start_dttm)]
  n_before = ord[!any_in & all_before, .N]
  n_after  = ord[!any_in & all_after,  .N]
  med_lead = ord[!any_in & all_before, median(days_before_first_adm, na.rm = TRUE)]
  message(sprintf("    out-of-window: %s before every admission (median %s days earlier), %s after every endpoint",
                  format(n_before, big.mark = ","),
                  if (is.finite(med_lead)) round(med_lead, 1) else NA,
                  format(n_after, big.mark = ",")))
  suppressWarnings(rm(cf, ord))
}

# do the in-window encounters even reach the Aim 3 LST + hospice denominator?
if (n_enc_comfort_window > 0L && "lst_any_flag" %in% names(cohort)) {
  cw_ids        = comfort_timing$hospitalization_id
  n_cw_lst      = cohort[hospitalization_id %in% cw_ids & lst_any_flag == 1L, .N]
  n_cw_lst_hosp = cohort[hospitalization_id %in% cw_ids & lst_any_flag == 1L &
                           !is.na(hospice_transition_dttm), .N]
  message(sprintf("    of those, %s are LST encounters and %s are LST + hospice (the Aim 3 denominator)",
                  format(n_cw_lst, big.mark = ","), format(n_cw_lst_hosp, big.mark = ",")))
  if (n_cw_lst_hosp == 0L)
    message("    QI: no in-window comfort order falls in the LST + hospice subset, so ",
            "code_status_sequence collapses to a single \"No comfort order\" level.")
}

if (n_comfort_classified == 0L) {
  message("    QI: no code_status values map to Comfort-focused. Comfort care may be ",
          "documented outside code_status at this site; Aim 3 comfort timing will be empty.")
} else if (n_enc_comfort_window == 0L) {
  message("    QI: comfort orders EXIST but none fall within [admission, endpoint] -> ",
          "windowing gap (likely pre-admission standing orders). Aim 3 comfort timing will be empty.")
} else if (n_enc_comfort_window < 0.02 * n_comfort_classified) {
  message("    QI: only a small fraction of comfort orders land in-window -> ",
          "review the endpoint window before interpreting Aim 3 comfort sequencing.")
}
suppressWarnings(rm(cs_comfort))

## Comfort status at admission, and comfort-order timing descriptives ----------
# comfort_preexisting: the code status IN EFFECT AT ADMISSION is comfort-focused,
# i.e. the most recent order at or before admission_dttm is comfort-focused. A
# later Full or Limited order supersedes an earlier one, so taking the LAST order
# at or before admission is what "in effect" means. This separates patients who
# ARRIVED comfort-focused (standing POLST / AND directives) from those who
# converted during the stay -- the distinction the Aim 3 sequence table could not
# make when every encounter collapsed to "No comfort order".
# NOTE: 04 builds preexisting_code_status on a wider admission + 24h window for
# risk adjustment, so the two variables are related but will not agree exactly.
cs_adm = merge(
  code_status[, .(patient_id, start_dttm, cs_category)],
  ref[, .(hospitalization_id, patient_id, admission_dttm)],
  by = "patient_id", allow.cartesian = TRUE
)
cs_adm = cs_adm[start_dttm <= admission_dttm]
setorder(cs_adm, hospitalization_id, start_dttm)
cs_adm = cs_adm[, .SD[.N], by = hospitalization_id]   # the order in effect at admission
comfort_preexist = cs_adm[, .(hospitalization_id,
                              comfort_preexisting = as.integer(cs_category == "Comfort-focused"))]
suppressWarnings(rm(cs_adm))

# Encounter-level descriptives over the FULL cohort. Timing is measured on the
# first comfort order placed DURING the encounter; both quantities are signed so
# that positive means "hours after admission" and "hours before the endpoint".
# The endpoint is death where the patient died, otherwise discharge, which is the
# clinically meaningful anchor for an end-of-life cohort.
cdesc = ref[, .(hospitalization_id, admission_dttm, endpoint_dttm)]
cdesc = merge(cdesc, comfort_preexist, by = "hospitalization_id", all.x = TRUE)
cdesc[is.na(comfort_preexisting), comfort_preexisting := 0L]
cdesc = merge(cdesc, comfort_timing, by = "hospitalization_id", all.x = TRUE)
cdesc[, comfort_during := as.integer(!is.na(first_comfort_dttm))]
cdesc[, comfort_at_adm_or_during := as.integer(comfort_preexisting == 1L | comfort_during == 1L)]
cdesc[, hours_adm_to_comfort  := as.numeric(difftime(first_comfort_dttm, admission_dttm,
                                                     units = "hours"))]
cdesc[, hours_comfort_to_endpoint := as.numeric(difftime(endpoint_dttm, first_comfort_dttm,
                                                         units = "hours"))]

.cd_count = function(label, k, denom) data.table(
  metric = label, n = as.integer(k),
  pct    = if (is.na(denom) || denom == 0) NA_real_ else round(100 * k / denom, 2),
  mean = NA_real_, median = NA_real_, sd = NA_real_,
  q1 = NA_real_, q3 = NA_real_, iqr = NA_real_)

.cd_stat = function(label, x) {
  x = x[is.finite(x)]
  if (!length(x)) return(.cd_count(label, 0L, NA_real_))
  qq = stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
  data.table(
    metric = label, n = length(x), pct = NA_real_,
    mean   = round(mean(x), 1),
    median = round(stats::median(x), 1),
    sd     = if (length(x) > 1L) round(stats::sd(x), 1) else NA_real_,
    q1     = round(qq[[1]], 1), q3 = round(qq[[2]], 1),
    iqr    = round(qq[[2]] - qq[[1]], 1))
}

n_enc_cd = nrow(cdesc)
comfort_desc = rbindlist(list(
  .cd_count("Cohort encounters", n_enc_cd, n_enc_cd),
  .cd_count("Comfort status in effect at admission",
            sum(cdesc$comfort_preexisting), n_enc_cd),
  .cd_count("Comfort order placed during admission",
            sum(cdesc$comfort_during), n_enc_cd),
  .cd_count("Comfort status at admission or during admission",
            sum(cdesc$comfort_at_adm_or_during), n_enc_cd),
  .cd_stat("Hours from admission to first comfort order",  cdesc$hours_adm_to_comfort),
  .cd_stat("Hours from first comfort order to endpoint (death or discharge)",
           cdesc$hours_comfort_to_endpoint)
))
comfort_desc[, site := site_lowercase]

# Small-cell suppression: mask person counts below threshold, and where the
# contributing n is suppressed blank the distributional statistics too (an IQR
# computed on fewer than `small_cell_threshold` encounters is itself disclosive).
if (isTRUE(suppress_small_cells_enabled)) {
  comfort_desc = suppress_small_cells(comfort_desc, count_cols = "n", pct_cols = "pct")
  for (cc in c("mean", "median", "sd", "q1", "q3", "iqr"))
    comfort_desc[suppressed == 1L, (cc) := NA_real_]
}

fwrite(comfort_desc, shareable_path(
  paste0("comfort_status_descriptives_", site_lowercase, "_", today, ".csv")))
message("  Comfort status descriptives:")
for (i in seq_len(nrow(comfort_desc))) {
  r = comfort_desc[i]
  if (!is.na(r$pct)) {
    message(sprintf("    %-48s n=%s (%s%%)", r$metric,
                    format(r$n, big.mark = ","), r$pct))
  } else {
    message(sprintf("    %-48s n=%s mean=%s median=%s sd=%s IQR %s to %s",
                    r$metric, format(r$n, big.mark = ","),
                    r$mean, r$median, r$sd, r$q1, r$q3))
  }
}
message("  Saved comfort_status_descriptives.")
suppressWarnings(rm(cdesc, comfort_desc))

cs_out = merge(cs_out, comfort_timing, by = "hospitalization_id", all.x = TRUE)
cs_out = merge(cs_out, limited_timing, by = "hospitalization_id", all.x = TRUE)
cs_out = merge(cs_out, nonfull_timing, by = "hospitalization_id", all.x = TRUE)

rm(code_status, cs, cs_72h, all_cs, prior_full, cs_flag,
   comfort_timing, limited_timing, nonfull_timing)

n1 = cs_out[flag_code_status_strict == 1L, .N]
message(sprintf("    Strict (Full->non-Full): %s | Lenient (any non-Full): %s",
                format(n1, big.mark = ","),
                format(cs_out[flag_code_status_lenient == 1L, .N], big.mark = ",")))


# 2. opioid infusion - initiation or escalation in final 48h -------------------

message("  Component 2: Opioid infusion pattern...")

# opioid infusions are continuous medications. NOTE: intermittent opioid boluses
# are deliberately EXCLUDED pending a clinician-reviewed MME conversion scheme for
# combining infusion rates and bolus doses on a common (potency-adjusted) scale.
opioid_cats = c("morphine", "hydromorphone", "fentanyl")

med_cont = data_list[["medication_admin_continuous"]] |>
  dplyr::filter(tolower(med_category) %in% opioid_cats) |>
  dplyr::select(hospitalization_id, med_category, admin_dttm, med_dose, med_dose_unit) |>
  dplyr::collect() |>
  as.data.table()

med_cont = merge(med_cont, ref[, .(hospitalization_id, ep_minus_24h, ep_minus_48h,
                                    ep_minus_72h, endpoint_dttm, admission_dttm)],
                 by = "hospitalization_id")

# new opioid infusion in final 48h: no prior administration of that med
opioid_prior = med_cont[admin_dttm < ep_minus_48h,
                        .(had_prior = 1L),
                        by = .(hospitalization_id, med_category)]

opioid_final = med_cont[admin_dttm >= ep_minus_48h & admin_dttm <= endpoint_dttm]

opioid_final = merge(opioid_final, opioid_prior,
                     by = c("hospitalization_id", "med_category"), all.x = TRUE)
opioid_final[is.na(had_prior), had_prior := 0L]

# flag: initiation (new med in final 48h)
opioid_init = opioid_final[
  had_prior == 0L,
  .(flag_opioid_initiation = 1L),
  by = hospitalization_id
]

# flag: escalation (dose increase of 25% or more in final 48h for existing med)
opioid_esc_prior = med_cont[
  admin_dttm < ep_minus_48h & !is.na(as.numeric(med_dose)),
  .(prior_max_dose = max(as.numeric(med_dose), na.rm = TRUE)),
  by = .(hospitalization_id, med_category)
]

opioid_esc_final = opioid_final[
  !is.na(as.numeric(med_dose)),
  .(final_max_dose = max(as.numeric(med_dose), na.rm = TRUE)),
  by = .(hospitalization_id, med_category)
]

opioid_esc = merge(opioid_esc_final, opioid_esc_prior,
                   by = c("hospitalization_id", "med_category"))
opioid_esc = opioid_esc[
  final_max_dose > prior_max_dose * 1.25, # e25% increase
  .(flag_opioid_escalation = 1L),
  by = hospitalization_id
]

# combine: either initiation or escalation
opioid_out = merge(
  data.table(hospitalization_id = hids),
  opioid_init, by = "hospitalization_id", all.x = TRUE
)
opioid_out = merge(opioid_out, opioid_esc, by = "hospitalization_id", all.x = TRUE)
opioid_out[is.na(flag_opioid_initiation), flag_opioid_initiation := 0L]
opioid_out[is.na(flag_opioid_escalation), flag_opioid_escalation := 0L]
opioid_out[, flag_opioid := fifelse(flag_opioid_initiation == 1L | flag_opioid_escalation == 1L, 1L, 0L)]

# Aim 3: windowed opioid administration counts (72-48h, 48-24h, 24-0h)
opioid_windows = med_cont[, .(
  opioid_72_48h = sum(admin_dttm >= ep_minus_72h & admin_dttm < ep_minus_48h),
  opioid_48_24h = sum(admin_dttm >= ep_minus_48h & admin_dttm < ep_minus_24h),
  opioid_24_0h  = sum(admin_dttm >= ep_minus_24h & admin_dttm <= endpoint_dttm)
), by = hospitalization_id]

opioid_out = merge(opioid_out, opioid_windows, by = "hospitalization_id", all.x = TRUE)
for (opi_col in c("opioid_72_48h", "opioid_48_24h", "opioid_24_0h"))
  set(opioid_out, which(is.na(opioid_out[[opi_col]])), opi_col, 0L)

# Aim 3: opioid trajectory classification
opioid_any_prior = med_cont[admin_dttm < ep_minus_48h,
                            .(had_prior_opioid = 1L), by = hospitalization_id]
opioid_any_final = med_cont[admin_dttm >= ep_minus_48h & admin_dttm <= endpoint_dttm,
                            .(had_final_opioid = 1L), by = hospitalization_id]

opioid_out = merge(opioid_out, opioid_any_prior, by = "hospitalization_id", all.x = TRUE)
opioid_out = merge(opioid_out, opioid_any_final, by = "hospitalization_id", all.x = TRUE)
opioid_out[is.na(had_prior_opioid), had_prior_opioid := 0L]
opioid_out[is.na(had_final_opioid), had_final_opioid := 0L]

opioid_out[, opioid_trajectory := fcase(
  had_prior_opioid == 0L & had_final_opioid == 1L, "New start",
  had_prior_opioid == 1L & had_final_opioid == 1L, "Continuation",
  had_prior_opioid == 1L & had_final_opioid == 0L, "Stopped",
  default = "None"
)]
opioid_out[, c("had_prior_opioid", "had_final_opioid") := NULL]

rm(med_cont, opioid_prior, opioid_final, opioid_init, opioid_esc_prior,
   opioid_esc_final, opioid_esc, opioid_windows, opioid_any_prior, opioid_any_final)

message(sprintf("    Opioid initiation: %s | Escalation: %s | Either: %s",
                format(opioid_out[flag_opioid_initiation == 1L, .N], big.mark = ","),
                format(opioid_out[flag_opioid_escalation == 1L, .N], big.mark = ","),
                format(opioid_out[flag_opioid == 1L, .N], big.mark = ",")))


# 3. cessation of labs - no routine labs in final 48h --------------------------

message("  Component 3: Routine lab cessation...")

routine_labs = c("wbc", 
                 "creatinine", 
                 "platelet_count", 
                 "sodium", 
                 "potassium", 
                 "glucose")

labs = data_list[["labs"]] |>
  dplyr::filter(tolower(lab_category) %in% routine_labs) |>
  dplyr::select(hospitalization_id, lab_category, lab_result_dttm) |>
  dplyr::collect() |>
  as.data.table()

labs = merge(labs, ref[, .(hospitalization_id, 
                           admission_dttm,
                           ep_minus_24h,
                           ep_minus_48h,
                           ep_minus_72h,
                           endpoint_dttm)],
             by = "hospitalization_id")

# labs before the final 48h window (must have had labs earlier)
labs_early = labs[
  lab_result_dttm >= admission_dttm & lab_result_dttm < ep_minus_48h,
  .(had_early_labs = 1L),
  by = hospitalization_id
]

# labs in the final 48h window
labs_late = labs[
  lab_result_dttm >= ep_minus_48h & lab_result_dttm <= endpoint_dttm,
  .(had_late_labs = 1L),
  by = hospitalization_id
]

labs_flag = merge(data.table(hospitalization_id = hids), labs_early,
                  by = "hospitalization_id", all.x = TRUE)
labs_flag = merge(labs_flag, labs_late, by = "hospitalization_id", all.x = TRUE)
labs_flag[is.na(had_early_labs), had_early_labs := 0L]
labs_flag[is.na(had_late_labs), had_late_labs := 0L]

# cessation: had labs before but none in final 48h
labs_flag[, flag_lab_cessation := fifelse(had_early_labs == 1L & had_late_labs == 0L, 1L, 0L)]

labs_out = labs_flag[, .(hospitalization_id, flag_lab_cessation)]

# Aim 3: windowed lab counts (72-48h, 48-24h, 24-0h before endpoint)
lab_windows = labs[, .(
  labs_72_48h = sum(lab_result_dttm >= ep_minus_72h & lab_result_dttm < ep_minus_48h),
  labs_48_24h = sum(lab_result_dttm >= ep_minus_48h & lab_result_dttm < ep_minus_24h),
  labs_24_0h  = sum(lab_result_dttm >= ep_minus_24h & lab_result_dttm <= endpoint_dttm)
), by = hospitalization_id]

labs_out = merge(labs_out, lab_windows, by = "hospitalization_id", all.x = TRUE)
for (lab_col in c("labs_72_48h", "labs_48_24h", "labs_24_0h"))
  set(labs_out, which(is.na(labs_out[[lab_col]])), lab_col, 0L)

rm(labs, labs_early, labs_late, labs_flag, lab_windows)

message(sprintf("    Lab cessation (had early, none in final 48h): %s",
                format(labs_out[flag_lab_cessation == 1L, .N], big.mark = ",")))



# 4. withdrawal of LST ---------------------------------------------------------

message("  Component 4: Withdrawal of life-sustaining treatment...")


## 4a. vasopressor discontinuation in final 48h --------------------------------

message("    4a: Vasopressor discontinuation...")

vaso_cats = c("norepinephrine", "vasopressin", "epinephrine", "phenylephrine")

vaso = data_list[["medication_admin_continuous"]] |>
  dplyr::filter(tolower(med_category) %in% vaso_cats) |>
  dplyr::select(hospitalization_id, med_category, admin_dttm, med_dose) |>
  dplyr::collect() |>
  as.data.table()

vaso = merge(vaso, ref[, .(hospitalization_id, ep_minus_48h, endpoint_dttm)],
             by = "hospitalization_id")

# had vasopressors before final 48h
vaso_prior = vaso[
  admin_dttm < ep_minus_48h,
  .(had_vaso_prior = 1L),
  by = hospitalization_id
]

# any vasopressor administration in final 48h
vaso_late = vaso[
  admin_dttm >= ep_minus_48h & admin_dttm <= endpoint_dttm
]

# last vasopressor admin before endpoint
vaso_last = vaso_late[, .(last_vaso_dttm = max(admin_dttm)), by = hospitalization_id]

# discontinued if: had vasopressors, but last admin is >6h before endpoint

vaso_disc = merge(vaso_prior, vaso_last, by = "hospitalization_id", all.x = TRUE)
vaso_disc = merge(vaso_disc, ref[, .(hospitalization_id, endpoint_dttm)],
                  by = "hospitalization_id")

vaso_disc[, hours_to_endpoint := as.numeric(
  difftime(endpoint_dttm, last_vaso_dttm, units = "hours")
)]

# flag: vasopressors stopped >6h before endpoint, OR had prior vasopressors
# but none at all in final 48h

vaso_disc[, flag_vaso_withdrawal := fcase(
  had_vaso_prior == 1L & is.na(last_vaso_dttm), 1L,   # had prior, none in final 48h
  had_vaso_prior == 1L & hours_to_endpoint >= 6, 1L,   # stopped e6h before endpoint
  default = 0L
)]

vaso_out = vaso_disc[, .(hospitalization_id, flag_vaso_withdrawal)]

# Aim 3 timing: last vasopressor admin overall (for all patients who ever had vasos)
vaso_last_overall = vaso[, .(last_vaso_dttm = max(admin_dttm)), by = hospitalization_id]
vaso_out = merge(vaso_out, vaso_last_overall, by = "hospitalization_id", all.x = TRUE)

rm(vaso, vaso_prior, vaso_late, vaso_last, vaso_disc, vaso_last_overall)

message(sprintf("      Vasopressor withdrawal: %s",
                format(vaso_out[flag_vaso_withdrawal == 1L, .N], big.mark = ",")))


## 4b. ventilation withdrawal --------------------------------------------------

# Discontinued IMV without extubation to room air
# Logic: had IMV, then device_category changes to something other than
# "Room Air" (suggesting comfort/terminal extubation)
#
# Tracheostomy (§3.4.2): patients with a tracheostomy at first IMV are flagged
# as potentially chronically ventilated and analyzed separately — their LST
# trajectory (and the meaning of "ventilator cessation") differs from acutely
# intubated patients.

message("    4b: Mechanical ventilation withdrawal...")

resp_cols = c("hospitalization_id", "recorded_dttm", "device_category")
# include tracheostomy field if present in the CLIF respiratory_support table
resp_avail = tryCatch(names(data_list[["respiratory_support"]]), error = function(e) character(0))
has_trach_field = "tracheostomy" %in% resp_avail
if (has_trach_field) resp_cols = c(resp_cols, "tracheostomy")

resp = data_list[["respiratory_support"]] |>
  dplyr::select(dplyr::all_of(resp_cols)) |>
  dplyr::collect() |>
  as.data.table()

resp = merge(resp, ref[, .(hospitalization_id, ep_minus_48h, endpoint_dttm)],
             by = "hospitalization_id")

resp[, dev_lower := tolower(trimws(device_category))]

# had IMV during hospitalization
had_imv = resp[dev_lower == "imv", .(had_imv = 1L), by = hospitalization_id]

# tracheostomy at first IMV → chronic-vent flag
if (has_trach_field) {
  resp[, trach_val := suppressWarnings(as.integer(tracheostomy))]
  imv_first = resp[dev_lower == "imv"][order(recorded_dttm), .SD[1], by = hospitalization_id]
  trach_flag = imv_first[, .(hospitalization_id,
                             trach_at_first_imv = fifelse(!is.na(trach_val) & trach_val == 1L, 1L, 0L))]
  message(sprintf("      Tracheostomy at first IMV: %s",
                  format(trach_flag[trach_at_first_imv == 1L, .N], big.mark = ",")))
  rm(imv_first)
} else {
  message("      No tracheostomy field in respiratory_support — chronic-vent flag set to 0.")
  trach_flag = had_imv[, .(hospitalization_id, trach_at_first_imv = 0L)]
}

# in final 48h: last respiratory device
resp_final = resp[recorded_dttm >= ep_minus_48h & recorded_dttm <= endpoint_dttm]
resp_last = resp_final[order(-recorded_dttm), .SD[1], by = hospitalization_id]

# withdrawal: had IMV, last device is NOT IMV and NOT room air
# (room air = successful extubation; anything else = terminal wean or comfort)

mv_flag = merge(had_imv, resp_last[, .(hospitalization_id, last_device = dev_lower)],
                by = "hospitalization_id", all.x = TRUE)

mv_flag[, flag_mv_withdrawal := fcase(
  had_imv == 1L & !is.na(last_device) & last_device != "imv" & last_device != "room air", 1L,
  had_imv == 1L & is.na(last_device), 1L,  # had IMV but no resp record in final 48h
  default = 0L
)]

mv_out = mv_flag[, .(hospitalization_id, flag_mv_withdrawal)]
mv_out = merge(mv_out, trach_flag, by = "hospitalization_id", all.x = TRUE)
mv_out[is.na(trach_at_first_imv), trach_at_first_imv := 0L]

# Aim 3 timing: last IMV record overall
imv_last_overall = resp[dev_lower == "imv",
                        .(last_imv_dttm = max(recorded_dttm)),
                        by = hospitalization_id]
mv_out = merge(mv_out, imv_last_overall, by = "hospitalization_id", all.x = TRUE)

rm(resp, had_imv, resp_final, resp_last, mv_flag, imv_last_overall, trach_flag)

message(sprintf("      MV withdrawal: %s",
                format(mv_out[flag_mv_withdrawal == 1L, .N], big.mark = ",")))


## 4c. (RRT timing removed — limited cross-site CRRT coverage) -----------------

## 4d. combined withdrawal flag ------------------------------------------------

wdl_out = merge(data.table(hospitalization_id = hids), vaso_out,
                by = "hospitalization_id", all.x = TRUE)
wdl_out = merge(wdl_out, mv_out, by = "hospitalization_id", all.x = TRUE)

wdl_out[is.na(flag_vaso_withdrawal), flag_vaso_withdrawal := 0L]
wdl_out[is.na(flag_mv_withdrawal), flag_mv_withdrawal := 0L]

wdl_out[, flag_withdrawal := fifelse(
  flag_vaso_withdrawal == 1L | flag_mv_withdrawal == 1L, 1L, 0L
)]


## 4e. Aim 3 derived timing variables -----------------------------------------

# bring in code status timing (comfort, limited, and any-non-full anchors)
wdl_out = merge(wdl_out,
                cs_out[, .(hospitalization_id, first_comfort_dttm,
                           first_limited_dttm, first_nonfull_dttm)],
                by = "hospitalization_id", all.x = TRUE)

# first LST withdrawal datetime: earliest withdrawal across types (vaso/imv)
wdl_out[, first_lst_withdrawal_dttm := as.POSIXct(NA)]

wdl_out[flag_vaso_withdrawal == 1L & !is.na(last_vaso_dttm),
        first_lst_withdrawal_dttm := last_vaso_dttm]

wdl_out[flag_mv_withdrawal == 1L & !is.na(last_imv_dttm),
        first_lst_withdrawal_dttm := fifelse(
          is.na(first_lst_withdrawal_dttm) | last_imv_dttm < first_lst_withdrawal_dttm,
          last_imv_dttm, first_lst_withdrawal_dttm
        )]

# --- Comfort-contextualized vs uncontextualized LST cessation (§3.5.5) --------
# A cessation is "comfort-contextualized" if it occurs within ±24h of a
# transition to limited-resuscitation or comfort-focused code status. Otherwise
# it is "uncontextualized" (ambiguous: may be liberation, improvement, or a
# comfort transition not captured in code status). Aim 3 focuses on the former.

# the relevant code status anchor for contextualization: earliest of comfort or limited
wdl_out[, cs_context_dttm := fcase(
  !is.na(first_comfort_dttm) & !is.na(first_limited_dttm),
    pmin(first_comfort_dttm, first_limited_dttm),
  !is.na(first_comfort_dttm), first_comfort_dttm,
  !is.na(first_limited_dttm), first_limited_dttm,
  default = as.POSIXct(NA)
)]

wdl_out[, cessation_context := fcase(
  is.na(first_lst_withdrawal_dttm), "No cessation",
  is.na(cs_context_dttm), "Uncontextualized",
  abs(as.numeric(difftime(first_lst_withdrawal_dttm, cs_context_dttm, units = "hours"))) <= 24,
    "Comfort-contextualized",
  default = "Uncontextualized"
)]

wdl_out[, flag_comfort_contextualized_cessation := fifelse(
  cessation_context == "Comfort-contextualized", 1L, 0L
)]

# hours from withdrawal to hospice and death
wdl_out = merge(wdl_out,
                ref[, .(hospitalization_id, hospice_transition_dttm, endpoint_dttm)],
                by = "hospitalization_id", all.x = TRUE)

wdl_out[, hours_first_withdrawal_to_hospice := as.numeric(
  difftime(hospice_transition_dttm, first_lst_withdrawal_dttm, units = "hours")
)]
wdl_out[, hours_first_withdrawal_to_death := as.numeric(
  difftime(endpoint_dttm, first_lst_withdrawal_dttm, units = "hours")
)]

# hours from last LST (vaso/imv) to hospice (regardless of withdrawal flag)
wdl_out[, hours_vaso_to_hospice := as.numeric(
  difftime(hospice_transition_dttm, last_vaso_dttm, units = "hours")
)]
wdl_out[, hours_imv_to_hospice := as.numeric(
  difftime(hospice_transition_dttm, last_imv_dttm, units = "hours")
)]

# code status timing relative to hospice and withdrawal
wdl_out[, hours_comfort_to_hospice := as.numeric(
  difftime(hospice_transition_dttm, first_comfort_dttm, units = "hours")
)]
wdl_out[, hours_comfort_to_withdrawal := as.numeric(
  difftime(first_lst_withdrawal_dttm, first_comfort_dttm, units = "hours")
)]

# sequencing classification
wdl_out[, code_status_sequence := fcase(
  is.na(first_comfort_dttm), "No comfort order",
  is.na(first_lst_withdrawal_dttm), "Comfort only (no LST withdrawal)",
  abs(as.numeric(difftime(first_comfort_dttm, first_lst_withdrawal_dttm, units = "hours"))) < 2,
    "Simultaneous (<2h)",
  first_comfort_dttm < first_lst_withdrawal_dttm, "Comfort first",
  default = "Withdrawal first"
)]

# --- Clinical inflection point (§3.5.6) ---------------------------------------
# Stage 1 anchor = earliest of:
#   (a) Full -> non-Full code status transition (first_nonfull_dttm)
#   (b) comfort-contextualized LST cessation
# Where comfort-focused status is not captured, first_nonfull_dttm still anchors
# on the Full -> Limited Resuscitation transition (the §3.5.4 fallback).

# code-status component of the anchor: prefer comfort, fall back to any non-Full
wdl_out[, cs_anchor_dttm := fcase(
  !is.na(first_comfort_dttm), first_comfort_dttm,
  !is.na(first_nonfull_dttm), first_nonfull_dttm,
  default = as.POSIXct(NA)
)]

# LST component: only comfort-contextualized cessations qualify as an anchor
wdl_out[, lst_anchor_dttm := fifelse(
  flag_comfort_contextualized_cessation == 1L, first_lst_withdrawal_dttm, as.POSIXct(NA)
)]

wdl_out[, inflection_point_dttm := fcase(
  !is.na(cs_anchor_dttm) & !is.na(lst_anchor_dttm), pmin(cs_anchor_dttm, lst_anchor_dttm),
  !is.na(cs_anchor_dttm), cs_anchor_dttm,
  !is.na(lst_anchor_dttm), lst_anchor_dttm,
  default = as.POSIXct(NA)
)]

# record which signal anchored the inflection point
wdl_out[, inflection_anchor := fcase(
  is.na(inflection_point_dttm), "No observed inflection",
  !is.na(first_comfort_dttm) & inflection_point_dttm == first_comfort_dttm, "Comfort code status",
  !is.na(lst_anchor_dttm) & inflection_point_dttm == lst_anchor_dttm, "Comfort-contextualized cessation",
  default = "Full->Limited transition (fallback)"
)]

# hours from inflection point to hospice and death
wdl_out[, hours_inflection_to_hospice := as.numeric(
  difftime(hospice_transition_dttm, inflection_point_dttm, units = "hours")
)]
wdl_out[, hours_inflection_to_death := as.numeric(
  difftime(endpoint_dttm, inflection_point_dttm, units = "hours")
)]

n_inflection = sum(!is.na(wdl_out$inflection_point_dttm))
message(sprintf("    Inflection point identified: %s encounters (%4.1f%%)",
                format(n_inflection, big.mark = ","),
                100 * n_inflection / nrow(wdl_out)))
message(sprintf("    Comfort-contextualized cessations: %s | Uncontextualized: %s",
                format(wdl_out[cessation_context == "Comfort-contextualized", .N], big.mark = ","),
                format(wdl_out[cessation_context == "Uncontextualized", .N], big.mark = ",")))

# drop merge-duplicated and intermediate columns before final assembly
wdl_out[, c("hospice_transition_dttm", "endpoint_dttm", "first_comfort_dttm",
            "cs_context_dttm", "cs_anchor_dttm", "lst_anchor_dttm") := NULL]

message(sprintf("    Aim 3 timing computed: %s with any withdrawal, %s with first_lst_withdrawal_dttm.",
                format(sum(wdl_out$flag_withdrawal == 1L), big.mark = ","),
                format(sum(!is.na(wdl_out$first_lst_withdrawal_dttm)), big.mark = ",")))


# 5. look at broader medication discontinuation --------------------------------

# is this useful?

message("  Exploratory 5: Broad medication discontinuation...")

# any continuous med that was running and then stopped in final 48h
all_cont = data_list[["medication_admin_continuous"]] |>
  dplyr::select(hospitalization_id, med_category, admin_dttm) |>
  dplyr::collect() |>
  as.data.table()

all_cont = merge(all_cont, ref[, .(hospitalization_id, ep_minus_48h, endpoint_dttm)],
                 by = "hospitalization_id")

med_prior = all_cont[
  admin_dttm < ep_minus_48h,
  .(n_meds_prior = uniqueN(med_category)),
  by = hospitalization_id
]

med_late = all_cont[
  admin_dttm >= ep_minus_48h & admin_dttm <= endpoint_dttm,
  .(n_meds_late = uniqueN(med_category)),
  by = hospitalization_id
]

med_disc = merge(med_prior, med_late, by = "hospitalization_id", all.x = TRUE)
med_disc[is.na(n_meds_late), n_meds_late := 0L]

# flag: e50% reduction in unique continuous meds
med_disc[, flag_med_discontinuation := fifelse(
  n_meds_prior >= 2 & n_meds_late <= n_meds_prior * 0.5, 1L, 0L
)]

med_disc_out = med_disc[, .(hospitalization_id, flag_med_discontinuation,
                            n_meds_prior, n_meds_late)]

rm(all_cont, med_prior, med_late, med_disc)

message(sprintf("    e50%% reduction in continuous meds: %s",
                format(med_disc_out[flag_med_discontinuation == 1L, .N], big.mark = ",")))


# 6. decreased frequency of vital signs ----------------------------------------

message("  Exploratory 6: Diminished vital signs frequency...")

vitals = data_list[["vitals"]] |>
  dplyr::filter(vital_category %in% c("heart_rate", 
                                      "respiratory_rate", 
                                      "map", 
                                      "spo2", 
                                      "temp_c")) |>
  dplyr::select(hospitalization_id, recorded_dttm) |>
  dplyr::collect() |>
  as.data.table()

vitals = merge(vitals, ref[, .(hospitalization_id, 
                               ep_minus_48h, 
                               ep_minus_24h, 
                               endpoint_dttm)],
               by = "hospitalization_id")

# count vitals in the 24-48h window vs final 24h window
vitals_penult = vitals[
  recorded_dttm >= ep_minus_48h & recorded_dttm < ep_minus_24h,
  .(n_vitals_penult = .N),
  by = hospitalization_id
]

vitals_final = vitals[
  recorded_dttm >= ep_minus_24h & recorded_dttm <= endpoint_dttm,
  .(n_vitals_final = .N),
  by = hospitalization_id
]

vitals_freq = merge(vitals_penult, vitals_final, by = "hospitalization_id", all = TRUE)
vitals_freq[is.na(n_vitals_penult), n_vitals_penult := 0L]
vitals_freq[is.na(n_vitals_final), n_vitals_final := 0L]

# flag: e50% reduction in vital sign frequency
vitals_freq[, flag_vitals_diminished := fifelse(
  n_vitals_penult >= 4 & n_vitals_final <= n_vitals_penult * 0.5, 1L, 0L
)]

vitals_out = vitals_freq[, .(hospitalization_id, flag_vitals_diminished,
                             n_vitals_penult, n_vitals_final)]

rm(vitals, vitals_penult, vitals_final, vitals_freq)

message(sprintf("    e50%% reduction in vital frequency: %s",
                format(vitals_out[flag_vitals_diminished == 1L, .N], big.mark = ",")))


# 7. assemble clinical inflection --------------------------------------------------

message("  Assembling composite heuristic...")

cip = data.table(hospitalization_id = hids)

# Note: wdl_out already carries first_limited_dttm, first_nonfull_dttm (merged
# from cs_out in section 4e). Drop them from cs_out here to avoid .x/.y suffixes;
# keep only the flags and the comfort timing that wdl_out consumed and dropped.
cs_out_slim = cs_out[, .(hospitalization_id, flag_code_status_strict,
                         flag_code_status_lenient, first_comfort_dttm)]

cip = merge(cip, cs_out_slim,  by = "hospitalization_id", all.x = TRUE)
# merged separately from cs_out: an encounter whose ONLY comfort order predates
# admission never appears in cs_out (which is built from final-72h orders), and
# that is exactly the standing-order population this flag exists to capture.
cip = merge(cip, comfort_preexist, by = "hospitalization_id", all.x = TRUE)
cip = merge(cip, opioid_out,   by = "hospitalization_id", all.x = TRUE)
cip = merge(cip, labs_out,     by = "hospitalization_id", all.x = TRUE)
cip = merge(cip, wdl_out,      by = "hospitalization_id", all.x = TRUE)
cip = merge(cip, med_disc_out, by = "hospitalization_id", all.x = TRUE)
cip = merge(cip, vitals_out,   by = "hospitalization_id", all.x = TRUE)

# fill NAs with 0
flag_cols = grep("^flag_", names(cip), value = TRUE)
for (col in flag_cols) set(cip, which(is.na(cip[[col]])), col, 0L)
if ("comfort_preexisting" %in% names(cip))
  set(cip, which(is.na(cip$comfort_preexisting)), "comfort_preexisting", 0L)

# protocol components (4 core)
core_flags = c("flag_code_status_strict", "flag_opioid",
               "flag_lab_cessation", "flag_withdrawal")

cip[, n_core_components := rowSums(.SD), .SDcols = core_flags]

# composite: clinical inflection if e2 core components
cip[, inflection_2of4 := fifelse(n_core_components >= 2L, 1L, 0L)]

# alternative: e1 core component (more sensitive)
cip[, inflection_1of4 := fifelse(n_core_components >= 1L, 1L, 0L)]

# additive score across all flags (core + exploratory)
all_flags = c(core_flags, "flag_med_discontinuation", "flag_vitals_diminished")
cip[, n_all_components := rowSums(.SD), .SDcols = all_flags]


# 8. summmary ------------------------------------------------------------------


message("\n========== clinical inflection SUMMARY ==========")
message(sprintf("Primary cohort encounters:   %s", format(n_cohort, big.mark = ",")))
message("")
message("Individual component prevalence:")
for (f in c(core_flags, "flag_opioid_initiation", "flag_opioid_escalation",
            "flag_vaso_withdrawal", "flag_mv_withdrawal", 
            "flag_med_discontinuation", "flag_vitals_diminished",
            "flag_code_status_lenient")) {
  if (f %in% names(cip)) {
    n = sum(cip[[f]] == 1L, na.rm = TRUE)
    message(sprintf("  %-40s %s (%4.1f%%)", f, format(n, big.mark = ","),
                    100 * n / n_cohort))
  }
}

message("")
message("Composite heuristic:")
message(sprintf("  e1 of 4 core components:  %s (%4.1f%%)",
                format(sum(cip$inflection_1of4), big.mark = ","),
                100 * sum(cip$inflection_1of4) / n_cohort))
message(sprintf("  e2 of 4 core components:  %s (%4.1f%%)",
                format(sum(cip$inflection_2of4), big.mark = ","),
                100 * sum(cip$inflection_2of4) / n_cohort))

message("")
message("Core component overlap (among those with e1):")
overlap = cip[inflection_1of4 == 1L]
if (nrow(overlap) > 0) {
  message(sprintf("  n with e1 component: %s", format(nrow(overlap), big.mark = ",")))
  for (f in core_flags) {
    message(sprintf("    %-38s %s (%4.1f%% of flagged)",
                    f, format(sum(overlap[[f]]), big.mark = ","),
                    100 * sum(overlap[[f]]) / nrow(overlap)))
  }
  # pairwise overlap
  message("")
  message("  Pairwise overlap (% of flagged encounters with BOTH):")
  pairs = combn(core_flags, 2, simplify = FALSE)
  for (pair in pairs) {
    n_both = sum(overlap[[pair[1]]] == 1L & overlap[[pair[2]]] == 1L)
    message(sprintf("    %-30s & %-30s %s (%4.1f%%)",
                    gsub("flag_", "", pair[1]),
                    gsub("flag_", "", pair[2]),
                    format(n_both, big.mark = ","),
                    100 * n_both / nrow(overlap)))
  }
}

# exploratory overlap with core
message("")
message("Exploratory flags among those clinical inflection (e2 core):")
rd = cip[inflection_2of4 == 1L]
if (nrow(rd) > 0) {
  for (f in c("flag_med_discontinuation", "flag_vitals_diminished", "flag_lab_cessation")) {
    message(sprintf("  %-40s %s (%4.1f%%)",
                    f, format(sum(rd[[f]]), big.mark = ","),
                    100 * sum(rd[[f]]) / nrow(rd)))
  }
}
message("=====================================================\n")

rm(overlap, rd)


# 9. save full-cohort output ----------------------------------------------------

save_output(cip, "clinical_inflection_flags", shareable = FALSE)

save_output(
  cip[, .(
    hospitalization_id,
    flag_code_status_strict, flag_code_status_lenient,
    first_comfort_dttm,
    flag_opioid, flag_opioid_initiation, flag_opioid_escalation,
    opioid_72_48h, opioid_48_24h, opioid_24_0h, opioid_trajectory,
    flag_lab_cessation,
    labs_72_48h, labs_48_24h, labs_24_0h,
    flag_withdrawal, flag_vaso_withdrawal, flag_mv_withdrawal,
    trach_at_first_imv,
    last_vaso_dttm, last_imv_dttm,
    first_lst_withdrawal_dttm,
    cessation_context, flag_comfort_contextualized_cessation,
    hours_first_withdrawal_to_hospice, hours_first_withdrawal_to_death,
    hours_vaso_to_hospice, hours_imv_to_hospice,
    hours_comfort_to_hospice, hours_comfort_to_withdrawal,
    first_limited_dttm, first_nonfull_dttm,
    code_status_sequence,
    inflection_point_dttm, inflection_anchor,
    hours_inflection_to_hospice, hours_inflection_to_death,
    flag_med_discontinuation, flag_vitals_diminished,
    n_core_components, inflection_2of4, inflection_1of4,
    n_all_components
  )],
  name      = "clinical_inflection",
  shareable = FALSE
)


# 10. subgroup-level summaries (shareable) -------------------------------------

# Runs the same summary statistics for the full cohort and each LST subgroup.
# Output: one CSV with all groups stacked, for cross-site compilation.

message("\n  Building subgroup clinical inflection summaries...")

# merge LST flags + outcome from cohort onto clinical inflection
lst_cols = intersect(
  c("hospitalization_id", "imv_flag", "vaso_flag",
    "lst_any_flag", "outcome_category"),
  names(cohort_all)
)
inflection_sg = merge(cip, cohort_all[in_cohort == 1L, ..lst_cols],
                 by = "hospitalization_id", all.x = TRUE)

# define subgroups
sg_defs = list(
  list(label = "full_cohort",  filter_expr = quote(TRUE)),
  list(label = "lst_any",      filter_expr = quote(lst_any_flag == 1L)),
  list(label = "lst_imv",      filter_expr = quote(imv_flag == 1L)),
  list(label = "lst_vaso",     filter_expr = quote(vaso_flag == 1L))
)

# summary function: component prevalence + composite for one subgroup
summarize_heuristic = function(dt, label) {
  n = nrow(dt)
  if (n == 0) return(NULL)

  flag_names = c(
    "flag_code_status_strict", "flag_code_status_lenient",
    "flag_opioid", "flag_opioid_initiation", "flag_opioid_escalation",
    "flag_lab_cessation",
    "flag_withdrawal", "flag_vaso_withdrawal", "flag_mv_withdrawal",
    
    "flag_med_discontinuation", "flag_vitals_diminished",
    "inflection_1of4", "inflection_2of4"
  )
  flag_names = intersect(flag_names, names(dt))

  rows = rbindlist(lapply(flag_names, \(f) {
    n_flagged = sum(dt[[f]] == 1L, na.rm = TRUE)
    data.table(
      component = f,
      n_flagged = n_flagged,
      pct       = round(100 * n_flagged / n, 2)
    )
  }))

  rows[, `:=`(subgroup = label, site = site_lowercase, n_subgroup = n)]
  setcolorder(rows, c("site", "subgroup", "n_subgroup", "component",
                       "n_flagged", "pct"))
  rows
}

heuristic_summary = rbindlist(lapply(sg_defs, \(sg) {
  sub = inflection_sg[eval(sg$filter_expr)]
  summarize_heuristic(sub, sg$label)
}))

# print to console
message("\n========== clinical inflection BY SUBGROUP ==========")
for (sg in unique(heuristic_summary$subgroup)) {
  sub = heuristic_summary[subgroup == sg]
  n_sub = sub$n_subgroup[1]
  message(sprintf("\n  %s (n = %s):", sg, format(n_sub, big.mark = ",")))
  for (i in seq_len(nrow(sub))) {
    message(sprintf("    %-40s %s (%5.2f%%)",
                    sub$component[i],
                    format(sub$n_flagged[i], big.mark = ","),
                    sub$pct[i]))
  }
}
message("==================================================\n")

# save for cross-site compilation
fwrite(
  heuristic_summary,
  shareable_path(
            paste0("clinical_inflection_by_subgroup_", site_lowercase, "_", today, ".csv"))
)
message("Subgroup clinical inflection summary saved to shareable output.")

rm(cs_out, opioid_out, labs_out, wdl_out, vaso_out, mv_out,
   med_disc_out, vitals_out, ref, flag_cols, core_flags,
   all_flags,
   inflection_sg, sg_defs, heuristic_summary)

gc()

message("02_clinical_inflection.R complete.")
