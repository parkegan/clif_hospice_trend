
# 01_cohort.R

# Primary analysis cohort for hospice discharge study

# Protocol v5.0 Sections 3.4–3.5

#   - Study period: 2016-01-01 through 2024-12-31
#   - Adults (≥18) with inpatient admissions
#   - Exclusions:
#       * Missing or invalid length of stay (missing dates, admit > discharge)
#       * Admitted from hospice (first ADT location = hospice)
#   - Exclusion tracking table exported to final_output/ for cross-site compilation
#   - Discharge disposition as multinomial outcome:
#       1. Routine discharge (home, SNF, LTAC, rehab — non-hospice, alive)
#       2. In-hospital death
#       3. Same-facility hospice → death
#       4. Same-facility hospice → alive (at 30 days)
#       5. External hospice discharge
#       6. Transfer to acute facility
#       7. Other / AMA

# Transfers to other acute care facilities are INCLUDED in the
# primary cohort to avoid collider bias (see protocol §1.7).


# Last updated bpe 27.05.2026


# 0. source setup --------------------------------------------------------------

source(file.path(here::here("code"), "00_setup.R"))


# 1. Adult pts within study window ---------------------------------------------

hosp = 
  data_list[["hospitalization"]] |>
  dplyr::filter(
    admission_dttm >= start_date,
    admission_dttm <= end_date,
    age_at_admission >= 18L
  ) |>
  dplyr::select(
    patient_id, hospitalization_id,
    admission_dttm, discharge_dttm,
    age_at_admission,
    admission_type_name, admission_type_category,
    discharge_name, discharge_category
  ) |>
  dplyr::collect() |>
  as.data.table()

n_total = nrow(hosp)
message(sprintf("Total adult hospitalizations in study window: %s", 
                format(n_total, big.mark = ",")))


## 1a. exclusion tracking — missing or invalid length of stay -----------------

# missing admission or discharge datetime → no calculable LOS
n_before_los = nrow(hosp)
hosp = hosp[!is.na(admission_dttm) & !is.na(discharge_dttm)]
n_excl_missing_dates = n_before_los - nrow(hosp)

# admission after discharge → invalid LOS
n_before_order = nrow(hosp)
hosp = hosp[admission_dttm <= discharge_dttm]
n_excl_date_order = n_before_order - nrow(hosp)

n_excl_los = n_excl_missing_dates + n_excl_date_order

n_base = nrow(hosp)
message(sprintf("After LOS exclusions: %s (removed %s missing dates, %s admit > discharge)", 
                format(n_base, big.mark = ","),
                format(n_excl_missing_dates, big.mark = ","),
                format(n_excl_date_order, big.mark = ",")))


# 2. ADT — identify hospice location and first ADT location -------------------

# ADT for relevant encounters
hosp_ids_arrow = 
  arrow::arrow_table(
  data.frame(hospitalization_id = unique(hosp$hospitalization_id))
)

adt = 
  data_list[["adt"]] |>
  dplyr::semi_join(hosp_ids_arrow, by = "hospitalization_id") |>
  dplyr::select(hospitalization_id, 
                hospital_id, 
                location_category, 
                in_dttm, 
                out_dttm) |>
  dplyr::collect() |>
  as.data.table()

rm(hosp_ids_arrow)


## 2a. ICU admission flag covariate --------------------------------------------
# we will need for aim 2/3 analyses and sensitivity cohorts

adt_icu = adt[tolower(location_category) == "icu"]

adt_icu = merge(
  adt_icu,
  hosp[, .(hospitalization_id, admission_dttm)],
  by = "hospitalization_id"
)

adt_icu[, hours_to_icu := as.numeric(difftime(in_dttm, admission_dttm, units = "hours"))]

# any ICU admission during hospitalization
icu_ever = adt_icu[
  hours_to_icu >= 0,
  .(icu_flag = 1L,
    icu_48h_flag = as.integer(any(hours_to_icu <= 48))),
  by = hospitalization_id
]

hosp = merge(hosp, icu_ever, by = "hospitalization_id", all.x = TRUE)
hosp[is.na(icu_flag), icu_flag := 0L]
hosp[is.na(icu_48h_flag), icu_48h_flag := 0L]

rm(adt_icu, icu_ever)

n_icu = sum(hosp$icu_flag == 1L)
message(sprintf("  Any ICU admission: %s", format(n_icu, big.mark = ",")))


## 2b. Any hospice ADT location -----------------------------------------------
# Used for same-facility hospice transition detection

# Note: OHSU stopped their inpatient hospice in 2022(?) 
# so will need buddy code site to have an inpatient hospice for checking this work

adt_hospice = adt[
  tolower(location_category) == "hospice",
  .(hospice_location_flag = 1L,
    hospice_in_dttm  = min(in_dttm),
    hospice_out_dttm = max(out_dttm)),
  by = hospitalization_id
]

# Hospital-level inpatient hospice unit flag (protocol v6.0 §3.5.2):
# a hospital has an inpatient hospice unit if ANY ADT record at that hospital
# has location_category = "hospice" during the study period. Used to restrict
# same-facility hospice variation analyses in Aim 2.
hosp_with_hospice_unit = adt[
  tolower(location_category) == "hospice",
  unique(hospital_id)
]
message(sprintf("  Hospitals with inpatient hospice unit: %d",
                length(hosp_with_hospice_unit)))


## 2c. first ADT location (exclude admissions straight to hospice) -------------

adt_first = adt[
  order(in_dttm),
  .SD[1],
  by = hospitalization_id
][, .(hospitalization_id, first_location = tolower(location_category))]


## 2d. hospital ID (for hospital-level analyses) -------------------------------

# use first hospital id (not currently accounting for readmissions)

# Note: There are some CLIF sites which are listed as individual facilities
# but per CMS are one big facility! This is taken into account in the 01 script


hospital_map = adt[
  !is.na(hospital_id) & nzchar(as.character(hospital_id))
][order(in_dttm),
  .(hospital_id = first(hospital_id)),
  by = hospitalization_id]

hosp = merge(hosp, hospital_map, by = "hospitalization_id", all.x = TRUE)

# Data-quality exclusion: an encounter with no resolvable hospital_id cannot be
# placed at a facility, so it cannot enter hospital-level (Aim 2) or CCN-level
# shareable analyses. Exclude these outright. NA also covers empty-string ids,
# which the hospital_map filter above already dropped via nzchar().
# TODO: revisit once the new ADT table lands — hospital_id should always populate.
n_before_hid  = nrow(hosp)
n_missing_hid = sum(is.na(hosp$hospital_id))
if (n_missing_hid > 0) {
  qi_flag(
    "Missing ADT hospital_id",
    sprintf(paste0(
      "%s of %s encounters (%.1f%%) have no ADT hospital_id and are being dropped\n",
      "(they cannot be placed at a facility, so cannot enter hospital- or CCN-level\n",
      "analyses). A non-trivial share usually means the ADT/location table is\n",
      "incomplete at this site -- review with your data team before trusting\n",
      "site- or hospital-level results."),
      format(n_missing_hid, big.mark = ","),
      format(n_before_hid,  big.mark = ","),
      100 * n_missing_hid / n_before_hid))

  # Breakdown: is the gap random or concentrated? Computed before the exclusion,
  # while hosp still holds the missing-hospital_id rows. A roughly even rate
  # across years/types points to a random coverage gap; a spike in one year or
  # encounter type points to a systematic extract problem (and possible bias).
  by_year = hosp[, .(encounters = .N, missing = sum(is.na(hospital_id)),
                     pct = round(100 * mean(is.na(hospital_id)), 1)),
                 by = .(admit_year = data.table::year(admission_dttm))][order(admit_year)]
  by_type = hosp[, .(encounters = .N, missing = sum(is.na(hospital_id)),
                     pct = round(100 * mean(is.na(hospital_id)), 1)),
                 by = .(admission_type_category)][order(-pct)]
  message("    Missing-hospital_id rate by admission year:")
  for (ln in capture.output(print(by_year, row.names = FALSE))) message("      ", ln)
  message("    Missing-hospital_id rate by admission type:")
  for (ln in capture.output(print(by_type, row.names = FALSE))) message("      ", ln)
  message("    (Even across years/types => random coverage gap; a spike => systematic.)")
  rm(by_year, by_type)
}

hosp = hosp[!is.na(hospital_id)]
n_excl_missing_hid = n_before_hid - nrow(hosp)

rm(hospital_map)



## 2e. ED admission flag covariate ---------------------------------------------

adt_ed = adt[
  tolower(location_category) == "ed",
  .(ed_flag = 1L),
  by = hospitalization_id
]

hosp = merge(hosp, adt_ed, by = "hospitalization_id", all.x = TRUE)
hosp[is.na(ed_flag), ed_flag := 0L]

rm(adt_ed)


## 2f. Life-sustaining therapy (LST) flags (protocol §3.4.2) ------------------

# Identifies patients who received IMV or vasopressors during their
# hospitalization. These flags define the Aim 3 LST subcohort and are used
# to run subgroup analyses across the pipeline.
#
# NOTE: RRT (renal replacement therapy) is NOT included as an LST type in this
# study. Only a minority of consortium sites have the crrt_therapy table CLIF'd,
# so RRT is deferred to a future analysis when coverage is broader.

message("  Identifying LST subcohort...")

# --- IMV: invasive mechanical ventilation ---
imv_ids = data_list[["respiratory_support"]] |>
  dplyr::filter(tolower(device_category) == "imv") |>
  dplyr::select(hospitalization_id) |>
  dplyr::distinct() |>
  dplyr::collect() |>
  as.data.table()
imv_ids[, imv_flag := 1L]

hosp = merge(hosp, imv_ids, by = "hospitalization_id", all.x = TRUE)
hosp[is.na(imv_flag), imv_flag := 0L]
rm(imv_ids)

# --- Vasopressors ---
vaso_cats = c("norepinephrine", "vasopressin", "epinephrine",
              "phenylephrine", "dopamine")

vaso_ids = data_list[["medication_admin_continuous"]] |>
  dplyr::filter(tolower(med_category) %in% vaso_cats) |>
  dplyr::select(hospitalization_id) |>
  dplyr::distinct() |>
  dplyr::collect() |>
  as.data.table()
vaso_ids[, vaso_flag := 1L]

hosp = merge(hosp, vaso_ids, by = "hospitalization_id", all.x = TRUE)
hosp[is.na(vaso_flag), vaso_flag := 0L]
rm(vaso_ids, vaso_cats)

# --- Any LST (union of IMV and vasopressors) ---
hosp[, lst_any_flag := fifelse(
  imv_flag == 1L | vaso_flag == 1L, 1L, 0L
)]

n_imv  = sum(hosp$imv_flag  == 1L)
n_vaso = sum(hosp$vaso_flag == 1L)
n_lst  = sum(hosp$lst_any_flag == 1L)

message(sprintf("  LST flags: IMV=%s | Vaso=%s | Any LST=%s",
                format(n_imv, big.mark = ","),
                format(n_vaso, big.mark = ","),
                format(n_lst, big.mark = ",")))


# 3. exclusion criteria --------------------------------------------------------

# Protocol 3.4.1: exclude patients admitted from hospice

hosp = merge(hosp, adt_first, by = "hospitalization_id", all.x = TRUE)

hosp[, excl_hospice_admit := fifelse(
  first_location == "hospice", 1L, 0L, na = 0L
)]

# 24-hour re-admit exception (§3.4.1): a patient discharged to hospice, then
# admitted to the SAME facility's hospice within 24h, is a same-facility
# transition (continuation of the same episode), NOT an admitted-from-hospice
# patient. Identify these by checking each hospice-admit encounter for a prior
# encounter at the same hospital with hospice discharge ending within 24h.

hospice_admit_ids = hosp[excl_hospice_admit == 1L, hospitalization_id]

if (length(hospice_admit_ids) > 0) {
  # candidate prior encounters: same patient, hospice discharge
  prior_disch = hosp[
    tolower(trimws(discharge_category)) == "hospice",
    .(patient_id, hospital_id, prior_discharge_dttm = discharge_dttm,
      prior_hid = hospitalization_id)
  ]

  # this-encounter admit info for the hospice admits
  this_admit = hosp[excl_hospice_admit == 1L,
                    .(hospitalization_id, patient_id, hospital_id, admission_dttm)]

  # join on patient + hospital, find prior discharge within 24h before admit
  readmit_check = merge(this_admit, prior_disch,
                        by = c("patient_id", "hospital_id"),
                        allow.cartesian = TRUE)
  readmit_check = readmit_check[
    prior_hid != hospitalization_id &
      prior_discharge_dttm <= admission_dttm &
      as.numeric(difftime(admission_dttm, prior_discharge_dttm, units = "hours")) <= 24
  ]

  reclassify_ids = unique(readmit_check$hospitalization_id)

  if (length(reclassify_ids) > 0) {
    hosp[hospitalization_id %in% reclassify_ids, excl_hospice_admit := 0L]
    message(sprintf("  24h re-admit exception: %s encounters reclassified as same-facility transitions",
                    format(length(reclassify_ids), big.mark = ",")))
  }
  rm(prior_disch, this_admit, readmit_check)
}

rm(adt_first)

n_excl_hospice = sum(hosp$excl_hospice_admit == 1L)

hosp[, in_cohort := fifelse(excl_hospice_admit == 0L, 1L, 0L)]

n_cohort = sum(hosp$in_cohort == 1L)

message(sprintf("  Excluded — admitted from hospice: %s", format(n_excl_hospice, big.mark = ",")))


# ---- minimum hospital contribution -------------------------------------------
# Hospitals contributing fewer than min_hospital_encounters to the cohort are
# removed from it entirely. Two reasons. Statistically, a hospital of a few dozen
# encounters still receives a random intercept in 04, along with a
# risk-standardized probability, a reliability, and a rank interval, and it
# inflates apparent between-hospital variation while carrying almost no
# information. For disclosure, per-hospital rates on a denominator that small are
# what cell suppression exists to prevent.
#
# This is a NETWORK constant. It must be identical at every site and is not a
# site-tunable setting.
#
# Applied at hospital_id grain, because that is the unit at which a facility
# contributes data and bears disclosure risk. Where several hospital_ids later
# collapse to one CCN (§8b), the threshold has already been applied to each
# component separately, which is the conservative direction.
min_hospital_encounters = 1000L

hosp[, excl_small_hospital := 0L]

hosp_contrib = hosp[in_cohort == 1L, .(n_encounters = .N),
                    by = hospital_id][order(n_encounters)]
small_hospitals = hosp_contrib[n_encounters < min_hospital_encounters, hospital_id]

if (length(small_hospitals) > 0) {
  hosp[in_cohort == 1L & hospital_id %in% small_hospitals, excl_small_hospital := 1L]
  hosp[excl_small_hospital == 1L, in_cohort := 0L]

  message(sprintf("  Excluded — %d hospital(s) contributing < %s cohort encounters:",
                  length(small_hospitals), format(min_hospital_encounters, big.mark = ",")))
  drop_tbl = hosp_contrib[n_encounters < min_hospital_encounters]
  for (i in seq_len(nrow(drop_tbl)))
    message(sprintf("      %-24s %s encounters", drop_tbl$hospital_id[i],
                    format(drop_tbl$n_encounters[i], big.mark = ",")))

  # named list is site-local and stays local: hospital_id is a facility name at
  # some sites, and every shareable output in this study is keyed on CCN.
  # The encounter count reaches the network through the exclusion table below.
  save_output(cbind(site = site_lowercase, drop_tbl,
                    threshold = min_hospital_encounters),
              "hospitals_excluded_small", shareable = FALSE)
  rm(drop_tbl)
}

n_excl_small_hosp = sum(hosp$excl_small_hospital == 1L)
n_cohort          = sum(hosp$in_cohort == 1L)
rm(hosp_contrib, small_hospitals)

message(sprintf("  Excluded — encounters at those hospitals: %s",
                format(n_excl_small_hosp, big.mark = ",")))
message(sprintf("PRIMARY COHORT: %s hospitalizations across %d hospitals",
                format(n_cohort, big.mark = ","),
                hosp[in_cohort == 1L, uniqueN(hospital_id)]))


# 4. Outcomes ------------------------------------------------------------------

# Protocol 3.5.2 multinomial outcomes:

#   1. Routine discharge (home, SNF, LTAC, rehab — non-hospice, alive)
#   2. In-hospital death (no hospice ADT)
#   3. Same-facility hospice → death (transitioned to hospice, died)
#   4. Same-facility hospice → alive at 30d (is this useful since most sites won't have state death data?)
#   5. External hospice discharge (disposition = hospice, no hospice ADT)
#   6. Transfer to another acute care facility
#   7. Other / AMA


## 4a. pull death dates from patient table -------------------------------------

# Note: OHSU and Minnesota have state death data that will capture outcomes
# after discharge. Other sites do not have this currently.

patient = data_list[["patient"]] |>
  dplyr::select(patient_id, death_dttm) |>
  dplyr::collect() |>
  as.data.table()

hosp = merge(hosp, patient, by = "patient_id", all.x = TRUE)
rm(patient)

# Death is the terminal event of an encounter: nothing in a stay can occur after
# the patient has died. A death timestamp that PRECEDES discharge means the
# discharge timestamp is spurious -- most often an organ-donor procurement
# admission recorded after death, or administrative discharge paperwork logged
# after an in-hospital death. Apply the earliest-timestamp rule: clamp discharge
# to the death time rather than dropping the row.
n_death_before_discharge = hosp[!is.na(death_dttm) & death_dttm < discharge_dttm, .N]
if (n_death_before_discharge > 0) {
  qi_flag(
    "Death precedes discharge -- discharge clamped to death",
    sprintf(paste0(
      "%s encounter(s): death_dttm is earlier than discharge_dttm. Death is terminal,\n",
      "so discharge is set to the death time (earliest-timestamp rule). This pattern\n",
      "comes from organ-donor readmissions recorded after death and from discharge\n",
      "paperwork timestamped after an in-hospital death. Review death-date hygiene if\n",
      "this count is large relative to the cohort."),
      format(n_death_before_discharge, big.mark = ",")))
  hosp[!is.na(death_dttm) & death_dttm < discharge_dttm, discharge_dttm := death_dttm]

  # an encounter whose clamped discharge now precedes admission was recorded
  # ENTIRELY after death (a post-mortem organ-donor admission); not a valid
  # inpatient stay for this cohort.
  n_post_mortem = hosp[discharge_dttm < admission_dttm, .N]
  if (n_post_mortem > 0) {
    message(sprintf("    ...of which %s were entirely post-mortem admissions; dropped.",
                    format(n_post_mortem, big.mark = ",")))
    hosp = hosp[discharge_dttm >= admission_dttm]
  }
}

# compute days from discharge to death (>= 0 after the clamp above)
hosp[, days_discharge_to_death := as.numeric(
  difftime(death_dttm, discharge_dttm, units = "days")
)]


## 4b. find any hospice ADT flags ----------------------------------------------

hosp = merge(hosp, adt_hospice, by = "hospitalization_id", all.x = TRUE)
hosp[is.na(hospice_location_flag), hospice_location_flag := 0L]

rm(adt_hospice, adt)


## 4c. classify outcomes -------------------------------------------------------

hosp[, discharge_cat_lower := tolower(trimws(discharge_category))]

hosp[, outcome_category := fcase(
  
  # 2. In-hospital death on acute status (no hospice location during stay)
  discharge_cat_lower == "expired" & hospice_location_flag == 0L,
  "in_hospital_death_acute",
  
  # 3. In-hospital death with hospice location → same-facility hospice death
  discharge_cat_lower == "expired" & hospice_location_flag == 1L,
  "same_facility_hospice_death",
  
  # 3. Discharged to hospice, had hospice ADT location (same-facility), died ≤30d
  discharge_cat_lower == "hospice" & hospice_location_flag == 1L &
    !is.na(death_dttm) & days_discharge_to_death >= 0 & days_discharge_to_death <= 30,
  "same_facility_hospice_death",
  
  # 4. Same-facility hospice, alive at 30 days or no death recorded
  discharge_cat_lower == "hospice" & hospice_location_flag == 1L &
    (is.na(death_dttm) | days_discharge_to_death > 30 | days_discharge_to_death < 0),
  "same_facility_hospice_alive",
  
  # 5. External hospice discharge (hospice disposition, no hospice ADT location)
  discharge_cat_lower == "hospice" & hospice_location_flag == 0L,
  "external_hospice_discharge",
  
  # 6. Transfer to another acute care facility. CLIF 2.1.0 codes hospital-to-
  #    hospital transfers as discharge_category = "Acute Care Hospital" (which
  #    contains no "transfer" substring); the grepl() is retained as a safety net
  #    for sites whose discharge_category still holds a raw transfer string.
  discharge_cat_lower == "acute care hospital" | grepl("transfer", discharge_cat_lower),
  "transfer_acute_facility",
  
  # 7. AMA / left against medical advice
  grepl("ama|against medical advice|left against", discharge_cat_lower),
  "other_ama",
  
  # 1. Routine discharge: alive, non-hospice, non-transfer
  discharge_cat_lower != "expired" & discharge_cat_lower != "hospice" &
    discharge_cat_lower != "acute care hospital" &
    !grepl("transfer", discharge_cat_lower) &
    !grepl("ama|against medical advice", discharge_cat_lower),
  "routine_discharge",
  
  # fallback
  default = "other_ama"
)]


## 4d. time from hospice transition to death (hours) --------------------------

# For same-facility hospice: use hospice ADT in_dttm
# For external hospice: use discharge_dttm as proxy for transition

hosp[, hospice_transition_dttm := fcase(
  hospice_location_flag == 1L, hospice_in_dttm,
  discharge_cat_lower == "hospice", discharge_dttm,
  default = as.POSIXct(NA)
)]

# a hospice transition cannot post-date death either -- same earliest-timestamp rule
hosp[!is.na(death_dttm) & !is.na(hospice_transition_dttm) &
       hospice_transition_dttm > death_dttm, hospice_transition_dttm := death_dttm]

hosp[, hours_transition_to_death := as.numeric(
  difftime(death_dttm, hospice_transition_dttm, units = "hours")
)]

# After the earliest-timestamp clamps above, no negative death interval should
# remain. This residual check catches anything the clamp did not (e.g. a death or
# transition timestamp that is itself corrupt), which would signal a parsing or
# timezone problem rather than the benign post-death-admission pattern.
n_neg_dtd = hosp[!is.na(days_discharge_to_death)  & days_discharge_to_death  < 0, .N]
n_neg_ttd = hosp[!is.na(hours_transition_to_death) & hours_transition_to_death < 0, .N]
if (n_neg_dtd > 0 || n_neg_ttd > 0) {
  qi_flag(
    "Residual negative death intervals after clamping",
    sprintf(paste0(
      "%s encounter(s) still show death before discharge and %s before the hospice\n",
      "transition AFTER the earliest-timestamp clamp. This is unexpected -- inspect\n",
      "death_dttm / discharge_dttm parsing and timezone handling before trusting\n",
      "Analysis 1.5 and the Aim 1 competing-risks outputs."),
      format(n_neg_dtd, big.mark = ","),
      format(n_neg_ttd, big.mark = ",")))
}

# protocol 3.5.3: short-survival = death within 48h of hospice transition
hosp[, short_interval_death := fifelse(
  !is.na(hours_transition_to_death) & hours_transition_to_death >= 0 &
    hours_transition_to_death <= 48,
  1L, 0L, na = 0L
)]


## 4e. convenience covariates (for downstream analyses) -----------------------

hosp[, age_65_flag := fifelse(age_at_admission >= 65L, 1L, 0L)]

hosp[, `:=`(
  admission_month  = floor_date(admission_dttm, "month"),
  discharge_month  = floor_date(discharge_dttm, "month"),
  los_days         = as.numeric(difftime(discharge_dttm, admission_dttm, units = "days"))
)]


# 5. cohort summary -----------------------------------------------------------

# hospital- and encounter-level flags must be set on hosp BEFORE splitting off
# the cohort, so both the full file and the primary cohort carry them.

# death-linked site flag (OHSU and UMN have state death registry data
# incorporated into their CLIF patient death_dttm, enabling external hospice
# death capture and sensitivity analyses per protocol §1.5, §2A.4)
death_linked_site = tolower(site_lowercase) %in% c("ohsu", "umn", "umm", "minnesota")
hosp[, death_linked_site := as.integer(death_linked_site)]
message(sprintf("  Death-linked site: %s", ifelse(death_linked_site, "YES", "NO")))

# hospital-level inpatient hospice unit flag (§3.5.2)
hosp[, has_hospice_unit := as.integer(hospital_id %in% hosp_with_hospice_unit)]


## hospice-pathway availability classification (protocol §3.4.4) ---------------
# Same-facility inpatient hospice exists only where a hospital runs a general
# inpatient hospice programme, bills it as a separate encounter, AND maps that
# encounter reliably into CLIF. The has_hospice_unit flag above cannot tell "no
# programme" apart from "programme not captured", so hospital variation in
# same-facility hospice use would conflate structure and data capture with
# clinical practice -- a caterpillar plot would rank billing, not care.
#
# Each hospital is therefore classified as:
#   a_available_captured -- same-facility hospice available and reliably captured
#   b_unavailable        -- no same-facility hospice pathway observed
#   c_uncertain          -- availability or capture uncertain
# The primary same-facility analysis is conditional on class (a). The broader
# network-wide outcome is any documented hospice disposition (§3.4.4).
#
# The protocol sources these classes from the site data-quality returns (§3.7),
# so a SITE ATTESTATION always wins. Absent an attestation file, the rules below
# produce a PROVISIONAL default that the site is expected to confirm.

pathway_min_n     = 10L   # same-facility transitions needed to call capture reliable
pathway_min_years = 0.5   # fraction of observed years that must show >= 1 transition

sf_outcomes = c("same_facility_hospice_death", "same_facility_hospice_alive")

hosp[, admit_year := data.table::year(admission_dttm)]
hosp[, hospital_id := as.character(hospital_id)]

# NB: min()/max() of an empty vector give Inf/-Inf (double), so a hospital with
# no same-facility transitions would return a different column type than one with
# some, and data.table rejects inconsistent types across groups. Branch on length
# so first/last_sf_year are integer for every hospital.
pathway = hosp[in_cohort == 1L, {
  sf     = outcome_category %in% sf_outcomes
  yrs_sf = admit_year[sf]
  .(n_admissions    = .N,
    n_same_facility = sum(sf),
    years_observed  = uniqueN(admit_year),
    years_with_sf   = uniqueN(yrs_sf),
    first_sf_year   = if (length(yrs_sf)) min(yrs_sf) else NA_integer_,
    last_sf_year    = if (length(yrs_sf)) max(yrs_sf) else NA_integer_)
}, by = hospital_id]

pathway[, sf_per_1000 := round(1000 * n_same_facility / n_admissions, 2)]

pathway[, pathway_class := fcase(
  n_same_facility == 0L,
    "b_unavailable",
  n_same_facility >= pathway_min_n &
    years_with_sf >= ceiling(pathway_min_years * years_observed),
    "a_available_captured",
  default = "c_uncertain"
)]

# A pathway that disappears mid-study (a unit that closes) is a structural break,
# not a stable hospital attribute. Flag it and force class (c): a single
# time-invariant class would misdescribe such a hospital, and conditioning the
# primary analysis on it would silently mix pre- and post-closure eras.
max_year = hosp[in_cohort == 1L, max(admit_year, na.rm = TRUE)]
pathway[, pathway_ended_midstudy := as.integer(
  n_same_facility > 0L & !is.na(last_sf_year) & last_sf_year < (max_year - 1L))]
pathway[pathway_ended_midstudy == 1L & pathway_class == "a_available_captured",
        pathway_class := "c_uncertain"]

pathway[, classification_source := "data_driven_provisional"]

# site attestation override (the §3.7 data-quality return)
pathway_map_path = file.path(project_location, "config", "hospital_pathway_map.csv")
if (file.exists(pathway_map_path)) {
  att = fread(pathway_map_path, colClasses = list(character = "hospital_id"))
  if (all(c("hospital_id", "pathway_class") %in% names(att))) {
    att[, hospital_id := as.character(hospital_id)]
    pathway[att, on = "hospital_id", `:=`(
      pathway_class         = i.pathway_class,
      classification_source = "site_attested")]
    message("  Hospice-pathway class: site attestation applied (hospital_pathway_map.csv).")
  } else {
    message("  WARNING: hospital_pathway_map.csv lacks hospital_id/pathway_class; ignoring.")
  }
} else {
  message("  Hospice-pathway class: PROVISIONAL (no hospital_pathway_map.csv supplied).")
  message("    Sites must confirm these against the §3.7 data-quality return.")
}

# carry the class onto encounters (update join: no copy, no reorder)
hosp[pathway, on = "hospital_id", `:=`(
  pathway_class          = i.pathway_class,
  pathway_ended_midstudy = i.pathway_ended_midstudy)]
hosp[is.na(pathway_class), pathway_class := "c_uncertain"]

# network-wide outcome (§3.4.4): any documented hospice disposition, which is
# defined for EVERY hospital regardless of whether it has a same-facility unit.
hosp[, any_hospice_disposition := as.integer(
  outcome_category %in% c(sf_outcomes, "external_hospice_discharge"))]

message("  Hospice-pathway availability (§3.4.4):")
for (i in seq_len(nrow(pathway))) {
  message(sprintf("    hospital %-12s %-22s n_same_facility=%-6s rate/1000=%-6s%s",
                  pathway$hospital_id[i], pathway$pathway_class[i],
                  format(pathway$n_same_facility[i], big.mark = ","),
                  pathway$sf_per_1000[i],
                  ifelse(pathway$pathway_ended_midstudy[i] == 1L,
                         " [pathway ended mid-study]", "")))
}

n_uncertain = pathway[pathway_class == "c_uncertain", .N]
if (n_uncertain > 0) {
  qi_flag(
    "Hospice-pathway availability uncertain",
    sprintf(paste0(
      "%d of %d hospitals fall in class (c): same-facility hospice availability or\n",
      "capture cannot be established from the data (too few transitions, intermittent\n",
      "capture, or a pathway that ended mid-study). Protocol §3.4.4 conditions the\n",
      "PRIMARY same-facility analysis on class (a) only, so these hospitals are\n",
      "excluded from it and appear only in the any-documented-hospice outcome.\n",
      "Confirm each hospital's true status in the §3.7 data-quality return and supply\n",
      "config/hospital_pathway_map.csv to override this provisional classification."),
      n_uncertain, nrow(pathway)))
}

save_output(pathway, "hospital_pathway_availability", shareable = TRUE)

hosp[, admit_year := NULL]

# pandemic period (§3.3)
hosp[, pandemic_period := assign_pandemic_period(admission_dttm)]
message("  Pandemic period distribution:")
pp_tab = hosp[in_cohort == 1L, .N, by = pandemic_period][order(pandemic_period)]
for (i in seq_len(nrow(pp_tab)))
  message(sprintf("    %-18s %s", pp_tab$pandemic_period[i], format(pp_tab$N[i], big.mark = ",")))

cohort = hosp[in_cohort == 1L]

message("\n========== COHORT SUMMARY ==========")
message(sprintf("Total adult hospitalizations in study window:    %s", format(n_total, big.mark = ",")))
message(sprintf("  Excluded — missing/invalid LOS data:           %s", format(n_excl_los, big.mark = ",")))
message(sprintf("    Missing admission or discharge datetime:      %s", format(n_excl_missing_dates, big.mark = ",")))
message(sprintf("    Admission datetime > discharge datetime:      %s", format(n_excl_date_order, big.mark = ",")))
message(sprintf("  Excluded — missing hospital_id:                 %s", format(n_excl_missing_hid, big.mark = ",")))
message(sprintf("  Excluded — admitted from hospice:               %s", format(n_excl_hospice, big.mark = ",")))
message(sprintf("  Excluded — hospital under %s encounters:      %s",
                format(min_hospital_encounters, big.mark = ","),
                format(n_excl_small_hosp, big.mark = ",")))
message(sprintf("  FINAL PRIMARY COHORT:                           %s", format(n_cohort, big.mark = ",")))
message(sprintf("    With any ICU admission:                       %s (%4.1f%%)",
                format(sum(cohort$icu_flag == 1L), big.mark = ","),
                100 * sum(cohort$icu_flag == 1L) / n_cohort))
message(sprintf("    With ED admission:                            %s (%4.1f%%)",
                format(sum(cohort$ed_flag == 1L), big.mark = ","),
                100 * sum(cohort$ed_flag == 1L) / n_cohort))
message("")
message("Outcome distribution (primary cohort):")
outcome_counts = cohort[, .N, by = outcome_category][order(-N)]
for (i in seq_len(nrow(outcome_counts))) {
  message(sprintf("  %-40s %s (%4.1f%%)",
                  outcome_counts$outcome_category[i],
                  format(outcome_counts$N[i], big.mark = ","),
                  100 * outcome_counts$N[i] / n_cohort))
}
message("====================================\n")


# 6. exclusion tracking table (shareable) --------------------------------------

# Structured exclusion table for compilation across all sites.
# Each row is one step in the sequential exclusion cascade.
# "remaining" is the count after that step's exclusion is applied.

exclusion_table = data.table(
  site           = site_lowercase,
  step           = 1:6,
  description    = c(
    "Adult (>=18) inpatient hospitalizations in study window",
    "Excluded: missing or invalid length of stay data",
    "Excluded: missing hospital_id (cannot place at a facility)",
    "Excluded: admitted from hospice",
    sprintf("Excluded: hospital contributed < %s cohort encounters",
            format(min_hospital_encounters, big.mark = ",")),
    "Final primary cohort"
  ),
  n_excluded     = c(0L, n_excl_los, n_excl_missing_hid, n_excl_hospice,
                     n_excl_small_hosp, 0L),
  n_remaining    = c(n_total,
                     n_total - n_excl_los,
                     n_total - n_excl_los - n_excl_missing_hid,
                     n_total - n_excl_los - n_excl_missing_hid - n_excl_hospice,
                     n_cohort, n_cohort)
)

# add detail sub-rows for the LOS exclusion
exclusion_detail = data.table(
  site           = site_lowercase,
  step           = c(2.1, 2.2),
  description    = c(
    "  — Missing admission or discharge datetime",
    "  — Admission datetime after discharge datetime"
  ),
  n_excluded     = c(n_excl_missing_dates, n_excl_date_order),
  n_remaining    = c(NA_integer_, NA_integer_)
)

exclusion_table = rbind(exclusion_table, exclusion_detail)
setorder(exclusion_table, step)

# print to console
message("\n========== EXCLUSION TABLE ==========")
for (i in seq_len(nrow(exclusion_table))) {
  row = exclusion_table[i]
  if (row$n_excluded == 0L && row$step %in% c(1, 5)) {
    message(sprintf("  %-55s n = %s",
                    row$description,
                    format(row$n_remaining, big.mark = ",")))
  } else {
    message(sprintf("  %-55s -%s",
                    row$description,
                    format(row$n_excluded, big.mark = ",")))
  }
}
message("=====================================\n")

# save to final_output for cross-site compilation
fwrite(
  exclusion_table,
  shareable_path(
            paste0("exclusion_table_", site_lowercase, "_", today, ".csv"))
)
message("Exclusion table saved to shareable output.")


# 6b. subgroup summary (shareable) --------------------------------------------

# Outcome distribution within the full cohort and each LST subgroup.
# One row per group × outcome_category for cross-site compilation.

subgroup_defs = list(
  list(label = "full_cohort",   filter_expr = quote(TRUE)),
  list(label = "lst_any",       filter_expr = quote(lst_any_flag == 1L)),
  list(label = "lst_imv",       filter_expr = quote(imv_flag == 1L)),
  list(label = "lst_vaso",      filter_expr = quote(vaso_flag == 1L))
)

subgroup_summary = rbindlist(lapply(subgroup_defs, \(sg) {
  sub = cohort[eval(sg$filter_expr)]
  n_sub = nrow(sub)
  if (n_sub == 0) return(NULL)

  oc = sub[, .N, by = outcome_category]
  oc[, `:=`(
    subgroup    = sg$label,
    site        = site_lowercase,
    n_subgroup  = n_sub,
    pct         = round(100 * N / n_sub, 2)
  )]

  # add summary row for hospice (all types combined)
  n_hospice_all = sum(sub$outcome_category %in% c(
    "same_facility_hospice_death", "same_facility_hospice_alive",
    "external_hospice_discharge"))
  n_short = sum(sub$short_interval_death == 1L, na.rm = TRUE)

  summary_rows = data.table(
    outcome_category = c("hospice_all", "short_survival_hospice"),
    N                = c(n_hospice_all, n_short),
    subgroup         = sg$label,
    site             = site_lowercase,
    n_subgroup       = n_sub,
    pct              = round(100 * c(n_hospice_all, n_short) / n_sub, 2)
  )

  rbind(oc, summary_rows)
}))

setcolorder(subgroup_summary, c("site", "subgroup", "n_subgroup",
                                 "outcome_category", "N", "pct"))

# print to console
message("\n========== SUBGROUP SUMMARY ==========")
for (sg in unique(subgroup_summary$subgroup)) {
  sub = subgroup_summary[subgroup == sg]
  n_sub = sub$n_subgroup[1]
  message(sprintf("\n  %s (n = %s):", sg, format(n_sub, big.mark = ",")))
  for (i in seq_len(nrow(sub))) {
    message(sprintf("    %-40s %s (%5.2f%%)",
                    sub$outcome_category[i],
                    format(sub$N[i], big.mark = ","),
                    sub$pct[i]))
  }
}
message("======================================\n")

# small-cell suppression before cross-site compilation
if (isTRUE(suppress_small_cells_enabled))
  subgroup_summary = suppress_small_cells(subgroup_summary,
                                          count_cols = "N", pct_cols = "pct")

# save for cross-site compilation
fwrite(
  subgroup_summary,
  shareable_path(
            paste0("subgroup_summary_", site_lowercase, "_", today, ".csv"))
)
message("Subgroup summary saved to shareable output.")


# 7. save cohort ---------------------------------------------------------------

# full hospitalization-level file with all flags (intermediate)
save_output(hosp, "cohort_all_hospitalizations", shareable = FALSE)

# primary cohort only
save_output(
  cohort[, .(
    patient_id, hospitalization_id, hospital_id,
    admission_dttm, discharge_dttm, age_at_admission,
    discharge_category, death_dttm,
    ed_flag, icu_flag, icu_48h_flag, age_65_flag,
    imv_flag, vaso_flag, lst_any_flag,
    hospice_location_flag, hospice_transition_dttm,
    outcome_category, hours_transition_to_death, short_interval_death,
    days_discharge_to_death,
    admission_month, discharge_month, los_days,
    pandemic_period, death_linked_site, has_hospice_unit
  )],
  name      = "cohort_primary",
  shareable = FALSE
)

# admitted-from-hospice cohort (for exploratory analysis E.2)
# These patients were excluded from the primary cohort but the protocol
# asks to analyze them as a comparison group.
hospice_admits = hosp[excl_hospice_admit == 1L]
n_hospice_admits = nrow(hospice_admits)
message(sprintf("  Admitted-from-hospice cohort: %s encounters", format(n_hospice_admits, big.mark = ",")))

if (n_hospice_admits > 0) {
  save_output(
    hospice_admits[, .(
      patient_id, hospitalization_id, hospital_id,
      admission_dttm, discharge_dttm, age_at_admission,
      discharge_category, death_dttm,
      ed_flag, icu_flag, icu_48h_flag, age_65_flag,
      imv_flag, vaso_flag, lst_any_flag,
      hospice_location_flag, hospice_transition_dttm,
      outcome_category, hours_transition_to_death, short_interval_death,
      days_discharge_to_death,
      admission_month, discharge_month, los_days,
      pandemic_period, death_linked_site, has_hospice_unit
    )],
    name      = "cohort_admitted_from_hospice",
    shareable = FALSE
  )
}

rm(hospice_admits)

rm(outcome_counts, exclusion_table, exclusion_detail,
   subgroup_summary, subgroup_defs,
   n_total, n_base, n_icu, n_excl_hospice, n_cohort, n_hospice_admits,
   n_excl_los, n_excl_missing_dates, n_excl_date_order,
   n_excl_missing_hid, n_before_hid, n_missing_hid,
   n_imv, n_vaso, n_lst, death_linked_site, hosp_with_hospice_unit)
gc()

# 8. shareable CCN-level summary (for ohsu to compile) -------------------------

# Aggregated to CMS Certification Number (CCN) so that hospitals sharing a CCN
# (e.g., consolidated campuses, "remote location" sub-units) are combined into
# a single row. This matches how CMS-derived datasets identify hospitals.

# Two outputs:
#  hospital_summary_local (intermediate, by hospital_id) — QC, not shared
#  hospital_summary       (shareable, by CCN) → upload_to_box/aim2/

message("Building hospital-level summaries...")

cohort[, hospital_id := as.character(hospital_id)]


## 8a. local hospital_id-level summary (QC, not shared) ----------------------

local_summary = cohort[, .(
  n_admissions               = .N,
  n_routine_discharge        = sum(outcome_category == "routine_discharge", na.rm = TRUE),
  n_in_hospital_death_acute  = sum(outcome_category == "in_hospital_death_acute", na.rm = TRUE),
  n_sf_hospice_death         = sum(outcome_category == "same_facility_hospice_death", na.rm = TRUE),
  n_sf_hospice_alive         = sum(outcome_category == "same_facility_hospice_alive", na.rm = TRUE),
  n_external_hospice         = sum(outcome_category == "external_hospice_discharge", na.rm = TRUE),
  n_transfer_acute           = sum(outcome_category == "transfer_acute_facility", na.rm = TRUE),
  n_other_ama                = sum(outcome_category == "other_ama", na.rm = TRUE),
  n_hospice_all = sum(outcome_category %in% c(
    "same_facility_hospice_death", "same_facility_hospice_alive",
    "external_hospice_discharge"), na.rm = TRUE),
  n_sf_hospice = sum(outcome_category %in% c(
    "same_facility_hospice_death", "same_facility_hospice_alive"), na.rm = TRUE),
  n_short_survival = sum(short_interval_death == 1L, na.rm = TRUE)
), by = hospital_id]

save_output(local_summary, "hospital_summary_by_hospital_id", shareable = FALSE)


## 8b. CCN-level aggregation -------------------------------------------------

if (is.null(hospital_medicare_map)) {
  warning("No CCN map available — skipping shareable CCN-level summary. ",
          "Add hospital_ccn_map.csv to <project_root>/reference/ and re-run.",
          call. = FALSE)
  
} else {
  
  # join CCN onto cohort
  cohort_ccn = merge(
    cohort,
    hospital_medicare_map[, .(hospital_id, medicare_id)],
    by = "hospital_id",
    all.x = TRUE
  )
  
  # check for unmapped hospital_ids
  unmapped_hids = unique(cohort_ccn[is.na(medicare_id), hospital_id])
  n_unmapped_admissions = sum(is.na(cohort_ccn$medicare_id))
  
  if (length(unmapped_hids) > 0) {
    warning(sprintf(
      "%d hospital_id(s) in cohort not found in CCN map — %s admissions excluded from shareable summary:\n  %s\n  Add these to reference/hospital_ccn_map.csv to include them.",
      length(unmapped_hids),
      format(n_unmapped_admissions, big.mark = ","),
      paste(unmapped_hids, collapse = ", ")
    ), call. = FALSE)
  }
  
  # aggregate at CCN level (drop unmapped rows)
  ccn_summary = cohort_ccn[!is.na(medicare_id), .(
    
    # volume
    n_admissions = .N,
    n_hospital_ids_combined = uniqueN(hospital_id),
    
    # outcome counts (multinomial)
    n_routine_discharge        = sum(outcome_category == "routine_discharge", na.rm = TRUE),
    n_in_hospital_death_acute  = sum(outcome_category == "in_hospital_death_acute", na.rm = TRUE),
    n_sf_hospice_death         = sum(outcome_category == "same_facility_hospice_death", na.rm = TRUE),
    n_sf_hospice_alive         = sum(outcome_category == "same_facility_hospice_alive", na.rm = TRUE),
    n_external_hospice         = sum(outcome_category == "external_hospice_discharge", na.rm = TRUE),
    n_transfer_acute           = sum(outcome_category == "transfer_acute_facility", na.rm = TRUE),
    n_other_ama                = sum(outcome_category == "other_ama", na.rm = TRUE),
    
    # hospice aggregates
    n_hospice_all = sum(outcome_category %in% c(
      "same_facility_hospice_death", "same_facility_hospice_alive",
      "external_hospice_discharge"), na.rm = TRUE),
    n_sf_hospice = sum(outcome_category %in% c(
      "same_facility_hospice_death", "same_facility_hospice_alive"), na.rm = TRUE),
    
    n_short_survival = sum(short_interval_death == 1L, na.rm = TRUE),
    
    # covariates
    mean_age     = round(mean(age_at_admission, na.rm = TRUE), 1),
    pct_age_65   = round(100 * mean(age_65_flag == 1L, na.rm = TRUE), 1),
    pct_icu      = round(100 * mean(icu_flag == 1L, na.rm = TRUE), 1),
    pct_icu_48h  = round(100 * mean(icu_48h_flag == 1L, na.rm = TRUE), 1),
    pct_ed_admit = round(100 * mean(ed_flag == 1L, na.rm = TRUE), 1),
    mean_los     = round(mean(los_days, na.rm = TRUE), 2),
    median_los   = round(median(los_days, na.rm = TRUE), 2)
    
  ), by = medicare_id]
  
  # site identifier
  ccn_summary[, site := site_lowercase]
  
  # computed rates (per 1,000 admissions)
  ccn_summary[, `:=`(
    rate_hospice_all = round(1000 * n_hospice_all / n_admissions, 2),
    rate_death_acute = round(1000 * n_in_hospital_death_acute / n_admissions, 2),
    rate_transfer    = round(1000 * n_transfer_acute / n_admissions, 2),
    rate_sf_hospice  = round(1000 * n_sf_hospice / n_admissions, 2),
    rate_ext_hospice = round(1000 * n_external_hospice / n_admissions, 2),
    rate_short_surv  = round(1000 * n_short_survival / n_admissions, 2),
    hospice_mortality_ratio = round(fifelse(
      n_in_hospital_death_acute > 0,
      n_hospice_all / n_in_hospital_death_acute,
      NA_real_), 3)
  )]
  
  # reorder columns: identifiers first
  id_cols = c("site", "medicare_id", "n_hospital_ids_combined")
  other_cols = setdiff(names(ccn_summary), id_cols)
  setcolorder(ccn_summary, c(id_cols, other_cols))
  
  message("CCN-level summary:")
  for (i in seq_len(nrow(ccn_summary))) {
    h = ccn_summary[i]
    message(sprintf("  CCN %s (%d hosp) | N=%s | Deaths=%s | Hospice=%s | Transfer=%s",
                    h$medicare_id, h$n_hospital_ids_combined,
                    format(h$n_admissions, big.mark = ","),
                    format(h$n_in_hospital_death_acute, big.mark = ","),
                    format(h$n_hospice_all, big.mark = ","),
                    format(h$n_transfer_acute, big.mark = ",")))
  }
  
  save_output(ccn_summary, "hospital_summary",
              shareable = TRUE, subfolder = "aim2")
  
  message("  Shareable CCN-level summary saved to shareable output.aim2/")
  
  
  ## 8c. Quarterly CCN-level summary (shareable) -----------------------------
  # Same outcome counts as the overall CCN summary, broken down by discharge
  # quarter for temporal trend analyses across the consortium.
  # Case-mix covariates are omitted (they're stable descriptors of the
  # hospital, not time-varying — they live in the overall summary).
  
  message("Building quarterly CCN-level summary...")
  
  # defensive: ensure discharge_quarter exists (in case cohort was loaded from
  # a parquet file saved before discharge_quarter was added to the save spec)
  if (!"discharge_quarter" %in% names(cohort)) {
    cohort[, discharge_quarter := floor_date(discharge_dttm, "quarter")]
  }
  
  # rebuild the CCN-joined cohort (we rm'd it above; cheap to redo)
  cohort_ccn_q = merge(
    cohort,
    hospital_medicare_map[, .(hospital_id, medicare_id)],
    by = "hospital_id",
    all.x = TRUE
  )
  
  ccn_quarterly = cohort_ccn_q[!is.na(medicare_id) & !is.na(discharge_quarter), .(
    
    n_admissions = .N,
    
    n_routine_discharge        = sum(outcome_category == "routine_discharge", na.rm = TRUE),
    n_in_hospital_death_acute  = sum(outcome_category == "in_hospital_death_acute", na.rm = TRUE),
    n_sf_hospice_death         = sum(outcome_category == "same_facility_hospice_death", na.rm = TRUE),
    n_sf_hospice_alive         = sum(outcome_category == "same_facility_hospice_alive", na.rm = TRUE),
    n_external_hospice         = sum(outcome_category == "external_hospice_discharge", na.rm = TRUE),
    n_transfer_acute           = sum(outcome_category == "transfer_acute_facility", na.rm = TRUE),
    n_other_ama                = sum(outcome_category == "other_ama", na.rm = TRUE),
    
    n_hospice_all = sum(outcome_category %in% c(
      "same_facility_hospice_death", "same_facility_hospice_alive",
      "external_hospice_discharge"), na.rm = TRUE),
    n_sf_hospice = sum(outcome_category %in% c(
      "same_facility_hospice_death", "same_facility_hospice_alive"), na.rm = TRUE),
    
    n_short_survival = sum(short_interval_death == 1L, na.rm = TRUE)
    
  ), by = .(medicare_id, discharge_quarter)]
  
  # site identifier
  ccn_quarterly[, site := site_lowercase]
  
  # quarter label (e.g., "2020 Q3") for human readability
  ccn_quarterly[, quarter_label := paste0(
    year(discharge_quarter), " Q", quarter(discharge_quarter)
  )]
  
  # computed rates
  ccn_quarterly[, `:=`(
    rate_hospice_all = round(1000 * n_hospice_all / n_admissions, 2),
    rate_death_acute = round(1000 * n_in_hospital_death_acute / n_admissions, 2),
    rate_transfer    = round(1000 * n_transfer_acute / n_admissions, 2),
    rate_sf_hospice  = round(1000 * n_sf_hospice / n_admissions, 2),
    rate_ext_hospice = round(1000 * n_external_hospice / n_admissions, 2),
    rate_short_surv  = round(1000 * n_short_survival / n_admissions, 2),
    hospice_mortality_ratio = round(fifelse(
      n_in_hospital_death_acute > 0,
      n_hospice_all / n_in_hospital_death_acute,
      NA_real_), 3)
  )]
  
  # reorder columns
  id_cols_q = c("site", "medicare_id", "discharge_quarter", "quarter_label")
  other_cols_q = setdiff(names(ccn_quarterly), id_cols_q)
  setcolorder(ccn_quarterly, c(id_cols_q, other_cols_q))
  setorder(ccn_quarterly, medicare_id, discharge_quarter)
  
  # report span
  message(sprintf("  Quarterly summary: %d CCN(s) × %d quarter(s) = %d rows",
                  uniqueN(ccn_quarterly$medicare_id),
                  uniqueN(ccn_quarterly$discharge_quarter),
                  nrow(ccn_quarterly)))
  message(sprintf("    Date range: %d Q%d to %d Q%d",
                  year(min(ccn_quarterly$discharge_quarter)),
                  quarter(min(ccn_quarterly$discharge_quarter)),
                  year(max(ccn_quarterly$discharge_quarter)),
                  quarter(max(ccn_quarterly$discharge_quarter))))
  
  save_output(ccn_quarterly, "hospital_summary_quarterly",
              shareable = TRUE, subfolder = "aim2")
  
  message("  Shareable quarterly CCN summary saved to shareable output.aim2/")
  
  rm(cohort_ccn, ccn_summary, cohort_ccn_q, ccn_quarterly)
}


rm(local_summary)
gc()

message("01_cohort.R complete.")