# fn_apache2.R
# Calculate APACHE II scores from CLIF tables.
#
# Components:
#   A. Acute Physiology Score (APS) - 12 variables, worst value in first 48h
#      1.  Temperature (�C)            - vitals: temp_c
#      2.  Mean arterial pressure       - vitals: map
#      3.  Heart rate                   - vitals: heart_rate
#      4.  Respiratory rate             - vitals: respiratory_rate
#      5.  Oxygenation (A-a or PaO2)    - labs: po2_arterial, pco2_arterial + respiratory_support: fio2
#      6.  Arterial pH                  - labs: ph_arterial
#      7.  Sodium                       - labs: sodium
#      8.  Potassium                    - labs: potassium
#      9.  Creatinine                   - labs: creatinine
#      10. Hemoglobin                   - labs: hemoglobin
#      11. WBC                          - labs: wbc
#      12. GCS                          - Not available without patient_assessments;
#                                         defaults to 15 (0 points). Fills in with assessments when available
#   B. Age points
#   C. Chronic health points (optional; derived from diagnosis codes if available)
#
# "Worst" = the value yielding the highest (most abnormal) APACHE II points.
#
# Usage:
#   source("fn_apache2.R")
#   apache_scores = calculate_apache2(data_list, hosp)
#
# Last updated bpe 02.20.2026


# =============================================================================
# APACHE II POINT FUNCTIONS  one per physiological variable
# =============================================================================
# Each function is vectorized: takes a numeric vector, returns integer points.

pts_temp = function(x) {
  fcase(
    x >= 41,                        4L,
    x >= 39   & x <= 40.9,          3L,
    x >= 38.5 & x <= 38.9,          1L,
    x >= 36   & x <= 38.4,          0L,
    x >= 34   & x <= 35.9,          1L,
    x >= 32   & x <= 33.9,          2L,
    x >= 30   & x <= 31.9,          3L,
    x <  30,                         4L,
    default = NA_integer_
  )
}

pts_map = function(x) {
  fcase(
    x >= 160,                        4L,
    x >= 130 & x <= 159,             3L,
    x >= 110 & x <= 129,             2L,
    x >= 70  & x <= 109,             0L,
    x >= 50  & x <= 69,              2L,
    x <  50,                          4L,
    default = NA_integer_
  )
}

pts_hr = function(x) {
  fcase(
    x >= 180,                        4L,
    x >= 140 & x <= 179,             3L,
    x >= 110 & x <= 139,             2L,
    x >= 70  & x <= 109,             0L,
    x >= 55  & x <= 69,              1L,
    x >= 40  & x <= 54,              2L,
    x <  40,                          4L,
    default = NA_integer_
  )
}

pts_rr = function(x) {
  fcase(
    x >= 50,                         4L,
    x >= 35 & x <= 49,              3L,
    x >= 25 & x <= 34,              1L,
    x >= 12 & x <= 24,              0L,
    x >= 10 & x <= 11,              1L,
    x >= 6  & x <= 9,               2L,
    x <  6,                          4L,
    default = NA_integer_
  )
}

pts_ph = function(x) {
  fcase(
    x >= 7.7,                        4L,
    x >= 7.6  & x < 7.7,            3L,
    x >= 7.5  & x < 7.6,            1L,
    x >= 7.33 & x < 7.5,            0L,
    x >= 7.25 & x < 7.33,           2L,
    x >= 7.15 & x < 7.25,           3L,
    x <  7.15,                        4L,
    default = NA_integer_
  )
}

pts_sodium = function(x) {
  fcase(
    x >= 180,                        4L,
    x >= 160 & x <= 179,            3L,
    x >= 155 & x <= 159,            2L,
    x >= 150 & x <= 154,            1L,
    x >= 130 & x <= 149,            0L,
    x >= 120 & x <= 129,            2L,
    x >= 111 & x <= 119,            3L,
    x <= 110,                        4L,
    default = NA_integer_
  )
}

pts_potassium = function(x) {
  fcase(
    x >= 7,                          4L,
    x >= 6   & x < 7,               3L,
    x >= 5.5 & x < 6,               1L,
    x >= 3.5 & x < 5.5,             0L,
    x >= 3   & x < 3.5,             1L,
    x >= 2.5 & x < 3,               2L,
    x <  2.5,                         4L,
    default = NA_integer_
  )
}

pts_creatinine = function(x, acute_renal_failure = FALSE) {
  base = fcase(
    x >= 3.5,                        4L,
    x >= 2   & x < 3.5,             3L,
    x >= 1.5 & x < 2,               2L,
    x >= 0.6 & x < 1.5,             0L,
    x <  0.6,                       2L,
    default = NA_integer_
  )
  # double points if acute renal failure
  if (isTRUE(acute_renal_failure)) base * 2L else base
}

pts_hemoglobin = function(x) {
  fcase(
    x >= 180,                        4L,
    x >= 150 & x < 180,              2L,
    x >= 138 & x < 150,              1L,
    x >= 90 & x < 138,               0L,
    x >= 60 & x < 90,                2L,
    x <  60,                         4L,
    default = NA_integer_
  )
}

pts_wbc = function(x) {
  fcase(
    x >= 40,                         4L,
    x >= 20 & x < 40,               2L,
    x >= 15 & x < 20,               1L,
    x >= 3  & x < 15,               0L,
    x >= 1  & x < 3,                2L,
    x <  1,                           4L,
    default = NA_integer_
  )
}

# oxygenation: depends on FiO2
# if FiO2 >= 0.5 - use A-a gradient = (FiO2 � 713) - PaCO2 - PaO2
# if FiO2 <  0.5 - use PaO2 directly
pts_oxygenation = function(po2_arterial, pco2_arterial, fio2) {
  aa_gradient = (fio2 * 713) - pco2_arterial - po2_arterial
  
  fcase(
    # high FiO2 path: use A-a gradient
    fio2 >= 0.5 & aa_gradient >= 500,                 4L,
    fio2 >= 0.5 & aa_gradient >= 350 & aa_gradient < 500, 3L,
    fio2 >= 0.5 & aa_gradient >= 200 & aa_gradient < 350, 2L,
    fio2 >= 0.5 & aa_gradient <  200,                 0L,
    
    # low FiO2 path: use PaO2
    fio2 < 0.5 & po2_arterial < 55,                       4L,
    fio2 < 0.5 & po2_arterial >= 55 & po2_arterial <= 60, 3L,
    fio2 < 0.5 & po2_arterial >  60 & po2_arterial <= 70, 1L,
    fio2 < 0.5 & po2_arterial >  70,                      0L,
    
    default = NA_integer_
  )
}

pts_gcs = function(gcs) {
  # APACHE II GCS component = 15 - GCS
  as.integer(15L - pmin(pmax(gcs, 3L), 15L))
}

pts_age = function(age) {
  fcase(
    age <= 44,              0L,
    age >= 45 & age <= 54,  2L,
    age >= 55 & age <= 64,  3L,
    age >= 65 & age <= 74,  5L,
    age >= 75,              6L,
    default = NA_integer_
  )
}


# =============================================================================
# MAIN CALCULATION FUNCTION
# =============================================================================

calculate_apache2 = function(data_list, hosp, window_hours = 48) {
  
  message("Calculating APACHE II scores (first ", window_hours, "h window)...")
  
  # reference table: hospitalization_id � admission_dttm, age
  ref = hosp[, .(hospitalization_id, admission_dttm, age_at_admission)]
  
  
  # ---------------------------------------------------------------------------
  # A. VITALS - temp_c, map, heart_rate, respiratory_rate
  # ---------------------------------------------------------------------------
  
  message("  Pulling vitals...")
  
  vital_cats = c("temp_c", "map", "heart_rate", "respiratory_rate")
  
  vitals = data_list[["vitals"]] |>
    dplyr::filter(vital_category %in% vital_cats) |>
    dplyr::select(hospitalization_id, vital_category, vital_value, recorded_dttm) |>
    dplyr::collect() |>
    as.data.table()
  
  vitals[, vital_num := as.numeric(vital_value)]
  vitals = vitals[!is.na(vital_num)]
  
  # restrict to first 48h
  vitals = merge(vitals, ref[, .(hospitalization_id, admission_dttm)], by = "hospitalization_id")
  vitals[, hours := as.numeric(difftime(recorded_dttm, admission_dttm, units = "hours"))]
  vitals = vitals[hours >= 0 & hours <= window_hours]
  
  # compute points per measurement
  vitals[, pts := fcase(
    vital_category == "temp_c",           pts_temp(vital_num),
    vital_category == "map",              pts_map(vital_num),
    vital_category == "heart_rate",       pts_hr(vital_num),
    vital_category == "respiratory_rate", pts_rr(vital_num),
    default = NA_integer_
  )]
  
  # worst (max points) per variable per encounter
  vital_worst = vitals[
    !is.na(pts),
    .(pts = max(pts)),
    by = .(hospitalization_id, vital_category)
  ]
  
  vital_wide = dcast(
    vital_worst,
    hospitalization_id ~ vital_category,
    value.var = "pts",
    fill = NA_integer_
  )
  
  setnames(vital_wide,
           old = c("temp_c", "map", "heart_rate", "respiratory_rate"),
           new = c("pts_temp", "pts_map", "pts_hr", "pts_rr"),
           skip_absent = TRUE
  )
  
  rm(vitals, vital_worst)
  
  
  # ---------------------------------------------------------------------------
  # B. LABS - pH, sodium, potassium, creatinine, hemoglobin, WBC, PaO2, PaCO2
  # ---------------------------------------------------------------------------
  
  message("  Pulling labs...")
  
  lab_cats = c(
    "ph_arterial", "sodium", "potassium", "creatinine",
    "hemoglobin", "wbc", "po2_arterial", "pco2_arterial"
  )
  
  labs = data_list[["labs"]] |>
    dplyr::filter(lab_category %in% lab_cats) |>
    dplyr::select(hospitalization_id, lab_category, lab_value, lab_result_dttm) |>
    dplyr::collect() |>
    as.data.table()
  
  labs[, lab_num := as.numeric(lab_value)]
  labs = labs[!is.na(lab_num)]
  
  labs = merge(labs, ref[, .(hospitalization_id, admission_dttm)], by = "hospitalization_id")
  labs[, hours := as.numeric(difftime(lab_result_dttm, admission_dttm, units = "hours"))]
  labs = labs[hours >= 0 & hours <= window_hours]
  
  # --- score non-oxygenation labs (worst per variable per encounter) ---
  
  score_labs = labs[!lab_category %in% c("po2_arterial", "pco2_arterial")]
  
  score_labs[, pts := fcase(
    lab_category == "ph_arterial", pts_ph(lab_num),
    lab_category == "sodium",      pts_sodium(lab_num),
    lab_category == "potassium",   pts_potassium(lab_num),
    lab_category == "creatinine",  pts_creatinine(lab_num),
    lab_category == "hemoglobin",  pts_hemoglobin(lab_num),
    lab_category == "wbc",         pts_wbc(lab_num),
    default = NA_integer_
  )]
  
  lab_worst = score_labs[
    !is.na(pts),
    .(pts = max(pts)),
    by = .(hospitalization_id, lab_category)
  ]
  
  lab_wide = dcast(
    lab_worst,
    hospitalization_id ~ lab_category,
    value.var = "pts",
    fill = NA_integer_
  )
  
  # rename to pts_ prefix
  lab_cols = intersect(lab_cats[!lab_cats %in% c("po2_arterial", "pco2_arterial")], names(lab_wide))
  setnames(lab_wide, old = lab_cols, new = paste0("pts_", lab_cols), skip_absent = TRUE)
  
  rm(score_labs, lab_worst)
  
  
  # ---------------------------------------------------------------------------
  # C. OXYGENATION - needs PaO2 + PaCO2 from labs, FiO2 from respiratory_support
  # ---------------------------------------------------------------------------
  
  message("  Calculating oxygenation component...")
  
  # pull FiO2 from respiratory_support
  resp = data_list[["respiratory_support"]] |>
    dplyr::select(hospitalization_id, recorded_dttm, any_of(c("fio2_set", "fio2"))) |>
    dplyr::collect() |>
    as.data.table()
  
  # harmonize FiO2 column name
  if ("fio2_set" %in% names(resp) && !"fio2" %in% names(resp)) {
    setnames(resp, "fio2_set", "fio2")
  }
  
  if ("fio2" %in% names(resp)) {
    resp[, fio2 := as.numeric(fio2)]
    # convert percentages to fractions if needed
    resp[fio2 > 1 & !is.na(fio2), fio2 := fio2 / 100]
    resp = resp[!is.na(fio2) & fio2 > 0]
    
    resp = merge(resp, ref[, .(hospitalization_id, admission_dttm)], by = "hospitalization_id")
    resp[, hours := as.numeric(difftime(recorded_dttm, admission_dttm, units = "hours"))]
    resp = resp[hours >= 0 & hours <= window_hours]
    setnames(resp, "recorded_dttm", "resp_dttm")
    
    # ABG values (PaO2 and PaCO2) within window
    abg = labs[lab_category %in% c("po2_arterial", "pco2_arterial")]
    abg_wide = dcast(
      abg,
      hospitalization_id + lab_result_dttm ~ lab_category,
      value.var = "lab_num",
      fun.aggregate = function(x) x[1]  # take first if duplicates at same timestamp
    )
    
    # require both PaO2 and PaCO2 at the same draw
    abg_wide = abg_wide[!is.na(po2_arterial) & !is.na(pco2_arterial)]
    
    if (nrow(abg_wide) > 0 && nrow(resp) > 0) {
      # rolling join: for each ABG, find nearest prior FiO2 (within 4h)
      setkey(resp, hospitalization_id, resp_dttm)
      setkey(abg_wide, hospitalization_id, lab_result_dttm)
      
      oxy = resp[abg_wide, roll = 4 * 3600, rollends = c(TRUE, TRUE)]
      
      oxy[, pts_oxy := pts_oxygenation(po2_arterial, pco2_arterial, fio2)]
      
      oxy_worst = oxy[
        !is.na(pts_oxy),
        .(pts_oxygenation = max(pts_oxy)),
        by = hospitalization_id
      ]
    } else {
      oxy_worst = data.table(
        hospitalization_id = character(),
        pts_oxygenation = integer()
      )
    }
  } else {
    message("    NOTE: FiO2 not found in respiratory_support; oxygenation scored as NA.")
    oxy_worst = data.table(
      hospitalization_id = character(),
      pts_oxygenation = integer()
    )
  }
  
  rm(resp, labs)
  if (exists("abg_wide")) rm(abg_wide)
  if (exists("oxy")) rm(oxy)
  
  
  # ---------------------------------------------------------------------------
  # D. GCS - from patient_assessments (assessment_category = "gcs_total")
  # ---------------------------------------------------------------------------
  
  message("  Pulling GCS from patient_assessments...")
  
  if ("patient_assessments" %in% names(data_list)) {
    
    gcs = data_list[["patient_assessments"]] |>
      dplyr::filter(tolower(assessment_category) == "gcs_total") |>
      dplyr::select(hospitalization_id, recorded_dttm, numerical_value) |>
      dplyr::collect() |>
      as.data.table()
    
    gcs[, gcs_val := as.numeric(numerical_value)]
    gcs = gcs[!is.na(gcs_val)]
    
    gcs = merge(gcs, ref[, .(hospitalization_id, admission_dttm)], by = "hospitalization_id")
    gcs[, hours := as.numeric(difftime(recorded_dttm, admission_dttm, units = "hours"))]
    gcs = gcs[hours >= 0 & hours <= window_hours]
    
    # worst GCS = lowest value (yields highest APACHE points)
    gcs_worst = gcs[, .(gcs_min = min(gcs_val)), by = hospitalization_id]
    gcs_worst[, pts_gcs_val := pts_gcs(as.integer(gcs_min))]
    
    message(sprintf("    GCS available for %s encounters.", format(nrow(gcs_worst), big.mark = ",")))
    rm(gcs)
    
  } else {
    message("    WARNING: patient_assessments not in data_list; GCS defaulting to 15 (0 pts).")
    gcs_worst = data.table(
      hospitalization_id = character(),
      gcs_min = numeric(),
      pts_gcs_val = integer()
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # E. AGE POINTS
  # ---------------------------------------------------------------------------
  
  age_pts = ref[, .(hospitalization_id, pts_age = pts_age(age_at_admission))]
  
  
  # ---------------------------------------------------------------------------
  # F. CHRONIC HEALTH POINTS (simplified)
  # ---------------------------------------------------------------------------
  # Full APACHE II awards 5 pts (emergency/non-op) or 2 pts (elective postop)
  # for severe chronic organ insufficiency or immunocompromised state.
  # This requires chart review; a rough approximation uses ICD codes.
  # Set to 0 by default; uncomment and adapt the block below if needed.
  
  chronic_pts = ref[, .(hospitalization_id, pts_chronic = 0L)]
  
  # --- OPTIONAL: ICD-based approximation ---
  # Uncomment to flag chronic organ insufficiency / immunocompromised from
  # hospital_diagnosis. These codes are intentionally conservative.
  #
  # chronic_icd10 = c(
  #   "K72",  "K70.4",           # hepatic failure / cirrhosis
  #   "I50.2", "I50.3", "I50.4", # NYHA class III/IV heart failure
  #   "J96.1",                    # chronic respiratory failure
  #   "N18.5", "N18.6",          # CKD stage 5 / ESRD
  #   "D80", "D81", "D82", "D83", "D84",  # immunodeficiency
  #   "B20"                       # HIV
  
  # )
  # chronic_dx = data_list[["hospital_diagnosis"]] |>
  #   dplyr::select(hospitalization_id, diagnosis_code) |>
  #   dplyr::collect() |>
  #   as.data.table()
  # chronic_dx[, code_clean := gsub("\\.", "", toupper(trimws(diagnosis_code)))]
  # chronic_flag = chronic_dx[
  #   grepl(paste0("^(", paste(gsub("\\.", "", chronic_icd10), collapse = "|"), ")"), code_clean),
  #   .(pts_chronic = 5L),     # conservative: assume emergency/non-op
  #   by = hospitalization_id
  # ]
  # chronic_pts = merge(
  #   ref[, .(hospitalization_id)], chronic_flag,
  #   by = "hospitalization_id", all.x = TRUE
  # )
  # chronic_pts[is.na(pts_chronic), pts_chronic := 0L]
  
  
  # ---------------------------------------------------------------------------
  # G. ASSEMBLE APACHE II SCORE
  # ---------------------------------------------------------------------------
  
  message("  Assembling scores...")
  
  apache = ref[, .(hospitalization_id)]
  
  # merge all component tables
  apache = merge(apache, vital_wide,  by = "hospitalization_id", all.x = TRUE)
  apache = merge(apache, lab_wide,    by = "hospitalization_id", all.x = TRUE)
  apache = merge(apache, oxy_worst,   by = "hospitalization_id", all.x = TRUE)
  apache = merge(apache, gcs_worst[, .(hospitalization_id, pts_gcs_val)],
                 by = "hospitalization_id", all.x = TRUE)
  apache = merge(apache, age_pts,     by = "hospitalization_id", all.x = TRUE)
  apache = merge(apache, chronic_pts, by = "hospitalization_id", all.x = TRUE)
  
  # GCS: use actual score where available, default to 0 (GCS 15) where missing
  apache[, pts_gcs := fifelse(!is.na(pts_gcs_val), pts_gcs_val, 0L)]
  apache[, pts_gcs_val := NULL]
  
  # identify all point columns
  pts_cols = grep("^pts_", names(apache), value = TRUE)
  
  # count how many components are available per encounter
  apache[, n_components := rowSums(!is.na(.SD)), .SDcols = pts_cols]
  
  # sum points (NA components contribute 0  conservative)
  apache[, apache2_score := rowSums(.SD, na.rm = TRUE), .SDcols = pts_cols]
  
  # flag completeness: 12 APS variables + age + chronic = 14 components
  apache[, apache2_complete := fifelse(n_components >= 10L, TRUE, FALSE)]
  
  message(sprintf(
    "  APACHE II calculated for %s encounters (%s with e10 components).",
    format(nrow(apache), big.mark = ","),
    format(sum(apache$apache2_complete, na.rm = TRUE), big.mark = ",")
  ))
  
  rm(vital_wide, lab_wide, oxy_worst, gcs_worst, age_pts, chronic_pts)
  gc()
  
  return(apache)
}
