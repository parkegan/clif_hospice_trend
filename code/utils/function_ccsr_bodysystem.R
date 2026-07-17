# function_ccsr_bodysystem.R
# Assign AHRQ CCSR principal-diagnosis BODY SYSTEM to ICD-10-CM codes.
#
# Method: use the *inpatient default* CCSR category — the single category AHRQ
# assigns when a code is the principal/first-listed inpatient diagnosis — and
# roll it up to one of 22 body systems via the 3-letter CCSR category prefix.
# This is the correct CCSR construct for classifying an admission diagnosis.
#
# Lookup file (built once from DXCCSR-Reference-File-v2026-1.xlsx), place at:
#   reference/ccsr_principal_lookup.csv
# columns: icd10cm_code, ccsr_category, ccsr_category_desc,
#          body_system, body_system_short, body_system_full
#
# NOTES
#  - "XXX" = ICD/CMS "unacceptable as principal diagnosis" (manifestation /
#    sequela / external-cause codes). When one appears as a principal dx it
#    signals a coding/selection quirk — kept as its own label, not silently
#    dropped, so you can see it.
#  - Codes absent from the lookup -> "Unmapped". A multi-year cohort WILL have
#    some (codes retired before FY2026, or invalid). Inspect the unmapped set;
#    if the rate is non-trivial, layer an older CCSR version for those years.
#
# Last updated bpe

library(data.table)

# ---- load the lookup --------------------------------------------------------
load_ccsr_lookup = function(path = NULL) {
  if (is.null(path)) {
    cand = c(
      tryCatch(here::here("reference", "ccsr_principal_lookup.csv"),
               error = function(e) NA_character_),
      if (exists("project_root"))
        file.path(project_root, "reference", "ccsr_principal_lookup.csv")
    )
    cand = cand[!is.na(cand) & file.exists(cand)]
    path = if (length(cand)) cand[1] else NA_character_
  }
  if (is.na(path) || !file.exists(path))
    stop("CCSR lookup not found. Place ccsr_principal_lookup.csv ",
         "in reference/ (or pass an explicit path).")
  lk = fread(path, colClasses = list(character = "icd10cm_code"))
  setkey(lk, icd10cm_code)
  lk[]
}

# ---- select the principal diagnosis per encounter ---------------------------
# Prefer an explicit primary/priority flag over data order. `dx_icd10` must have
# hospitalization_id + code_clean; pass the detected flag column name(s).
select_principal_dx = function(dx_icd10, flag_col = NULL) {
  d = dx_icd10
  if (!is.null(flag_col) && flag_col %in% names(d)) {
    v = d[[flag_col]]
    if (grepl("rank|priority|sequence|seq|position", flag_col, ignore.case = TRUE)) {
      d = copy(d); d[, .prio := suppressWarnings(as.numeric(get(flag_col)))]
      flagged = d[!is.na(.prio)]
      if (nrow(flagged)) {
        idx = flagged[, .I[which.min(.prio)], by = hospitalization_id]$V1
        return(flagged[idx, .(hospitalization_id, principal_code = code_clean)])
      }
    } else {
      is_primary = tolower(trimws(as.character(v))) %in%
        c("1","true","t","y","yes","primary","principal")
      flagged = d[is_primary]
      if (nrow(flagged)) {
        idx = flagged[, .I[1L], by = hospitalization_id]$V1
        return(flagged[idx, .(hospitalization_id, principal_code = code_clean)])
      }
    }
  }
  # fallback: first-listed code in data order (last resort)
  idx = d[, .I[1L], by = hospitalization_id]$V1
  d[idx, .(hospitalization_id, principal_code = code_clean)]
}

# ---- assign CCSR body system ------------------------------------------------
# codes: character vector of principal ICD-10 codes (any format; cleaned here).
# returns a data.table aligned row-for-row with `codes`.
assign_ccsr_bodysystem = function(codes, lookup = load_ccsr_lookup()) {
  dt = data.table(icd10cm_code = toupper(gsub("\\.", "", trimws(codes))))
  out = lookup[dt, on = "icd10cm_code"]   # left join, preserves order/length
  out[is.na(body_system), `:=`(
    ccsr_category      = "UNMAPPED",
    ccsr_category_desc = "Code not in CCSR (retired/invalid)",
    body_system        = "UNM",
    body_system_short  = "Unmapped",
    body_system_full   = "Code not found in CCSR lookup"
  )]
  out[]
}
