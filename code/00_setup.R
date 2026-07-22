# 00_setup.R
# This script loads packages, detects system resources, reads config,
# loads and validates CLIF tables, and defines shared utility functions.
# Last updated bpe 22.07.2026


# install and load packages ----------------------------------------------------

packages_to_install = c(
  "lubridate", 
  "rprojroot", 
  "tidytable", 
  "collapse",
  "janitor", 
  "arrow", 
  "readr", 
  "yaml", 
  "here"
)

packages_to_load = c(
  "rprojroot", 
  "data.table", 
  "lubridate", 
  "tidytable", 
  "collapse",
  "stringr", 
  "arrow", 
  "here", 
  "dplyr", 
  "rlang", 
  "tidyr"
)

invisible(lapply(packages_to_install, \(p) {
  if (!requireNamespace(p, quietly = TRUE))
    try(install.packages(p, dependencies = TRUE), silent = TRUE)
}))

invisible(lapply(packages_to_load, \(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

options(dplyr.summarise.inform = FALSE)

rm(packages_to_install, packages_to_load)
gc()


# threads and ram --------------------------------------------------------------

os_type   = Sys.info()[["sysname"]]
all_cores = parallel::detectCores(logical = TRUE)
all_cores = if (is.na(all_cores)) 1L else as.integer(all_cores)

get_ram_gb = function() {
  tryCatch({
    if (os_type == "Darwin") {
      bytes = as.numeric(suppressWarnings(
        system("sysctl -n hw.memsize", intern = TRUE)
      ))
      if (length(bytes) && !is.na(bytes)) bytes / 1024^3 else NA_real_
      
    } else if (file.exists("/proc/meminfo")) {
      kb = as.numeric(suppressWarnings(
        system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)
      ))
      if (length(kb) && !is.na(kb)) kb / 1024^2 else NA_real_
      
    } else {
      mem = ps::ps_system_memory()
      key = intersect(c("available", "avail"), names(mem))
      if (length(key)) mem[[key[1]]] / 1024^3 else NA_real_
    }
  }, error = \(e) NA_real_)
}

avail_ram_gb = get_ram_gb()

## choose threads conservatively -----------------------------------------------

reserve_cores  = 1L
gb_per_thread  = 0.50
max_by_cores   = max(1L, all_cores - reserve_cores)
max_by_memory  = if (is.finite(avail_ram_gb)) max(1L, floor(avail_ram_gb / gb_per_thread)) else max_by_cores
n_threads      = as.integer(min(max_by_cores, max_by_memory, 8L))
n_threads      = max(1L, n_threads)


## apply thread settings -------------------------------------------------------

data.table::setDTthreads(threads = n_threads)
collapse::set_collapse(nthreads  = n_threads)
options(arrow.use_threads        = TRUE)
Sys.setenv(ARROW_NUM_THREADS     = n_threads)
options(mc.cores                 = n_threads)

message(sprintf(
  "Env OK | OS=%s | Cores=%d | Threads=%d | Avail RAM\u2248%s GB",
  os_type, all_cores, n_threads,
  if (is.finite(avail_ram_gb)) round(avail_ram_gb, 1) else "NA"))


# site details -----------------------------------------------------------------

project_root = find_root(rprojroot::has_dir("config"))
config       = yaml::read_yaml(file.path(project_root, "config", "config.yaml"))

allowed_files = c("parquet", "csv", "fst")


## user enters site details ----------------------------------------------------

site_raw       = config$site_lowercase 
site_lowercase = if (is.null(site_raw) || site_raw == "") "site" else tolower(site_raw)

tz_local         = config$time_zone 
file_type        = tolower(config$file_type)
tables_location  = config$tables_location

project_location = config$project_location
if (is.null(project_location) || project_location == "") {
  project_location = project_root
}

message("Site: ", site_lowercase)
message("Project location: ", project_location)


## hospital → CCN map ----------------------------------------------------------

# Loads a multi-site CSV that maps each CLIF hospital_id to its CMS Certification
# Number (CCN). Filters to the current site. Hospitals sharing a CCN
# (e.g., consolidated campuses, remote locations) will be aggregated together
# in the cohort summary to match how CMS/Vizient identify hospitals.

# Default location: <project_root>/reference/hospital_ccn_map.csv
# Override in config.yaml with: hospital_ccn_map_file: "/path/to/file.csv"

# CSV columns required:
#   site             - matches site_lowercase / institution (case-insensitive)
#   hospital_id      - matches hospital_id values in CLIF adt table
#   hospital_ccn_id  - 6-digit CMS Certification Number (zero-padded)

# Network standard population for the Aim 2A.3 risk-standardized composition
# (federation). A version-tagged, NON-identifying table of RAW covariate profiles
# (age, sex_category, race_eth, dx_category, elixhauser_vw, apache2, pandemic_period,
# and where used preexisting_code_status / admission_source) that EVERY site
# standardizes onto, so hospital estimates are comparable across sites. Distributed
# by the coordinating center. Until the group finalizes it, the file is absent and
# 04 falls back to site-local standardization (single-site behavior is unchanged).
# CCSR -> acuity diagnosis-category crosswalk (Aim 2A case-mix adjustment). Maps
# the 530-level CCSR category to an 11-level acuity grouping tailored to end-of-
# life disposition (cancer, dementia, frailty, chronic organ failure, and so on),
# which is more prognostically relevant than body system. If the file is absent,
# 04 falls back to the body-system grouping from function_ccsr_bodysystem.R.
ccsr_acuity_path = config$ccsr_acuity_crosswalk_file %||%
  file.path(project_root, "reference", "crosswalk_ccsr_to_acuity.csv")

load_ccsr_acuity_crosswalk = function() {
  if (!file.exists(ccsr_acuity_path)) {
    message("CCSR->acuity crosswalk: not supplied -> dx_category stays body-system.")
    return(NULL)
  }
  cw = tryCatch(fread(ccsr_acuity_path), error = function(e) {
    warning("Could not read CCSR->acuity crosswalk: ", conditionMessage(e), call. = FALSE); NULL })
  if (!is.null(cw)) {
    n_flag = if ("review_flag" %in% names(cw)) sum(cw$review_flag == 1L, na.rm = TRUE) else NA
    message(sprintf("CCSR->acuity crosswalk: %s rows, %s acuity categories%s.",
                    format(nrow(cw), big.mark = ","),
                    if ("acuity_category" %in% names(cw)) length(unique(cw$acuity_category)) else "?",
                    if (!is.na(n_flag) && n_flag > 0) sprintf(", %d still flagged for review", n_flag) else ""))
  }
  cw
}

# ---- shared demographic recodes & small-cell handling (federation) -----------
# Demographics MUST be categorized identically at every site for cross-site
# compilation, so the recodes live here and are called from the reporting scripts.
# Two small-cell strategies from the analysis discussion are BOTH provided, to be
# assessed once we see how many sites hit sparse cells:
#   (A) a collapsed race/ethnicity "mashup" (fewer, larger cells), and
#   (B) suppression of counts below a threshold (for coordinating-center noise
#       injection with denominator compensation).

# Small-cell disclosure control. These are deliberately HARD-CODED constants and
# are NOT read from config.yaml. Disclosure protection has to be identical at
# every site for the compiled tables to be comparable, and config.yaml is
# site-specific and deliberately not version controlled, so a site could weaken
# the threshold or switch suppression off entirely with no record of the change.
# Adjusting these requires a code edit, which is reviewed and tracked in git for
# the whole consortium. A site whose local IRB requires a HIGHER threshold should
# raise the value here rather than add a config key.
small_cell_threshold         = 5L
suppress_small_cells_enabled = TRUE

# Age bands. The former 18-39 group was too sparse in the hospice cells, so the
# non-elderly are one 18-64 band; 65-79 and 80+ carry the end-of-life signal.
recode_age_group = function(dt) {
  dt[, age_group := fcase(
    age_at_admission < 65L, "18-64",
    age_at_admission < 80L, "65-79",
    default = "80+")]
  dt
}

# Sex as a binary: Female vs not-female. Unknown and missing are grouped with the
# non-female side because Unknown is itself a very small cell.
recode_sex_female = function(dt) {
  dt[, sex_female := fifelse(tolower(trimws(sex_category)) == "female",
                             "Female", "Not female")]
  dt
}

# Race for the granular Table-1 block, with the small Unknown and Other cells
# merged into one "Other/Unknown" level.
recode_race_report = function(dt) {
  dt[, race_report := fifelse(
    tolower(trimws(race_category)) %in% c("unknown", "other", "") |
      is.na(race_category),
    "Other/Unknown", as.character(race_category))]
  dt
}

# Combined race/ethnicity (the MODEL covariate in 04, and reported in Table 1):
# one variable merging race and ethnicity, with the smallest cells (Other +
# Unknown) collapsed together.
recode_race_eth = function(dt) {
  dt[, race_eth := fcase(
    tolower(ethnicity_category) %in% c("hispanic", "hispanic or latino"), "Hispanic",
    tolower(race_category) %in% c("white", "caucasian"), "White non-Hispanic",
    tolower(race_category) %in% c("black", "african american",
                                  "black or african american"), "Black non-Hispanic",
    tolower(race_category) %in% c("asian"), "Asian",
    default = "Other/Unknown")]
  dt
}

# (A) Binary non-Hispanic-white indicator ("nhw01"): the maximally collapsed
# race/ethnicity option, for assessing whether the finer scheme is even needed.
# 1 = non-Hispanic white; 0 = everyone else (Hispanic of any race, or non-white).
# Emitted both as a 0/1 integer (nhw01) and a labeled factor (nhw_group) for
# readable Table-1 rows.
recode_nonhispanic_white = function(dt) {
  dt[, nhw01 := as.integer(
    tolower(race_category) %in% c("white", "caucasian") &
      !tolower(ethnicity_category) %in% c("hispanic", "hispanic or latino"))]
  dt[, nhw_group := fifelse(nhw01 == 1L, "Non-Hispanic White", "Other")]
  dt
}

# (B) small-cell suppression for shareable COUNT tables. A non-zero count below
# the threshold is set to NA and flagged in a `suppressed` column; group totals
# and denominators are left intact so the coordinating center can inject
# 0..(threshold-1) with denominator compensation. True zeros are NOT suppressed
# (a zero is not disclosive). Continuous-summary rows (whose n is the group total,
# far above the threshold) are untouched.
suppress_small_cells = function(dt, count_cols = "n", pct_cols = "pct",
                                threshold = small_cell_threshold) {
  out = data.table::copy(as.data.table(dt))
  # guard: a caller may RAISE the threshold but never lower it below the study
  # minimum, so no call site can weaken disclosure control locally
  threshold = max(as.integer(threshold), small_cell_threshold, na.rm = TRUE)
  count_cols = intersect(count_cols, names(out))
  if (!length(count_cols)) return(out)
  is_small = function(cc) !is.na(out[[cc]]) & out[[cc]] > 0 & out[[cc]] < threshold
  # evaluate every column BEFORE any masking, so the flags are not affected by the
  # type conversion performed in the loop below
  small_list = lapply(count_cols, is_small)
  # row-level flag: any count cell in the row is a small non-zero cell
  out[, suppressed := as.integer(Reduce(`|`, small_list))]
  # Mask each small count cell with a "<threshold" LABEL rather than a blank, so a
  # reader can tell the value was computed and is simply small, not missing or
  # failed. This converts the count column to character; that is intentional, and
  # happens immediately before the table is written, so no downstream arithmetic
  # is affected. True zeros and counts at or above the threshold keep their value.
  lab = paste0("<", threshold)
  for (i in seq_along(count_cols)) {
    cc  = count_cols[[i]]
    val = out[[cc]]
    chr = rep(NA_character_, length(val))
    keep = !is.na(val)
    chr[keep] = format(val[keep], scientific = FALSE, trim = TRUE)
    chr[small_list[[i]]] = lab
    data.table::set(out, j = cc, value = chr)
  }
  # A percent derived from a suppressed count is itself disclosive (it can be
  # inverted with the denominator), so blank it rather than labelling it.
  for (pc in intersect(pct_cols, names(out))) out[suppressed == 1L, (pc) := NA_real_]
  out
}


standard_pop_path = config$network_standard_population_file %||%
  file.path(project_root, "config", "network_standard_population.csv")
standard_pop_version = config$network_standard_population_version %||% "site_local"

load_standard_population = function() {
  if (!file.exists(standard_pop_path)) {
    message("Network standard population: not supplied -> 04 will standardize site-locally.")
    return(NULL)
  }
  sp = tryCatch(fread(standard_pop_path), error = function(e) {
    warning("Could not read network standard population: ", conditionMessage(e), call. = FALSE)
    NULL
  })
  if (!is.null(sp))
    message(sprintf("Network standard population: %s rows, version '%s'.",
                    format(nrow(sp), big.mark = ","), standard_pop_version))
  sp
}

ccn_map_path = config$hospital_ccn_map_file %||%
  file.path(project_root, "reference", "hospital_ccn_map.csv")

hospital_medicare_map = NULL

if (!file.exists(ccn_map_path)) {
  warning("Hospital CCN map not found at: ", ccn_map_path,
          "\n  CCN-level aggregation will be skipped.",
          "\n  Place hospital_ccn_map.csv in <project_root>/reference/ ",
          "or set 'hospital_ccn_map_file' in config.yaml.", call. = FALSE)
  
} else {
  
  ccn_map_full = fread(
    ccn_map_path,
    colClasses = c(hospital_ccn_id = "character",
                   hospital_id     = "character")
  )
  
  # required columns
  required_cols = c("site", "hospital_id", "hospital_ccn_id")
  missing_cols = setdiff(required_cols, names(ccn_map_full))
  if (length(missing_cols)) {
    stop("Hospital CCN map is missing required columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  # standardize column names
  setnames(ccn_map_full, "hospital_ccn_id", "medicare_id")
  
  # ensure ccn_note column exists
  if (!"ccn_note" %in% names(ccn_map_full)) {
    ccn_map_full[, ccn_note := NA_character_]
  }
  
  # zero-pad CCNs to 6 digits (preserve leading zeros for CA, AK, etc.)
  ccn_map_full[, medicare_id := stringr::str_pad(
    medicare_id, width = 6, side = "left", pad = "0"
  )]
  
  # case-insensitive site name match
  target_site = tolower(trimws(site_lowercase))
  ccn_map_full[, site_norm := tolower(trimws(site))]
  
  hospital_medicare_map = ccn_map_full[
    site_norm == target_site,
    .(hospital_id, medicare_id, ccn_note)
  ]
  
  if (nrow(hospital_medicare_map) == 0) {
    available = sort(unique(ccn_map_full$site))
    warning("No CCN map entries for site '", site_lowercase, "'.",
            "\n  Available sites in map: ", paste(available, collapse = ", "),
            "\n  Set 'site_lowercase' or 'institution' in config.yaml ",
            "to one of these values (case-insensitive).", call. = FALSE)
    hospital_medicare_map = NULL
    
  } else {
    n_hosp = nrow(hospital_medicare_map)
    n_ccn  = uniqueN(hospital_medicare_map$medicare_id)
    message(sprintf("CCN map loaded: %d hospital_id(s) → %d unique CCN(s) for site '%s'",
                    n_hosp, n_ccn, site_lowercase))
    
    # report shared CCNs (multiple hospital_ids with same CCN)
    shared = hospital_medicare_map[, .N, by = medicare_id][N > 1]
    if (nrow(shared) > 0) {
      message(sprintf("  %d shared CCN(s) — these hospital_ids will be aggregated:",
                      nrow(shared)))
      for (i in seq_len(nrow(shared))) {
        hids = hospital_medicare_map[medicare_id == shared$medicare_id[i], hospital_id]
        message(sprintf("    %s ← %s", shared$medicare_id[i],
                        paste(hids, collapse = ", ")))
      }
    }
  }
  
  rm(ccn_map_full)
}


## file_type must be one of c(parquet, csv, fst) -------------------------------

if (!file_type %in% allowed_files) {
  stop(
    "Invalid '", file_type,
    "'. Expected one of: ",
    paste(allowed_files, collapse = ", "),
    call. = FALSE
  )
}


## create output directories ---------------------------------------------------

# Site-named output directories for the federated workflow:
#   {site}_private_data   — patient-level intermediate files (NOT shared)
#   {site}_upload_to_box  — deidentified shareable files for cross-site compilation

private_dir   = file.path(project_location, paste0(site_lowercase, "_private_data"))
shareable_dir = file.path(project_location, paste0(site_lowercase, "_upload_to_box"))

for (d in c(private_dir, shareable_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Created: ", d)
  }
}

message(sprintf("Private data:   %s", private_dir))
message(sprintf("Shareable data: %s", shareable_dir))


## dates -----------------------------------------------------------------------

start_date = with_tz(ymd_hms("2016-01-01 00:00:00", tz = tz_local), "UTC")
end_date   = with_tz(ymd_hms("2024-12-31 23:59:59", tz = tz_local), "UTC")


## pandemic period (protocol v6.0 §3.3) ----------------------------------------

# Three calendar-year eras (NOT US PHE dates):
#   Pre-pandemic:   2016-01-01 to 2019-12-31
#   Acute pandemic: 2020-01-01 to 2021-12-31
#   Post-pandemic:  2022-01-01 to 2024-12-31

assign_pandemic_period = function(dttm) {
  yr = year(dttm)
  factor(
    fcase(
      yr <= 2019, "Pre-pandemic",
      yr <= 2021, "Acute pandemic",
      default    = "Post-pandemic"
    ),
    levels = c("Pre-pandemic", "Acute pandemic", "Post-pandemic")
  )
}

message("Pandemic eras: Pre (2016-2019) | Acute (2020-2021) | Post (2022-2024)")
today      = format(Sys.Date(), "%y%m%d")


## standardized save function --------------------------------------------------

# shareable = FALSE → patient-level data → {site}_private_data/
# shareable = TRUE  → deidentified data  → {site}_upload_to_box/
#
# IMPORTANT: shareable files must NOT contain patient_id, hospitalization_id,
# or any other direct identifiers. Only aggregated counts, rates, model
# coefficients, and hospital-level (by CCN) summaries are shareable.

save_output = function(obj, name, shareable = FALSE, subfolder = NULL, format = "csv") {
  
  base_dir = if (shareable) shareable_dir else private_dir
  if (!is.null(subfolder)) base_dir = file.path(base_dir, subfolder)
  
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
  
  path = file.path(base_dir, paste0(name, "_", site_lowercase, "_", today, ".", format))
  
  switch(format,
         csv     = fwrite(obj, path),
         fst     = fst::write_fst(obj, path),
         parquet = arrow::write_parquet(obj, path),
         stop("Unsupported format: ", format, call. = FALSE)
  )
  
  message("Saved: ", path)
  invisible(path)
}

# helper for direct fwrite/ggsave calls in aim scripts
private_path  = function(filename) file.path(private_dir, filename)
shareable_path = function(filename) file.path(shareable_dir, filename)
# newest date-stamped output matching `pattern` in `dir` (files end _YYMMDD.csv,
# so a plain descending sort is chronological). Errors if none found; callers
# treating a file as optional should use list.files() + a length() guard instead.
latest_output = function(pattern, dir = private_dir) {
  f = list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(f)) stop("No files match '", pattern, "' in ", dir)
  f[order(f, decreasing = TRUE)][1]
}

# ---- QI data-quality flags --------------------------------------------------
# Surface data-quality problems the pipeline otherwise drops silently, so every
# site sees them at run time. Prints a visible banner to the console; it does not
# stop the run.
qi_flag = function(title, detail) {
  bar = strrep("=", 74)
  message("\n  ", bar)
  message("  [QI DATA-QUALITY FLAG] ", title)
  for (ln in strsplit(detail, "\n", fixed = TRUE)[[1]]) message("    ", ln)
  message("  ", bar, "\n")
}



## colorblind-safe palettes and standard theme ---------------------------------

# Okabe-Ito palette (widely recommended for accessibility)
# Named variants for specific use cases across all analysis scripts.

pal_okabe_ito = c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#000000"   # black
)

# outcome groups (used in KM, temporal trends, trajectory comparisons)
pal_outcome = c(
  "Same-facility hospice" = "#0072B2",   # blue
  "In-hospital death"     = "#D55E00",   # vermillion
  "Acute death"           = "#D55E00",   # vermillion (alias)
  "Routine discharge"     = "#009E73",   # bluish green
  "Transfer"              = "#E69F00",   # orange
  "External hospice"      = "#56B4E9",   # sky blue
  "Other/AMA"             = "#CC79A7"    # reddish purple
)

# LST groups
pal_lst = c(
  "LST"     = "#D55E00",   # vermillion
  "Non-LST" = "#0072B2"    # blue
)

# pandemic periods (3 eras per §3.3)
pal_pandemic = c(
  "Pre-pandemic"   = "#0072B2",   # blue
  "Acute pandemic" = "#D55E00",   # vermillion
  "Post-pandemic"  = "#009E73"    # bluish green
)

# for multi-site overlay (up to 12 sites)
pal_sites = c(
  "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9",
  "#CC79A7", "#F0E442", "#000000", "#999999", "#661100",
  "#332288", "#88CCEE"
)

# standard ggplot theme for all figures
theme_clif = function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}


## safe division and log -------------------------------------------------------

safe_div  = function(num, den) fifelse(!is.na(num) & !is.na(den) & den > 0, num / den, NA_real_)
safe_log2 = function(x) fifelse(!is.na(x) & x > 0, log2(x), NA_real_)

## standardize for trajectory fitting ------------------------------------------

standardize = function(x) {
  mu = fmean(x, na.rm = TRUE)
  s  = fsd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}


# clif tables needed for project -----------------------------------------------

required_tables = c(
  "patient",
  "hospitalization",
  "adt",
  "labs",
  "vitals",
  "hospital_diagnosis",
  "respiratory_support",
  "medication_admin_continuous",
  "code_status"
)

## check required tables against available tables ------------------------------

clif_table_filenames = list.files(
  path       = tables_location,
  pattern    = paste0("^clif_.*\\.", file_type, "$"),
  full.names = TRUE
)

clif_table_basenames =
  basename(clif_table_filenames) |>
  str_remove(paste0("\\.", file_type, "$")) |>
  str_remove("^clif_") |>
  str_remove("(_\\d{4}(_\\d{4})*)$")

table_file_map     = setNames(clif_table_filenames, clif_table_basenames)
missing_tables     = setdiff(required_tables, clif_table_basenames)
required_filenames = table_file_map[required_tables]

if (length(missing_tables)) {
  stop("Missing required tables: ", paste(missing_tables, collapse = ", "), call. = FALSE)
}
message("All required tables are present.")

rm(allowed_files, missing_tables)
gc()


# load tables ------------------------------------------------------------------

data_list = switch(file_type,
                   parquet = lapply(required_filenames, arrow::open_dataset),
                   csv     = lapply(required_filenames, arrow::read_csv_arrow),
                   fst     = lapply(required_filenames, \(f) arrow::arrow_table(fst::read_fst(f, as.data.table = TRUE))),
                   stop("Unsupported file format", call. = FALSE)
)
names(data_list) = names(required_filenames)


# normalize identifier key types -----------------------------------------------
# Parquet id columns may be utf8 at one site and large_utf8 at another; arrow
# refuses to join keys whose string types differ (e.g. adt hospitalization_id as
# large_string vs a utf8 lookup table). Cast every id key to plain utf8 at load
# so downstream joins are type-consistent regardless of how a site wrote its
# parquet, and so collected id columns land as character for data.table merges.

cast_id_keys = function(tbl) {
  nm = names(tbl)
  if ("hospitalization_id" %in% nm)
    tbl = dplyr::mutate(tbl, hospitalization_id = cast(hospitalization_id, arrow::string()))
  if ("patient_id" %in% nm)
    tbl = dplyr::mutate(tbl, patient_id = cast(patient_id, arrow::string()))
  if ("hospital_id" %in% nm)
    tbl = dplyr::mutate(tbl, hospital_id = cast(hospital_id, arrow::string()))
  tbl
}

data_list = lapply(data_list, cast_id_keys)


# validate table contents ------------------------------------------------------

validate_table = function(tbl, table_name, req_vars = NULL, req_values = list()) {
  
  req_vars = req_vars %||% character(0)
  tbl_cols = names(tbl)
  problems = character()
  
  # check required variables
  missing_vars = setdiff(req_vars, tbl_cols)
  if (length(missing_vars)) {
    problems = c(problems, sprintf("Missing required vars: %s", paste(missing_vars, collapse = ", ")))
  }
  
  # check required values
  for (var in names(req_values)) {
    if (!var %in% tbl_cols) {
      problems = c(problems, sprintf("Missing '%s' needed for value checks.", var))
      next
    }
    
    present_vals = tryCatch({
      sampled = tbl |>
        dplyr::select(!!rlang::sym(var)) |>
        dplyr::distinct() |>
        head(5000) |>
        dplyr::collect()
      tolower(trimws(na.omit(as.character(sampled[[var]]))))
    }, error = \(e) character(0))
    
    expected = tolower(trimws(unique(req_values[[var]])))
    missing  = setdiff(expected, present_vals)
    
    if (length(missing)) {
      problems = c(problems, sprintf(
        "Variable '%s' is missing expected values: %s", var, paste(missing, collapse = ", ")
      ))
    }
  }
  
  if (length(problems)) {
    return(sprintf("Table '%s':\n- %s", table_name, paste(problems, collapse = "\n- ")))
  }
  invisible(NULL)
}

validate_all_tables = function(data_list, specs) {
  all_problems = vapply(specs, \(spec) {
    if (!spec$table_name %in% names(data_list)) {
      return(sprintf("Table '%s' is missing entirely.", spec$table_name))
    }
    result = validate_table(
      tbl        = data_list[[spec$table_name]],
      table_name = spec$table_name,
      req_vars   = spec$req_vars,
      req_values = spec$req_values
    )
    result %||% ""
  }, character(1))
  
  all_problems = all_problems[nzchar(all_problems)]
  
  if (length(all_problems)) {
    stop("Validation errors found:\n", paste(all_problems, collapse = "\n\n"), call. = FALSE)
  }
  message("Validation passed: all required tables present with needed values.")
}


## validation specs ------------------------------------------------------------

req_vitals = c(
  "heart_rate", 
  "respiratory_rate", 
  "map", 
  "spo2", 
  "temp_c",
  "height_cm", 
  "weight_kg"
)

req_labs = c(
  "crp", 
  "wbc", 
  "creatinine", 
  "platelet_count",
  "lymphocytes_absolute", 
  "ldh", 
  "ferritin",
  "hemoglobin", 
  "neutrophils_absolute",
  "alt", 
  "ast", 
  "bilirubin_total", 
  "ph_arterial", 
  "sodium", 
  "potassium", 
  "po2_arterial", 
  "pco2_arterial")

validation_specs = list(
  list(
    table_name = "patient",
    req_vars   = c("patient_id", 
                   "race_category", 
                   "ethnicity_category",
                   "sex_category", 
                   "death_dttm")
  ),
  list(
    table_name = "hospitalization",
    req_vars   = c("patient_id", 
                   "hospitalization_id", 
                   "age_at_admission",
                   "admission_dttm", 
                   "discharge_dttm")
  ),
  list(
    table_name = "adt",
    req_vars   = c("hospitalization_id", 
                   "hospital_id",
                   "location_category", 
                   "in_dttm", 
                   "out_dttm")
  ),
  list(
    table_name = "medication_admin_continuous",
    req_vars   = c("hospitalization_id", 
                   "med_group", 
                   "med_category"),
    req_values = list(
      med_group    = "vasoactives",
      med_category = c("norepinephrine", 
                       "vasopressin", 
                       "epinephrine", 
                       "phenylephrine")
    )
  ),
  list(
    table_name = "respiratory_support",
    req_vars   = c("hospitalization_id", 
                   "device_category", 
                   "recorded_dttm"),
    req_values = list(device_category = c("IMV", 
                                          "NIPPV", 
                                          "High Flow NC"))
  ),
  list(
    table_name = "vitals",
    req_vars   = c("hospitalization_id", 
                   "vital_category", 
                   "vital_value", 
                   "recorded_dttm"),
    req_values = list(vital_category = req_vitals)
  ),
  list(
    table_name = "labs",
    req_vars   = c("hospitalization_id", 
                   "lab_category", 
                   "lab_value", 
                   "lab_result_dttm"),
    req_values = list(lab_category = req_labs)
  ),
  list(
    table_name = "hospital_diagnosis",
    req_vars   = c("hospitalization_id", 
                   "diagnosis_code", 
                   "diagnosis_code_format")
  ),
  list(
    table_name = "code_status",
    req_vars   = c("patient_id", 
                   "start_dttm", 
                   "code_status_category")
  ))


## run validation --------------------------------------------------------------

validate_all_tables(data_list, validation_specs)

## optional: patient_assessments (GCS for APACHE II) ---------------------------
# Used only to source GCS for APACHE II, which defaults GCS to 15 when this table
# is absent or unusable. A missing or malformed patient_assessments table must
# therefore NOT stop the pipeline; it is optional, like crrt_therapy below.
has_patient_assessments = "patient_assessments" %in% clif_table_basenames

if (has_patient_assessments) {
  
  pa_file = table_file_map[["patient_assessments"]]
  
  data_list[["patient_assessments"]] = cast_id_keys(switch(
    file_type,
    parquet = arrow::open_dataset(pa_file),
    csv     = arrow::read_csv_arrow(pa_file),
    fst     = arrow::arrow_table(fst::read_fst(pa_file, as.data.table = TRUE))
  ))
  
  pa_validation = validate_table(
    tbl        = data_list[["patient_assessments"]],
    table_name = "patient_assessments",
    req_vars   = c("hospitalization_id", "recorded_dttm",
                   "assessment_category", "numerical_value"),
    req_values = list(assessment_category = "gcs_total")
  )
  
  if (!is.null(pa_validation)) {
    warning("patient_assessments present but failed validation:\n",
            pa_validation,
            "\n  GCS will be unavailable; APACHE II will default GCS to 15.", call. = FALSE)
    has_patient_assessments = FALSE
    data_list[["patient_assessments"]] = NULL
  } else {
    message("Optional table patient_assessments loaded and validated.")
  }
  
} else {
  message("Optional table patient_assessments not found; APACHE II will default GCS to 15.")
}

## optional: crrt_therapy ------------------------------------------------------

has_crrt = "crrt_therapy" %in% clif_table_basenames

if (has_crrt) {
  
  crrt_file = table_file_map[["crrt_therapy"]]
  
  data_list[["crrt_therapy"]] = cast_id_keys(switch(
    file_type,
    parquet = arrow::open_dataset(crrt_file),
    csv     = arrow::read_csv_arrow(crrt_file),
    fst     = arrow::arrow_table(fst::read_fst(crrt_file, as.data.table = TRUE))
  ))
  
  # validate crrt_therapy
  crrt_validation = validate_table(
    tbl        = data_list[["crrt_therapy"]],
    table_name = "crrt_therapy",
    req_vars   = c("hospitalization_id", "recorded_dttm", "crrt_mode_category")
  )
  
  if (!is.null(crrt_validation)) {
    warning("crrt_therapy table present but failed validation:\n",
            crrt_validation,
            "\n  CRRT discontinuation flags will be skipped.", call. = FALSE)
    has_crrt = FALSE
    data_list[["crrt_therapy"]] = NULL
  } else {
    message("Optional table crrt_therapy loaded and validated.")
  }
  
} else {
  message("Optional table crrt_therapy not found — CRRT flags will be skipped.")
}

# clean up ---------------------------------------------------------------------

rm(validation_specs, req_labs, req_vitals,
   clif_table_basenames, clif_table_filenames,
   required_filenames, required_tables,
   table_file_map, validate_table, validate_all_tables,
   cast_id_keys)
gc()

message("00_setup.R complete.")
