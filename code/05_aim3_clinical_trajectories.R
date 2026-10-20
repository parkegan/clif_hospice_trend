

# 05_aim3_clinical_trajectories.R
# Aim 3: Clinical trajectories preceding hospice transition
# Protocol section 3.6.4

# Requires updated 02_clinical_inflection.R

# Analyses:
#   3.1  Overall trajectory characterization (hospice transitions)
#   3.2  Transitions from life-sustaining therapies (LST subcohort)
#   3.3  Hospital-level variation in LST-to-hospice patterns
#   3.4  LST vs non-LST comparison

# Outputs (shareable CSVs in final_output/):

#   aim3_trajectory_summary_{site}.csv
#   aim3_lst_disposition_{site}.csv
#   aim3_lst_timing_{site}.csv
#   aim3_code_status_sequence_{site}.csv
#   aim3_hospital_lst_hospice_rate_{site}.csv
#   aim3_lst_vs_nonlst_km_{site}.csv
#   aim3_lst_vs_nonlst_trajectories_{site}.csv


# Outputs (local PDFs):

#   aim3_lab_trajectory_{site}.pdf
#   aim3_opioid_trajectory_{site}.pdf
#   aim3_lst_timing_boxplot_{site}.pdf
#   aim3_event_sequence_{site}.pdf
#   aim3_lst_caterpillar_{site}.pdf
#   aim3_lst_vs_nonlst_km_{site}.pdf


# Last updated 26.10.19: small cell suppression in the trajectory summary now
# applies to COUNT rows only. The median_count and median_hours rows carry a
# continuous summary (a median) in the n column, not a patient count, and the
# UMN 261015 bundle shipped medians of 2 to 4 masked as "<5". Those rows now
# pass through unmasked; the flag and category count rows are suppressed
# exactly as before.
#
# Previous update bpe 27.05.2026


# 0. setup ---------------------------------------------------------------------

if (!exists("data_list")) source(here::here("code", "00_setup.R"))

if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")

library(ggplot2)
library(survival)
library(patchwork)
library(scales)


# 1. load data -----------------------------------------------------------------

message("\n===== 05_aim3: Clinical Trajectories Preceding Hospice Transition =====\n")

cohort_all = fread(latest_output("cohort_all_hospitalizations"),
                   colClasses = list(character = c("patient_id", "hospitalization_id")))
cohort_all[, hospitalization_id := as.character(hospitalization_id)]

# dying heuristic flags
dying = fread(latest_output("clinical_inflection_flags"))
dying[, hospitalization_id := as.character(hospitalization_id)]

# covariates (from 06)
cov_file = list.files(
  private_dir,
  pattern = "covariates_aim2a3", full.names = TRUE
)
if (length(cov_file) > 0) {
  covariates = fread(cov_file[order(cov_file, decreasing = TRUE)][1])
  covariates[, hospitalization_id := as.character(hospitalization_id)]
} else {
  covariates = NULL
}

# primary cohort
cohort = cohort_all[in_cohort == 1L]

# ensure datetimes
dt_cols = c("admission_dttm", "discharge_dttm", "death_dttm", "hospice_transition_dttm")
for (col in dt_cols) {
  if (col %in% names(cohort) && !inherits(cohort[[col]], "POSIXct")) {
    set(cohort, j = col, value = as.POSIXct(cohort[[col]], tz = "UTC"))
  }
}

# endpoint
cohort[, endpoint_dttm := fifelse(
  !is.na(death_dttm) & death_dttm <= discharge_dttm + lubridate::days(30),
  death_dttm, discharge_dttm
)]

# merge dying heuristic
cohort = merge(cohort, dying, by = "hospitalization_id", all.x = TRUE)

# merge covariates
if (!is.null(covariates)) {
  cov_cols = setdiff(names(covariates), names(cohort))
  cov_cols = c("hospitalization_id", cov_cols)
  cohort = merge(cohort, covariates[, ..cov_cols], by = "hospitalization_id", all.x = TRUE)
}

# fill NA heuristic flags with 0
flag_cols = grep("^flag_|^inflection|^n_core|^n_all", names(cohort), value = TRUE)
for (col in flag_cols) set(cohort, which(is.na(cohort[[col]])), col, 0L)

# defensive LST flag check
for (col in c("imv_flag", "vaso_flag", "lst_any_flag"))
  if (!col %in% names(cohort)) set(cohort, j = col, value = 0L)

# hospital unit
if (!is.null(hospital_medicare_map)) {
  if (!"hospital_unit" %in% names(cohort)) {
    cohort = merge(cohort,
                   hospital_medicare_map[, .(hospital_id, medicare_id)],
                   by = "hospital_id", all.x = TRUE)
    cohort[is.na(medicare_id), medicare_id := hospital_id]
    cohort[, hospital_unit := medicare_id]
  }
} else {
  if (!"hospital_unit" %in% names(cohort))
    cohort[, hospital_unit := as.character(hospital_id)]
}

# derived flags
cohort[, any_hospice := fifelse(
  outcome_category %in% c("same_facility_hospice_death",
                           "same_facility_hospice_alive",
                           "external_hospice_discharge"), 1L, 0L
)]
cohort[, died := fifelse(
  outcome_category %in% c("in_hospital_death_acute",
                           "same_facility_hospice_death",
                           "external_hospice_death",
                           "post_discharge_death_non_hospice"), 1L, 0L
)]

n_cohort = nrow(cohort)
n_lst = sum(cohort$lst_any_flag == 1L)
n_hospice = sum(cohort$any_hospice == 1L)

message(sprintf("Cohort: %s | LST subcohort: %s | Any hospice: %s",
                format(n_cohort, big.mark = ","),
                format(n_lst, big.mark = ","),
                format(n_hospice, big.mark = ",")))

rm(dying, covariates)


# 2. Verify Aim 3 timing variables from 02_clinical_inflection -----------------


message("\n--- Verifying Aim 3 timing variables from dying heuristic ---")

aim3_timing_cols = c(
  "last_vaso_dttm", "last_imv_dttm",
  "first_lst_withdrawal_dttm",
  "hours_first_withdrawal_to_hospice", "hours_first_withdrawal_to_death",
  "hours_vaso_to_hospice", "hours_imv_to_hospice",
  "hours_comfort_to_hospice", "hours_comfort_to_withdrawal",
  "first_comfort_dttm", "code_status_sequence",
  "labs_72_48h", "labs_48_24h", "labs_24_0h",
  "opioid_72_48h", "opioid_48_24h", "opioid_24_0h", "opioid_trajectory"
)

present = intersect(aim3_timing_cols, names(cohort))
missing = setdiff(aim3_timing_cols, names(cohort))

if (length(missing) > 0) {
  message(sprintf("  WARNING: %d Aim 3 timing columns missing from dying heuristic output.",
                  length(missing)))
  message(sprintf("    Missing: %s", paste(missing, collapse = ", ")))
  message("  These analyses may be incomplete. Re-run 02_clinical_inflection.R")
  message("  with the updated version that computes timing variables.")

  # fill numeric missings with NA, integer counts with 0
  for (col in missing) {
    if (grepl("hours_|dttm", col)) {
      set(cohort, j = col, value = NA_real_)
    } else if (grepl("labs_|opioid_[0-9]", col)) {
      set(cohort, j = col, value = 0L)
    } else if (col == "code_status_sequence") {
      set(cohort, j = col, value = "Unknown")
    } else if (col == "opioid_trajectory") {
      set(cohort, j = col, value = "None")
    } else {
      set(cohort, j = col, value = NA_real_)
    }
  }
} else {
  message(sprintf("  All %d Aim 3 timing columns present.", length(present)))
}

# fill NA windowed counts with 0
for (col in c("labs_72_48h", "labs_48_24h", "labs_24_0h",
              "opioid_72_48h", "opioid_48_24h", "opioid_24_0h")) {
  if (col %in% names(cohort))
    set(cohort, which(is.na(cohort[[col]])), col, 0L)
}
if ("opioid_trajectory" %in% names(cohort))
  cohort[is.na(opioid_trajectory), opioid_trajectory := "None"]



# 3. Analysis 3.1: Overall trajectories ----------------------------------------

message("\n--- Analysis 3.1: Overall trajectory characterization ---")

# restrict to same-facility hospice transitions
hospice_sf = cohort[outcome_category %in% c("same_facility_hospice_death",
                                              "same_facility_hospice_alive")]
hospice_sf[, survival_group := fifelse(
  short_interval_death == 1L, "Short-survival (<48h)", "Standard"
)]

n_hospice_sf = nrow(hospice_sf)
message(sprintf("  Same-facility hospice transitions: %s", format(n_hospice_sf, big.mark = ",")))

# component prevalence by survival group
traj_components = c(
  "flag_code_status_strict", "flag_code_status_lenient",
  "flag_opioid", "flag_opioid_initiation", "flag_opioid_escalation",
  "flag_lab_cessation",
  "flag_withdrawal", "flag_vaso_withdrawal", "flag_mv_withdrawal",
  # protocol v7.0: the composites are DEMOTED to diagnostics. The four core
  # signals above are reported separately; these two are retained only as a
  # data-quality cross-check, tagged signal_role = "diagnostic_composite".
  "inflection_2of4", "inflection_1of4"
)
traj_components = intersect(traj_components, names(hospice_sf))

# also include non-hospice deaths for comparison
acute_deaths = cohort[outcome_category == "in_hospital_death_acute"]
routine = cohort[outcome_category == "routine_discharge"]

trajectory_summary = rbindlist(lapply(
  list(
    list(dt = hospice_sf, group = "Same-facility hospice (all)"),
    list(dt = hospice_sf[survival_group == "Short-survival (<48h)"],
         group = "Hospice: short-survival"),
    list(dt = hospice_sf[survival_group == "Standard"],
         group = "Hospice: standard"),
    list(dt = acute_deaths, group = "Acute in-hospital death"),
    list(dt = routine, group = "Routine discharge")
  ), \(x) {
    dt = x$dt; grp = x$group; n = nrow(dt)
    if (n == 0) return(NULL)

    # component prevalence
    comp_rows = rbindlist(lapply(traj_components, \(f) {
      data.table(
        variable = f, level = "flagged",
        n = sum(dt[[f]] == 1L, na.rm = TRUE),
        pct = round(100 * mean(dt[[f]] == 1L, na.rm = TRUE), 2)
      )
    }))

    # opioid trajectory distribution
    opi_tab = dt[, .N, by = opioid_trajectory]
    opi_rows = opi_tab[, .(variable = "opioid_trajectory", level = opioid_trajectory,
                           n = N, pct = round(100 * N / n, 2))]

    # code status sequence distribution (LST patients only)
    if ("code_status_sequence" %in% names(dt)) {
      cs_tab = dt[lst_any_flag == 1L & !is.na(code_status_sequence),
                  .N, by = code_status_sequence]
      if (nrow(cs_tab) > 0) {
        n_cs = sum(cs_tab$N)
        cs_rows = cs_tab[, .(variable = "code_status_sequence",
                             level = code_status_sequence,
                             n = N, pct = round(100 * N / n_cs, 2))]
      } else {
        cs_rows = NULL
      }
    } else {
      cs_rows = NULL
    }

    # lab trajectory medians
    lab_rows = data.table(
      variable = c("labs_72_48h", "labs_48_24h", "labs_24_0h"),
      level    = "median_count",
      n        = c(median(dt$labs_72_48h), median(dt$labs_48_24h), median(dt$labs_24_0h)),
      pct      = NA_real_
    )

    # opioid count trajectory medians
    opi_count_rows = data.table(
      variable = c("opioid_72_48h", "opioid_48_24h", "opioid_24_0h"),
      level    = "median_count",
      n        = c(median(dt$opioid_72_48h), median(dt$opioid_48_24h), median(dt$opioid_24_0h)),
      pct      = NA_real_
    )

    # timing of code status change relative to hospice (median hours)
    timing_rows = NULL
    if (sum(!is.na(dt$hours_comfort_to_hospice)) >= 5) {
      timing_rows = data.table(
        variable = "hours_comfort_to_hospice",
        level    = "median_hours",
        n        = round(median(dt$hours_comfort_to_hospice, na.rm = TRUE), 1),
        pct      = NA_real_
      )
    }

    out = rbind(comp_rows, opi_rows, cs_rows, lab_rows, opi_count_rows, timing_rows,
                fill = TRUE)
    out[, `:=`(comparison_group = grp, site = site_lowercase, n_group = n)]
    out
  }
))

trajectory_summary[, signal_role := fifelse(grepl("^inflection_", variable),
                                            "diagnostic_composite", "signal")]
setcolorder(trajectory_summary,
            c("site", "comparison_group", "n_group", "variable", "signal_role",
              "level", "n", "pct"))

# print key results -- the four core signals reported SEPARATELY (protocol v7.0),
# with the composite shown only as a trailing diagnostic.
core_signals = c(`code-status change` = "flag_code_status_strict",
                 `opioid`             = "flag_opioid",
                 `lab cessation`      = "flag_lab_cessation",
                 `LST withdrawal`     = "flag_withdrawal")
pct_of = function(sub, v) { r = sub[variable == v]; if (nrow(r)) sprintf("%s", r$pct[1]) else "N/A" }
message("\n===== DOCUMENTED CARE-TRANSITION SIGNALS (reported separately) =====")
for (grp in unique(trajectory_summary$comparison_group)) {
  sub = trajectory_summary[comparison_group == grp]
  n_g = sub$n_group[1]
  message(sprintf("  %s (n = %s):", grp, format(n_g, big.mark = ",")))
  message(sprintf("    code-status change %s%% | opioid %s%% | lab cessation %s%% | LST withdrawal %s%%",
                  pct_of(sub, "flag_code_status_strict"), pct_of(sub, "flag_opioid"),
                  pct_of(sub, "flag_lab_cessation"),      pct_of(sub, "flag_withdrawal")))
  message(sprintf("    [diagnostic] >=2/4 composite: %s%%", pct_of(sub, "inflection_2of4")))
}
message("====================================================================\n")

# FIXED 2026.10.19: suppression applies to COUNT rows only. The median_count
# and median_hours rows carry a MEDIAN in the n column (median labs or opioid
# administrations per window, median hours from comfort order to hospice),
# not a patient count: a median of 3 labs is not a small cell, and masking it
# as "<5" destroyed the trajectory signal at UMN (261015). The count rows
# (flag prevalences and category distributions) are suppressed exactly as
# before; the median rows pass through unmasked with suppressed = 0. The row
# order of the written file is unchanged.
if (isTRUE(suppress_small_cells_enabled)) {
  trajectory_summary[, .row_order := .I]
  is_median = grepl("^median", trajectory_summary$level)
  ts_counts = suppress_small_cells(trajectory_summary[!is_median],
                                   count_cols = "n", pct_cols = "pct")
  ts_median = copy(trajectory_summary[is_median])
  ts_median[, n := as.character(n)]   # per element, so 3 prints "3", not "3.0"
  ts_median[, suppressed := 0L]
  trajectory_summary = rbind(ts_counts, ts_median, use.names = TRUE)
  setorder(trajectory_summary, .row_order)
  trajectory_summary[, .row_order := NULL]
}
fwrite(
  trajectory_summary,
  shareable_path(
            paste0("aim3_trajectory_summary_", site_lowercase, "_", today, ".csv"))
)
message("  Saved aim3_trajectory_summary.")


## 3.1b: Lab trajectory plot ---------------------------------------------------

message("  Generating lab trajectory plot...")

lab_plot_dt = rbindlist(lapply(
  list(
    list(dt = hospice_sf, group = "Same-facility hospice"),
    list(dt = acute_deaths, group = "Acute death"),
    list(dt = routine, group = "Routine discharge")
  ), \(x) {
    dt = x$dt; n = nrow(dt)
    if (n < 10) return(NULL)
    data.table(
      group  = x$group,
      window = factor(c("72-48h", "48-24h", "24-0h"),
                      levels = c("72-48h", "48-24h", "24-0h")),
      median_labs = c(
        median(dt$labs_72_48h), median(dt$labs_48_24h), median(dt$labs_24_0h)
      ),
      mean_labs = c(
        mean(dt$labs_72_48h), mean(dt$labs_48_24h), mean(dt$labs_24_0h)
      )
    )
  }
))

if (nrow(lab_plot_dt) > 0) {
  p_labs = ggplot(lab_plot_dt, aes(x = window, y = mean_labs,
                                    color = group, group = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "Same-facility hospice" = pal_outcome[["Same-facility hospice"]],
      "Acute death"           = pal_outcome[["Acute death"]],
      "Routine discharge"     = pal_outcome[["Routine discharge"]]
    )) +
    labs(
      x     = "Time window before endpoint",
      y     = "Mean routine lab count",
      title = "Laboratory ordering trajectory in final 72 hours",
      subtitle = paste0(toupper(site_lowercase), " | By disposition group"),
      color = NULL
    ) +
    theme_clif() +
    theme(
    )

  ggsave(
    shareable_path(
              paste0("aim3_lab_trajectory_", site_lowercase, "_", today, ".pdf")),
    plot = p_labs, width = 8, height = 6, dpi = 300
  )
  message("  Saved aim3_lab_trajectory.")

  # companion data for cross-site re-plotting
  lab_plot_dt[, site := site_lowercase]
  fwrite(lab_plot_dt, shareable_path(
    paste0("aim3_lab_trajectory_data_", site_lowercase, "_", today, ".csv")))
}


## 3.1c: Opioid trajectory plot -----------------------------------------------

message("  Generating opioid trajectory plot...")

opi_plot_dt = rbindlist(lapply(
  list(
    list(dt = hospice_sf, group = "Same-facility hospice"),
    list(dt = acute_deaths, group = "Acute death")
  ), \(x) {
    dt = x$dt; n = nrow(dt)
    if (n < 10) return(NULL)
    data.table(
      group  = x$group,
      window = factor(c("72-48h", "48-24h", "24-0h"),
                      levels = c("72-48h", "48-24h", "24-0h")),
      mean_count = c(
        mean(dt$opioid_72_48h), mean(dt$opioid_48_24h), mean(dt$opioid_24_0h)
      )
    )
  }
))

if (nrow(opi_plot_dt) > 0) {
  p_opioid = ggplot(opi_plot_dt, aes(x = window, y = mean_count,
                                      color = group, group = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "Same-facility hospice" = pal_outcome[["Same-facility hospice"]],
      "Acute death"           = pal_outcome[["Acute death"]]
    )) +
    labs(
      x     = "Time window before endpoint",
      y     = "Mean opioid administrations",
      title = "Opioid administration trajectory in final 72 hours",
      subtitle = paste0(toupper(site_lowercase), " | Hospice vs acute death"),
      color = NULL
    ) +
    theme_clif() +
    theme(
    )

  ggsave(
    shareable_path(
              paste0("aim3_opioid_trajectory_", site_lowercase, "_", today, ".pdf")),
    plot = p_opioid, width = 8, height = 6, dpi = 300
  )
  message("  Saved aim3_opioid_trajectory.")

  # companion data for cross-site re-plotting
  opi_plot_dt[, site := site_lowercase]
  fwrite(opi_plot_dt, shareable_path(
    paste0("aim3_opioid_trajectory_data_", site_lowercase, "_", today, ".csv")))
}

rm(hospice_sf, acute_deaths, routine, lab_plot_dt, opi_plot_dt)


# 4. Analysis 3.2: Transitions from life-sustaining therapies ------------------

message("\n--- Analysis 3.2: LST transitions ---")

lst = cohort[lst_any_flag == 1L]
n_lst = nrow(lst)
message(sprintf("  LST subcohort: %s encounters.", format(n_lst, big.mark = ",")))


## 4a. LST disposition flow ----------------------------------------------------

# classify LST type
lst[, lst_type := fcase(
  imv_flag == 1L & vaso_flag == 1L, "IMV+Vaso",
  imv_flag == 1L,                   "IMV only",
  vaso_flag == 1L,                  "Vaso only",
  default = "Other"
)]

lst_disposition = lst[, .(N = .N), by = .(lst_type, outcome_category)]
lst_disposition[, pct := round(100 * N / sum(N), 2), by = lst_type]
lst_disposition[, `:=`(site = site_lowercase, n_lst = n_lst)]

# also compute any-LST level
lst_any_disp = lst[, .(N = .N), by = outcome_category]
lst_any_disp[, `:=`(
  lst_type = "Any LST",
  pct = round(100 * N / n_lst, 2),
  site = site_lowercase,
  n_lst = n_lst
)]

lst_disposition = rbind(lst_disposition, lst_any_disp)
setcolorder(lst_disposition, c("site", "lst_type", "n_lst", "outcome_category", "N", "pct"))

message("\n  LST disposition (any LST):")
any_disp = lst_disposition[lst_type == "Any LST"][order(-N)]
for (i in seq_len(nrow(any_disp))) {
  message(sprintf("    %-40s %s (%4.1f%%)",
                  any_disp$outcome_category[i],
                  format(any_disp$N[i], big.mark = ","),
                  any_disp$pct[i]))
}

fwrite(
  lst_disposition,
  shareable_path(
            paste0("aim3_lst_disposition_", site_lowercase, "_", today, ".csv"))
)
message("  Saved aim3_lst_disposition.")


## 4b. LST withdrawal → hospice timing -----------------------------------------

message("\n  Computing LST withdrawal timing intervals...")

lst_hospice = lst[any_hospice == 1L & !is.na(hospice_transition_dttm)]
n_lst_hospice = nrow(lst_hospice)
message(sprintf("    LST patients with hospice transition: %s", format(n_lst_hospice, big.mark = ",")))

compute_timing_stats = function(x, label) {
  x = x[!is.na(x) & is.finite(x)]
  if (length(x) < 5) return(NULL)
  data.table(
    measure = label,
    n       = length(x),
    median  = round(median(x), 1),
    q25     = round(quantile(x, 0.25), 1),
    q75     = round(quantile(x, 0.75), 1),
    mean    = round(mean(x), 1),
    sd      = round(sd(x), 1)
  )
}

# by LST type
lst_types_for_timing = c("IMV only", "Vaso only", "IMV+Vaso", "Any LST")

# for "Any LST" use all lst_hospice rows
lst_timing_results = rbindlist(lapply(lst_types_for_timing, \(lt) {
  sub = if (lt == "Any LST") lst_hospice else lst_hospice[lst_type == lt]
  if (nrow(sub) < 5) return(NULL)

  rows = rbindlist(Filter(Negate(is.null), list(
    compute_timing_stats(sub$hours_first_withdrawal_to_hospice,
                         "Hours: LST withdrawal -> hospice"),
    compute_timing_stats(sub$hours_transition_to_death,
                         "Hours: hospice -> death"),
    compute_timing_stats(sub$hours_first_withdrawal_to_death,
                         "Hours: LST withdrawal -> death"),
    compute_timing_stats(sub$hours_comfort_to_hospice,
                         "Hours: comfort order -> hospice"),
    compute_timing_stats(sub$hours_comfort_to_withdrawal,
                         "Hours: comfort order -> LST withdrawal")
  )))

  if (nrow(rows) == 0) return(NULL)
  rows[, `:=`(lst_type = lt, site = site_lowercase)]
  rows
}))

if (nrow(lst_timing_results) > 0) {
  setcolorder(lst_timing_results, c("site", "lst_type", "measure", "n",
                                     "median", "q25", "q75", "mean", "sd"))

  message("\n  Key timing intervals (Any LST):")
  any_t = lst_timing_results[lst_type == "Any LST"]
  for (i in seq_len(nrow(any_t))) {
    message(sprintf("    %-45s median = %6.1f h (IQR %6.1f–%6.1f, n = %s)",
                    any_t$measure[i], any_t$median[i],
                    any_t$q25[i], any_t$q75[i],
                    format(any_t$n[i], big.mark = ",")))
  }

  fwrite(
    lst_timing_results,
    shareable_path(
              paste0("aim3_lst_timing_", site_lowercase, "_", today, ".csv"))
  )
  message("  Saved aim3_lst_timing.")
}


## 4b-viz: Timing box plot -----------------------------------------------------

if (nrow(lst_timing_results) > 0 && n_lst_hospice >= 20) {

  # box plot: hours from withdrawal to hospice, by LST type
  box_dt = lst_hospice[!is.na(hours_first_withdrawal_to_hospice) &
                        is.finite(hours_first_withdrawal_to_hospice) &
                        hours_first_withdrawal_to_hospice >= 0 &
                        hours_first_withdrawal_to_hospice <= 168]  # cap at 7d

  if (nrow(box_dt) >= 20) {
    p_box = ggplot(box_dt, aes(x = lst_type, y = hours_first_withdrawal_to_hospice)) +
      geom_boxplot(fill = "#D6E8F7", color = pal_outcome[["Same-facility hospice"]], outlier.size = 0.8,
                   outlier.alpha = 0.3, linewidth = 0.4) +
      coord_flip() +
      labs(
        x     = NULL,
        y     = "Hours from LST withdrawal to hospice transition",
        title = "LST withdrawal to hospice transition interval by therapy type",
        subtitle = paste0(toupper(site_lowercase), " | Capped at 7 days")
      ) +
      theme_clif() +
      theme(
      )

    ggsave(
      shareable_path(
                paste0("aim3_lst_timing_boxplot_", site_lowercase, "_", today, ".pdf")),
      plot = p_box, width = 9, height = 6, dpi = 300
    )
    message("  Saved aim3_lst_timing_boxplot.")

    # companion data: boxplot summary stats (not individual patient data)
    box_summary = box_dt[, .(
      n       = .N,
      median  = round(median(hours_first_withdrawal_to_hospice), 1),
      q25     = round(quantile(hours_first_withdrawal_to_hospice, 0.25), 1),
      q75     = round(quantile(hours_first_withdrawal_to_hospice, 0.75), 1),
      whisker_low  = round(quantile(hours_first_withdrawal_to_hospice, 0.05), 1),
      whisker_high = round(quantile(hours_first_withdrawal_to_hospice, 0.95), 1)
    ), by = lst_type]
    box_summary[, site := site_lowercase]
    fwrite(box_summary, shareable_path(
      paste0("aim3_lst_timing_boxplot_data_", site_lowercase, "_", today, ".csv")))
  }
}


## 4c. Code status sequencing --------------------------------------------------

message("\n  Code status sequencing relative to LST withdrawal...")

cs_seq = lst_hospice[, .N, by = code_status_sequence]
cs_seq[, pct := round(100 * N / sum(N), 2)]
cs_seq[, `:=`(site = site_lowercase, n_lst_hospice = n_lst_hospice)]

fwrite(
  cs_seq,
  shareable_path(
            paste0("aim3_code_status_sequence_", site_lowercase, "_", today, ".csv"))
)

message("  Code status sequence (LST hospice patients):")
for (i in seq_len(nrow(cs_seq))) {
  message(sprintf("    %-35s %s (%4.1f%%)",
                  cs_seq$code_status_sequence[i],
                  format(cs_seq$N[i], big.mark = ","),
                  cs_seq$pct[i]))
}


## 4c-viz: Event sequence stacked bar ------------------------------------------

if (nrow(cs_seq) > 0 && sum(cs_seq$N) >= 20) {

  cs_seq_plot = cs_seq[code_status_sequence != "No comfort order"]
  cs_seq_plot[, code_status_sequence := factor(code_status_sequence)]

  p_seq = ggplot(cs_seq_plot, aes(x = "", y = N, fill = code_status_sequence)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Comfort first"                    = pal_okabe_ito[5],
      "Withdrawal first"                 = pal_okabe_ito[6],
      "Simultaneous (<2h)"               = pal_okabe_ito[1],
      "Comfort only (no LST withdrawal)" = pal_okabe_ito[3]
    )) +
    labs(
      x     = NULL,
      y     = "Number of patients",
      title = "Sequencing of comfort care orders and LST withdrawal",
      subtitle = paste0(toupper(site_lowercase),
                        " | LST patients with hospice transition (n = ",
                        format(n_lst_hospice, big.mark = ","), ")"),
      fill  = NULL
    ) +
    theme_clif() +
    theme(
      axis.text.y      = element_blank()
    )

  ggsave(
    shareable_path(
              paste0("aim3_event_sequence_", site_lowercase, "_", today, ".pdf")),
    plot = p_seq, width = 9, height = 4, dpi = 300
  )
  message("  Saved aim3_event_sequence.")
}


# 5. Analysis 3.3: Hospital-level variation in LST-to-hospice ------------------

message("\n--- Analysis 3.3: Hospital-level LST-to-hospice variation ---")

# among LST patients who died (acute or hospice)
lst_died = lst[died == 1L]
n_lst_died = nrow(lst_died)

if (n_lst_died >= 20) {

  hosp_lst_rate = lst_died[, .(
    n_lst_died       = .N,
    n_hospice_died   = sum(any_hospice == 1L),
    hospice_rate_pct = round(100 * mean(any_hospice == 1L), 2)
  ), by = hospital_unit]

  hosp_lst_rate[, `:=`(site = site_lowercase)]

  # Spearman correlation with overall hospital hospice rate
  # (need overall hospice rate — compute from full cohort)
  overall_hosp_rate = cohort[, .(
    overall_hospice_rate = round(100 * mean(any_hospice == 1L), 2)
  ), by = hospital_unit]

  hosp_lst_rate = merge(hosp_lst_rate, overall_hosp_rate,
                        by = "hospital_unit", all.x = TRUE)

  if (nrow(hosp_lst_rate) >= 5) {
    rho = cor(hosp_lst_rate$hospice_rate_pct,
              hosp_lst_rate$overall_hospice_rate,
              method = "spearman", use = "complete.obs")
    message(sprintf("  Spearman: LST hospice rate vs overall hospice rate: rho = %.3f (n = %d)",
                    rho, nrow(hosp_lst_rate)))
  }

  fwrite(
    hosp_lst_rate,
    shareable_path(
              paste0("aim3_hospital_lst_hospice_rate_", site_lowercase, "_", today, ".csv"))
  )
  message("  Saved aim3_hospital_lst_hospice_rate.")


  ## Caterpillar plot
  if (nrow(hosp_lst_rate) >= 5) {

    hosp_lst_rate[, rank := rank(hospice_rate_pct)]
    overall_median = median(hosp_lst_rate$hospice_rate_pct)

    p_cat = ggplot(hosp_lst_rate, aes(x = rank, y = hospice_rate_pct)) +
      geom_point(size = 2.5, color = pal_outcome[["Same-facility hospice"]]) +
      geom_hline(yintercept = overall_median, linetype = "dashed",
                 color = "grey50", linewidth = 0.5) +
      annotate("text", x = max(hosp_lst_rate$rank), y = overall_median + 2,
               label = sprintf("Median: %.1f%%", overall_median),
               hjust = 1, size = 3, color = "grey40") +
      labs(
        x     = "Hospital rank",
        y     = "Hospice transition rate among LST deaths (%)",
        title = "Hospital-level variation: LST-to-hospice transition rate",
        subtitle = paste0(toupper(site_lowercase), " | ",
                          nrow(hosp_lst_rate), " hospitals")
      ) +
      theme_clif() +
      theme(
      )

    ggsave(
      shareable_path(
                paste0("aim3_lst_caterpillar_", site_lowercase, "_", today, ".pdf")),
      plot = p_cat, width = 9, height = 6, dpi = 300
    )
    message("  Saved aim3_lst_caterpillar.")
  }

  rm(overall_hosp_rate)

} else {
  message("  Insufficient LST deaths for hospital-level analysis.")
  hosp_lst_rate = NULL
}


# 6. Analysis 3.4: LST vs non-LST comparison -----------------------------------

message("\n--- Analysis 3.4: LST vs non-LST comparison ---")

# same-facility hospice transitions, split by LST status
sf_hospice = cohort[outcome_category %in% c("same_facility_hospice_death",
                                              "same_facility_hospice_alive")]

sf_hospice[, lst_group := fifelse(lst_any_flag == 1L, "LST", "Non-LST")]

n_lst_h = sum(sf_hospice$lst_group == "LST")
n_nonlst_h = sum(sf_hospice$lst_group == "Non-LST")
message(sprintf("  Hospice transitions: LST = %s, Non-LST = %s",
                format(n_lst_h, big.mark = ","),
                format(n_nonlst_h, big.mark = ",")))


## 6a. KM curves ---------------------------------------------------------------

if (n_lst_h >= 10 && n_nonlst_h >= 10) {

  sf_hospice[, ttd_hours := fifelse(
    !is.na(hours_transition_to_death) & hours_transition_to_death >= 0,
    pmin(hours_transition_to_death, 720), 720
  )]
  sf_hospice[, event := fifelse(
    !is.na(hours_transition_to_death) &
      hours_transition_to_death >= 0 &
      hours_transition_to_death <= 720,
    1L, 0L
  )]

  km_fit = survfit(Surv(ttd_hours, event) ~ lst_group, data = sf_hospice)

  # log-rank test
  lr_test = survdiff(Surv(ttd_hours, event) ~ lst_group, data = sf_hospice)
  lr_p = 1 - pchisq(lr_test$chisq, df = 1)

  message(sprintf("  Log-rank p-value: %.4g", lr_p))

  # extract KM tables
  km_summary = summary(km_fit, times = c(24, 48, 72, 168, 336))

  km_data = data.table(
    site        = site_lowercase,
    lst_group   = rep(c("LST", "Non-LST"), each = length(km_summary$time) / 2),
    time        = km_summary$time,
    n_risk      = km_summary$n.risk,
    n_event     = km_summary$n.event,
    surv        = round(km_summary$surv, 4),
    lower       = round(km_summary$lower, 4),
    upper       = round(km_summary$upper, 4)
  )

  # summary stats by group
  km_stats = sf_hospice[event == 1L, .(
    n              = .N,
    median_hours   = round(median(hours_transition_to_death, na.rm = TRUE), 1),
    q25_hours      = round(quantile(hours_transition_to_death, 0.25, na.rm = TRUE), 1),
    q75_hours      = round(quantile(hours_transition_to_death, 0.75, na.rm = TRUE), 1),
    pct_died_48h   = round(100 * mean(hours_transition_to_death <= 48, na.rm = TRUE), 2),
    pct_died_7d    = round(100 * mean(hours_transition_to_death <= 168, na.rm = TRUE), 2)
  ), by = lst_group]
  km_stats[, `:=`(site = site_lowercase, logrank_p = lr_p)]

  km_export = rbind(
    km_data[, .(site, lst_group, time, n_risk, n_event, surv, lower, upper)],
    fill = TRUE
  )

  fwrite(km_export,
         shareable_path(
                   paste0("aim3_lst_vs_nonlst_km_", site_lowercase, "_", today, ".csv")))
  message("  Saved aim3_lst_vs_nonlst_km.")

  # plot
  km_plot_data = data.table(
    time = c(0, km_fit$time),
    surv = c(1, 1, km_fit$surv[-1]),
    group = c("LST", "Non-LST",
              rep(c("LST", "Non-LST"), each = (length(km_fit$time) - 1) / 2 ))
  )

  # simpler approach: plot from survfit directly
  p_km34 = ggplot() +
    theme_minimal(base_size = 12)

  # extract per-stratum data
  strata_names = names(km_fit$strata)
  strata_n = km_fit$strata

  idx = 1
  km_lines = list()
  for (s in seq_along(strata_names)) {
    end_idx = idx + strata_n[s] - 1
    km_lines[[s]] = data.table(
      time_days = c(0, km_fit$time[idx:end_idx]) / 24,
      surv      = c(1, km_fit$surv[idx:end_idx]),
      group     = gsub("lst_group=", "", strata_names[s])
    )
    idx = end_idx + 1
  }
  km_lines = rbindlist(km_lines)

  p_km34 = ggplot(km_lines, aes(x = time_days, y = surv, color = group)) +
    geom_step(linewidth = 0.8) +
    scale_color_manual(values = c("LST" = pal_lst[["LST"]], "Non-LST" = pal_lst[["Non-LST"]])) +
    scale_x_continuous(breaks = c(0, 2, 7, 14, 21, 30), limits = c(0, 30)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    annotate("text", x = 25, y = 0.15,
             label = sprintf("Log-rank p = %.3g", lr_p),
             size = 3.5, color = "grey30") +
    labs(
      x     = "Days from hospice transition",
      y     = "Survival probability",
      title = "Time to death: LST vs non-LST hospice transitions",
      subtitle = paste0(toupper(site_lowercase),
                        " | LST n = ", format(n_lst_h, big.mark = ","),
                        ", Non-LST n = ", format(n_nonlst_h, big.mark = ",")),
      color = NULL
    ) +
    theme_clif() +
    theme(
    )

  ggsave(
    shareable_path(
              paste0("aim3_lst_vs_nonlst_km_", site_lowercase, "_", today, ".pdf")),
    plot = p_km34, width = 10, height = 7, dpi = 300
  )
  message("  Saved aim3_lst_vs_nonlst_km plot.")
}


# 6b. Clinical trajectory comparison -------------------------------------------

message("  Building LST vs non-LST trajectory comparison...")

traj_compare_components = c(
  "flag_code_status_strict", "flag_opioid", "flag_lab_cessation",
  "flag_withdrawal", "inflection_2of4", "inflection_1of4"
)
traj_compare_components = intersect(traj_compare_components, names(sf_hospice))

traj_compare = rbindlist(lapply(c("LST", "Non-LST"), \(grp) {
  sub = sf_hospice[lst_group == grp]
  n_sub = nrow(sub)
  if (n_sub < 5) return(NULL)

  comp = rbindlist(lapply(traj_compare_components, \(f) {
    data.table(
      variable = f, n = sum(sub[[f]] == 1L, na.rm = TRUE),
      pct = round(100 * mean(sub[[f]] == 1L, na.rm = TRUE), 2)
    )
  }))

  # add timing variables
  timing = rbindlist(Filter(Negate(is.null), list(
    if (sum(!is.na(sub$hours_comfort_to_hospice)) >= 5)
      data.table(variable = "hours_comfort_to_hospice (median)",
                 n = round(median(sub$hours_comfort_to_hospice, na.rm = TRUE), 1),
                 pct = NA_real_),
    if (sum(!is.na(sub$hours_transition_to_death)) >= 5)
      data.table(variable = "hours_transition_to_death (median)",
                 n = round(median(sub$hours_transition_to_death, na.rm = TRUE), 1),
                 pct = NA_real_),
    data.table(variable = "short_interval_death",
               n = sum(sub$short_interval_death == 1L, na.rm = TRUE),
               pct = round(100 * mean(sub$short_interval_death == 1L, na.rm = TRUE), 2))
  )))

  # lab trajectory
  lab_traj = data.table(
    variable = c("labs_72_48h (mean)", "labs_48_24h (mean)", "labs_24_0h (mean)"),
    n = c(round(mean(sub$labs_72_48h), 1),
          round(mean(sub$labs_48_24h), 1),
          round(mean(sub$labs_24_0h), 1)),
    pct = NA_real_
  )

  # opioid trajectory
  opi_traj = data.table(
    variable = c("opioid_72_48h (mean)", "opioid_48_24h (mean)", "opioid_24_0h (mean)"),
    n = c(round(mean(sub$opioid_72_48h), 1),
          round(mean(sub$opioid_48_24h), 1),
          round(mean(sub$opioid_24_0h), 1)),
    pct = NA_real_
  )

  out = rbind(comp, timing, lab_traj, opi_traj, fill = TRUE)
  out[, `:=`(lst_group = grp, site = site_lowercase, n_group = n_sub)]
  out
}))

if (nrow(traj_compare) > 0) {
  setcolorder(traj_compare, c("site", "lst_group", "n_group", "variable", "n", "pct"))

  fwrite(
    traj_compare,
    shareable_path(
              paste0("aim3_lst_vs_nonlst_trajectories_", site_lowercase, "_", today, ".csv"))
  )
  message("  Saved aim3_lst_vs_nonlst_trajectories.")
}


# 6b. Death-linked site sensitivity: external hospice trajectories

# At umn and ohsu, compare clinical trajectories for same-facility vs external hospice deaths

is_death_linked = "death_linked_site" %in% names(cohort) &&
                  any(cohort$death_linked_site == 1L, na.rm = TRUE)

if (is_death_linked) {

  message("\n--- Death-linked sensitivity: External hospice trajectories ---")

  # identify external hospice deaths
  cohort[, external_hospice_death := fifelse(
    outcome_category == "external_hospice_discharge" &
      !is.na(days_discharge_to_death) &
      days_discharge_to_death <= 30,
    1L, 0L
  )]

  n_ext_dead = sum(cohort$external_hospice_death == 1L)
  n_sf_dead  = sum(cohort$outcome_category == "same_facility_hospice_death")

  message(sprintf("  Same-facility hospice deaths: %s | External hospice deaths: %s",
                  format(n_sf_dead, big.mark = ","),
                  format(n_ext_dead, big.mark = ",")))

  if (n_ext_dead >= 10 && n_sf_dead >= 10) {

    # compare dying heuristic prevalence
    compare_cols = c(
      "flag_code_status_strict", "flag_opioid", "flag_lab_cessation",
      "flag_withdrawal", "inflection_2of4", "inflection_1of4"
    )
    compare_cols = intersect(compare_cols, names(cohort))

    sf_dead  = cohort[outcome_category == "same_facility_hospice_death"]
    ext_dead = cohort[external_hospice_death == 1L]

    pathway_compare = rbindlist(lapply(
      list(
        list(dt = sf_dead,  group = "Same-facility hospice death"),
        list(dt = ext_dead, group = "External hospice death")
      ), \(x) {
        n = nrow(x$dt)
        comp = rbindlist(lapply(compare_cols, \(f) {
          data.table(
            variable = f,
            n = sum(x$dt[[f]] == 1L, na.rm = TRUE),
            pct = round(100 * mean(x$dt[[f]] == 1L, na.rm = TRUE), 2)
          )
        }))

        # lab and opioid trajectories
        lab_rows = data.table(
          variable = c("labs_72_48h (mean)", "labs_48_24h (mean)", "labs_24_0h (mean)"),
          n = c(round(mean(x$dt$labs_72_48h), 1),
                round(mean(x$dt$labs_48_24h), 1),
                round(mean(x$dt$labs_24_0h), 1)),
          pct = NA_real_
        )

        opi_rows = data.table(
          variable = c("opioid_72_48h (mean)", "opioid_48_24h (mean)", "opioid_24_0h (mean)"),
          n = c(round(mean(x$dt$opioid_72_48h), 1),
                round(mean(x$dt$opioid_48_24h), 1),
                round(mean(x$dt$opioid_24_0h), 1)),
          pct = NA_real_
        )

        # time-to-death
        ttd_row = NULL
        if (x$group == "Same-facility hospice death" &&
            sum(!is.na(x$dt$hours_transition_to_death)) >= 5) {
          ttd_row = data.table(
            variable = "hours_transition_to_death (median)",
            n = round(median(x$dt$hours_transition_to_death, na.rm = TRUE), 1),
            pct = NA_real_
          )
        } else if (x$group == "External hospice death" &&
                   sum(!is.na(x$dt$days_discharge_to_death)) >= 5) {
          ttd_row = data.table(
            variable = "hours_discharge_to_death (median)",
            n = round(median(x$dt$days_discharge_to_death * 24, na.rm = TRUE), 1),
            pct = NA_real_
          )
        }

        out = rbind(comp, lab_rows, opi_rows, ttd_row, fill = TRUE)
        out[, `:=`(pathway = x$group, site = site_lowercase, n_group = n)]
        out
      }
    ))

    setcolorder(pathway_compare,
                c("site", "pathway", "n_group", "variable", "n", "pct"))

    # print summary -- core signals separately; composite as diagnostic only
    message("\n  Documented signals by pathway (reported separately):")
    for (pw in unique(pathway_compare$pathway)) {
      sub = pathway_compare[pathway == pw]
      n_g = sub$n_group[1]
      pv = function(v) { r = sub[variable == v]; if (nrow(r)) r$pct[1] else "N/A" }
      message(sprintf("    %s (n = %s): code-status %s%% | opioid %s%% | lab cess %s%% | LST wdl %s%%  [dx >=2/4: %s%%]",
                      pw, format(n_g, big.mark = ","),
                      pv("flag_code_status_strict"), pv("flag_opioid"),
                      pv("flag_lab_cessation"), pv("flag_withdrawal"), pv("inflection_2of4")))
    }

    fwrite(
      pathway_compare,
      shareable_path(
                paste0("aim3_sensitivity_pathway_comparison_", site_lowercase, "_", today, ".csv"))
    )
    message("  Saved aim3_sensitivity_pathway_comparison.")

    rm(sf_dead, ext_dead, pathway_compare)

  } else {
    message("  Insufficient external hospice deaths for comparison — skipping.")
  }
} else {
  message("\n--- Death-linked sensitivity: Skipped (not a death-linked site) ---")
}


# 7. cleanup -------------------------------------------------------------------


rm(cohort, cohort_all, lst, lst_hospice, lst_died, sf_hospice,
   lst_disposition, lst_any_disp, lst_timing_results,
   cs_seq, hosp_lst_rate, trajectory_summary, traj_compare)

suppressWarnings(rm(
  p_labs, p_opioid, p_box, p_seq, p_cat, p_km34,
  km_fit, km_data, km_stats, km_export, km_lines, km_summary,
  lr_test, lr_p, box_dt, cs_seq_plot,
  n_cohort, n_lst, n_hospice, n_lst_hospice, n_lst_died,
  n_lst_h, n_nonlst_h, flag_cols, traj_components, traj_compare_components,
  is_death_linked, n_ext_dead, n_sf_dead
))

gc()

message("\n===== 05_aim3_clinical_trajectories.R complete. =====\n")
