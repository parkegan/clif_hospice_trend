# 05_aim1.R
# Aim 1: Characterize hospice discharge patterns and timing
# Protocol v5.0 Section 3.6.1
#
#   1.1  Disposition distribution (overall, by age group)
#   1.2  Time-to-death distribution (KM curves, median/IQR, proportions at cutpoints)
#   1.3  Comparison to NHPCO national benchmarks
#   1.4  Clinical characteristics comparison (Table 1)
#
# Outputs (shareable CSVs in final_output/ for cross-site compilation):
#   - aim1_disposition_summary_{site}_{date}.csv
#   - aim1_time_to_death_{site}_{date}.csv
#   - aim1_table1_{site}_{date}.csv
#   - aim1_km_data_{site}_{date}.csv
#
# Outputs (local PDFs in final_output/):
#   - aim1_km_curve_{site}_{date}.pdf
#   - aim1_ttd_histogram_{site}_{date}.pdf
#
# Last updated bpe 16.07.2026


# 0. setup ---------------------------------------------------------------------

if (!exists("data_list")) source(here::here("code", "00_setup.R"))

if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")

library(ggplot2)
library(survival)
library(patchwork)
library(scales)


# 1. load data -----------------------------------------------------------------

message("\n===== 05_aim1.R: Aim 1 — Hospice Discharge Patterns & Timing =====\n")

# cohort with outcome_category and flags
cohort_all = fread(latest_output("cohort_all_hospitalizations"),
                   colClasses = list(character = c("patient_id", "hospitalization_id")))
cohort_all[, hospitalization_id := as.character(hospitalization_id)]
cohort_all[, patient_id := as.character(patient_id)]

# primary cohort
cohort = cohort_all[in_cohort == 1L]

# dying heuristic flags
dying = fread(latest_output("clinical_inflection_flags"))
dying[, hospitalization_id := as.character(hospitalization_id)]

# patient demographics
patient = data_list[["patient"]] |>
  dplyr::select(patient_id, race_category, ethnicity_category, sex_category) |>
  dplyr::collect() |>
  as.data.table()
patient[, patient_id := as.character(patient_id)]

# merge demographics onto cohort
cohort = merge(cohort, patient, by = "patient_id", all.x = TRUE)

# merge dying heuristic (keep only composite flags for Table 1)
dying_slim = dying[, .(
  hospitalization_id,
  inflection_2of4, inflection_1of4,
  n_core_components,
  flag_code_status_strict, flag_opioid, flag_lab_cessation, flag_withdrawal
)]
cohort = merge(cohort, dying_slim, by = "hospitalization_id", all.x = TRUE)

# fill NA heuristic flags with 0 (patients not in dying heuristic dataset)
heur_cols = c("inflection_2of4", "inflection_1of4",
              "n_core_components",
              "flag_code_status_strict", "flag_opioid",
              "flag_lab_cessation", "flag_withdrawal")
for (col in heur_cols) {
  if (col %in% names(cohort)) set(cohort, which(is.na(cohort[[col]])), col, 0L)
}

# defensive LST flag check
for (col in c("imv_flag", "vaso_flag", "lst_any_flag")) {
  if (!col %in% names(cohort)) {
    message(sprintf("  NOTE: %s not in cohort — defaulting to 0.", col))
    set(cohort, j = col, value = 0L)
  }
}

n_cohort = nrow(cohort)
message(sprintf("Primary cohort loaded: %s encounters", format(n_cohort, big.mark = ",")))


# 2. derived variables ----------------------------------------------------------

# age group
cohort[, age_group := fcase(
  age_at_admission < 40, "18-39",
  age_at_admission < 65, "40-64",
  age_at_admission < 80, "65-79",
  default = "80+"
)]

# analysis groups for Table 1
# Compare: those who died (acute death vs hospice death) and
# among hospice, short-survival (<48h) vs standard
cohort[, table1_group := fcase(
  outcome_category == "in_hospital_death_acute",          "Acute in-hospital death",
  outcome_category == "same_facility_hospice_death" &
    short_interval_death == 1L,                           "Hospice death (short-survival)",
  outcome_category == "same_facility_hospice_death" &
    (short_interval_death == 0L | is.na(short_interval_death)),
                                                          "Hospice death (standard)",
  outcome_category == "same_facility_hospice_alive",      "Hospice alive at 30d",
  outcome_category == "external_hospice_discharge",       "External hospice discharge",
  outcome_category == "routine_discharge",                "Routine discharge",
  outcome_category == "transfer_acute_facility",          "Transfer",
  default = "Other/AMA"
)]

# hospice transition flag (same-facility or external)
cohort[, any_hospice := fifelse(
  outcome_category %in% c("same_facility_hospice_death",
                           "same_facility_hospice_alive",
                           "external_hospice_discharge"),
  1L, 0L
)]

# death flag (broad)
cohort[, died := fifelse(
  outcome_category %in% c("in_hospital_death_acute",
                           "same_facility_hospice_death",
                           "external_hospice_death",
                           "post_discharge_death_non_hospice"),
  1L, 0L
)]

rm(patient, dying, dying_slim)


# 3. subgroup definitions (reused throughout) -----------------------------------

subgroup_defs = list(
  list(label = "full_cohort", filter_expr = quote(TRUE)),
  list(label = "icu_only",    filter_expr = quote(icu_flag == 1L)),
  list(label = "lst_any",     filter_expr = quote(lst_any_flag == 1L)),
  list(label = "lst_imv",     filter_expr = quote(imv_flag == 1L)),
  list(label = "lst_vaso",    filter_expr = quote(vaso_flag == 1L))
)


# =============================================================================
# Analysis 1.1: Disposition distribution
# =============================================================================

message("\n--- Analysis 1.1: Disposition distribution ---")

disposition_summary = rbindlist(lapply(subgroup_defs, \(sg) {
  sub = cohort[eval(sg$filter_expr)]
  n_sub = nrow(sub)
  if (n_sub == 0) return(NULL)

  # overall
  overall = sub[, .N, by = outcome_category]
  overall[, `:=`(stratification = "overall", stratum = "all")]

  # by age group
  by_age = sub[, .N, by = .(outcome_category, age_group)]
  setnames(by_age, "age_group", "stratum")
  by_age[, stratification := "age_group"]

  out = rbind(overall, by_age, fill = TRUE)
  out[, `:=`(
    subgroup   = sg$label,
    site       = site_lowercase,
    n_subgroup = n_sub,
    pct        = round(100 * N / n_sub, 2)
  )]
  out
}))

setcolorder(disposition_summary,
            c("site", "subgroup", "n_subgroup", "stratification", "stratum",
              "outcome_category", "N", "pct"))

fwrite(
  disposition_summary,
  shareable_path(
            paste0("aim1_disposition_summary_", site_lowercase, "_", today, ".csv"))
)
message("  Saved aim1_disposition_summary.")


# =============================================================================
# Analysis 1.2 & 1.3: Time-to-death distribution
# =============================================================================

message("\n--- Analysis 1.2/1.3: Time-to-death distribution ---")

# KM analysis for patients with a hospice transition (same-facility).
# Time: hours from hospice transition to death.
# Censoring: patients alive at 30 days (720 hours) post-transition.
# External hospice patients excluded from primary KM (no reliable transition time).

compute_cif_stats = function(dt, label) {
  # Aim 1.2 (protocol v7.0): competing-risks cumulative incidence, NOT Kaplan-
  # Meier. For same-facility hospice transitions, death in hospice competes with
  # live discharge -- a patient discharged alive from the unit can no longer die
  # in it, so live discharge is a competing event, not censoring. KM (1 - S)
  # OVERSTATES the probability of death; Aalen-Johansen is the correct estimator.
  # NOTE: "transfer to another hospice setting" (the protocol's third competing
  # event) is not separately coded in the current cohort, so it is folded into
  # live discharge. hospice_out_dttm supplies the live-discharge time; given the
  # negative-transition-time history, check its quality before trusting that CIF.
  hospice = dt[outcome_category %in% c("same_facility_hospice_death",
                                       "same_facility_hospice_alive")]
  hospice = hospice[!is.na(hospice_transition_dttm)]
  n_hospice = nrow(hospice)
  if (n_hospice < 5) {
    message(sprintf("    %s: only %d hospice transitions -- skipping", label, n_hospice))
    return(NULL)
  }

  # competing event times from hospice transition (time zero), horizon 720h
  hospice[, t_death := fifelse(!is.na(hours_transition_to_death) &
                                 hours_transition_to_death >= 0 &
                                 hours_transition_to_death <= 720,
                               hours_transition_to_death, Inf)]
  hospice[, t_live := {
    tl = as.numeric(difftime(hospice_out_dttm, hospice_transition_dttm, units = "hours"))
    fifelse(!is.na(tl) & tl > 0 & tl <= 720, tl, Inf)
  }]
  # earliest observed event wins; if both Inf the patient is censored at 720h
  hospice[, cr_event := fcase(
    is.finite(t_death) & t_death <= t_live,  "death",
    is.finite(t_live)  & t_live  <  t_death, "live_discharge",
    default = "censor")]
  hospice[, cr_time := fcase(
    cr_event == "death",          t_death,
    cr_event == "live_discharge", t_live,
    default = 720)]
  hospice[, cr_event := factor(cr_event, levels = c("censor", "death", "live_discharge"))]

  # Aalen-Johansen cumulative incidence via the multistate survfit
  aj = survfit(Surv(cr_time, cr_event) ~ 1, data = hospice)
  states    = aj$states
  i_death   = match("death", states)
  i_live    = match("live_discharge", states)
  i_still   = match("(s0)", states)

  cutpoints_h     = c(24, 48, 72, 168, 336)
  cutpoint_labels = c("24h", "48h", "72h", "7d", "14d")
  ss  = summary(aj, times = cutpoints_h)
  pst = ss$pstate
  cif_at_cut = data.table(
    window       = cutpoint_labels,
    hours        = cutpoints_h,
    cif_death    = round(100 * pst[, i_death], 2),
    cif_live     = round(100 * pst[, i_live], 2),
    still_in_pct = round(100 * pst[, i_still], 2)
  )

  full = summary(aj)
  cif_table = data.table(
    time      = full$time,
    cif_death = round(full$pstate[, i_death], 4),
    cif_live  = round(full$pstate[, i_live], 4),
    still_in  = round(full$pstate[, i_still], 4),
    subgroup  = label, site = site_lowercase
  )

  n_death = hospice[cr_event == "death", .N]
  n_live  = hospice[cr_event == "live_discharge", .N]
  n_cens  = hospice[cr_event == "censor", .N]

  # median time to death only exists if the death CIF reaches 0.5 (often it does
  # not under competing risks); report NA otherwise
  d_reach   = cif_table[cif_death >= 0.5, time]
  med_death = if (length(d_reach)) min(d_reach) else NA_real_

  summary_row = data.table(
    subgroup = label, site = site_lowercase,
    n_hospice = n_hospice, n_death = n_death,
    n_live_discharge = n_live, n_censored = n_cens,
    median_hours_to_death = round(med_death, 1),
    cif_death_24h = cif_at_cut[window == "24h", cif_death],
    cif_death_48h = cif_at_cut[window == "48h", cif_death],
    cif_death_72h = cif_at_cut[window == "72h", cif_death],
    cif_death_7d  = cif_at_cut[window == "7d",  cif_death],
    cif_death_14d = cif_at_cut[window == "14d", cif_death]
  )

  list(summary = summary_row, cif_table = cif_table, aj_fit = aj,
       hospice_dt = hospice, cif_at_cut = cif_at_cut)
}

# run for each subgroup
ttd_results = lapply(subgroup_defs, \(sg) {
  sub = cohort[eval(sg$filter_expr)]
  compute_cif_stats(sub, sg$label)
})
names(ttd_results) = sapply(subgroup_defs, \(sg) sg$label)
ttd_results = Filter(Negate(is.null), ttd_results)

ttd_summary  = rbindlist(lapply(ttd_results, \(r) r$summary))
cif_data_all = rbindlist(lapply(ttd_results, \(r) r$cif_table))

message("\n===== CUMULATIVE INCIDENCE OF DEATH IN HOSPICE (Aalen-Johansen) =====")
for (i in seq_len(nrow(ttd_summary))) {
  row = ttd_summary[i]
  message(sprintf("\n  %s (n = %s; %s deaths, %s live discharges, %s censored):",
                  row$subgroup, format(row$n_hospice, big.mark = ","),
                  row$n_death, row$n_live_discharge, row$n_censored))
  message(sprintf("    Death CIF  24h: %5.1f%% | 48h: %5.1f%% | 72h: %5.1f%% | 7d: %5.1f%% | 14d: %5.1f%%",
                  row$cif_death_24h, row$cif_death_48h, row$cif_death_72h,
                  row$cif_death_7d, row$cif_death_14d))
}
message("====================================================================\n")

fwrite(ttd_summary,  shareable_path(
  paste0("aim1_cif_summary_", site_lowercase, "_", today, ".csv")))
fwrite(cif_data_all, shareable_path(
  paste0("aim1_cif_data_", site_lowercase, "_", today, ".csv")))
message("  Saved aim1_cif_summary and aim1_cif_data.")

# Gray's test (protocol v7.0 replaces the KM log-rank under competing risks):
# equality of the death CIF across pandemic eras, on the full-cohort transitions.
if (requireNamespace("cmprsk", quietly = TRUE) &&
    !is.null(ttd_results[["full_cohort"]])) {
  g = ttd_results[["full_cohort"]]$hospice_dt
  if ("pandemic_period" %in% names(g) && uniqueN(g$pandemic_period) >= 2 && nrow(g) >= 20) {
    fstatus = as.integer(g$cr_event) - 1L   # 0 = censor, 1 = death, 2 = live discharge
    ci = tryCatch(cmprsk::cuminc(ftime = g$cr_time, fstatus = fstatus,
                                 group = as.character(g$pandemic_period), cencode = 0),
                  error = \(e) { message("  Gray's test failed: ", e$message); NULL })
    if (!is.null(ci) && !is.null(ci$Tests)) {
      gt = as.data.table(ci$Tests, keep.rownames = "cause")
      gt[, `:=`(site = site_lowercase, comparison = "pandemic_period")]
      # cause "1" = death in hospice
      pv = suppressWarnings(as.numeric(gt[cause == "1", 3][1]))
      message(sprintf("  Gray's test, death CIF across pandemic eras: p = %s",
                      if (is.na(pv)) "NA" else formatC(pv, format = "f", digits = 4)))
      fwrite(gt, shareable_path(
        paste0("aim1_cif_grays_test_", site_lowercase, "_", today, ".csv")))
    }
  } else {
    message("  Gray's test skipped (need >= 2 eras and >= 20 transitions).")
  }
}


# =============================================================================
# Analysis 1.2/1.3 Visualizations: KM curve + histogram
# =============================================================================

message("\n--- Generating KM curves ---")

# KM curve for full cohort (primary visualization)
if ("full_cohort" %in% names(ttd_results)) {

  aj_fc      = ttd_results[["full_cohort"]]$aj_fit
  hospice_fc = ttd_results[["full_cohort"]]$hospice_dt
  cif_fc     = ttd_results[["full_cohort"]]$cif_table

  # cumulative incidence of DEATH in hospice (rising curve). Under competing
  # risks this is the correct quantity; 1 - KM would sit above it.
  cif_plot_dt = rbind(
    data.table(time = 0, cif = 0, cause = "Death in hospice"),
    cif_fc[, .(time, cif = cif_death, cause = "Death in hospice")],
    data.table(time = 0, cif = 0, cause = "Live discharge"),
    cif_fc[, .(time, cif = cif_live,  cause = "Live discharge")]
  )
  cif_plot_dt[, time_days := time / 24]

  p_km = ggplot(cif_plot_dt, aes(x = time_days, y = cif, color = cause)) +
    geom_step(linewidth = 0.8) +
    geom_vline(xintercept = c(2, 7, 14), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    # NHPCO benchmark: ~28% die within 7 days -> death CIF ~0.28 at 7d
    annotate("point", x = 7, y = 0.28, shape = 4, size = 3,
             color = pal_outcome[["In-hospital death"]], stroke = 1.5) +
    annotate("text", x = 8, y = 0.28, label = "NHPCO ref:\n28% died by 7d",
             hjust = 0, size = 3, color = pal_outcome[["In-hospital death"]]) +
    scale_x_continuous(breaks = c(0, 2, 7, 14, 21, 30), limits = c(0, 30)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_color_manual(values = c("Death in hospice" = pal_outcome[["In-hospital death"]],
                                  "Live discharge"   = pal_outcome[["Routine discharge"]])) +
    labs(
      x = "Days from hospice transition",
      y = "Cumulative incidence",
      color = NULL,
      title = "Cumulative incidence after same-facility hospice transition (Aalen-Johansen)",
      subtitle = paste0(toupper(site_lowercase),
                        " | n = ", format(nrow(hospice_fc), big.mark = ","),
                        " | competing risks: death in hospice vs live discharge")
    ) +
    theme_clif()

  # LST subgroup overlay: death CIF only, dashed
  lst_cif = ttd_results[["lst_any"]]
  if (!is.null(lst_cif) && lst_cif$summary$n_hospice >= 10) {
    lst_dt = rbind(data.table(time = 0, cif_death = 0),
                   lst_cif$cif_table[, .(time, cif_death)])
    lst_dt[, time_days := time / 24]
    p_km = p_km +
      geom_step(data = lst_dt, aes(x = time_days, y = cif_death),
                inherit.aes = FALSE, linewidth = 0.7,
                color = pal_lst[["LST"]], linetype = "dashed") +
      annotate("text", x = 28, y = max(lst_dt$cif_death[lst_dt$time_days <= 28], 0.1),
               label = sprintf("LST death CIF\n(n=%s)",
                               format(lst_cif$summary$n_hospice, big.mark = ",")),
               hjust = 1, size = 3, color = pal_lst[["LST"]])
  }

  ggsave(
    shareable_path(
              paste0("aim1_cif_curve_", site_lowercase, "_", today, ".pdf")),
    plot = p_km, width = 10, height = 7, dpi = 300
  )
  message("  Saved aim1_cif_curve (Aalen-Johansen cumulative incidence).")


  # histogram of time-to-death (Analysis 1.2/1.3)
  hist_dt = hospice_fc[cr_event == "death" & cr_time > 0][, ttd_hours := cr_time]

  if (nrow(hist_dt) >= 10) {
    hist_dt[, ttd_days := ttd_hours / 24]

    p_hist = ggplot(hist_dt, aes(x = ttd_days)) +
      geom_histogram(binwidth = 1, fill = pal_outcome[["Same-facility hospice"]], alpha = 0.7,
                     color = "white", linewidth = 0.2) +
      geom_vline(xintercept = 2, linetype = "dashed", color = pal_outcome[["In-hospital death"]],
                 linewidth = 0.5) +
      geom_vline(xintercept = 7, linetype = "dashed", color = pal_outcome[["In-hospital death"]],
                 linewidth = 0.5) +
      # NHPCO national median ~18 days
      geom_vline(xintercept = 18, linetype = "dotted", color = pal_outcome[["Routine discharge"]],
                 linewidth = 0.7) +
      annotate("text", x = 2.5, y = Inf, label = "48h", vjust = 1.5,
               size = 3, color = pal_outcome[["In-hospital death"]]) +
      annotate("text", x = 7.5, y = Inf, label = "7d", vjust = 1.5,
               size = 3, color = pal_outcome[["In-hospital death"]]) +
      annotate("text", x = 18.5, y = Inf,
               label = "NHPCO median\n(~18 days)", vjust = 1.5,
               size = 3, color = pal_outcome[["Routine discharge"]]) +
      scale_x_continuous(breaks = c(0, 2, 7, 14, 18, 21, 28, 30)) +
      labs(
        x     = "Days from hospice transition to death",
        y     = "Number of patients",
        title = "Distribution of Time to Death After Same-Facility Hospice Transition",
        subtitle = paste0(toupper(site_lowercase),
                          " | n = ", format(nrow(hist_dt), big.mark = ","),
                          " deaths | Bins = 1 day")
      ) +
      theme_clif() +
      theme(
      )

    ggsave(
      shareable_path(
                paste0("aim1_ttd_histogram_", site_lowercase, "_", today, ".pdf")),
      plot = p_hist, width = 10, height = 6, dpi = 300
    )

    # companion data: binned counts for cross-site re-plotting
    hist_bins = hist_dt[, .(n = .N), by = .(day = floor(ttd_days))]
    hist_bins[, `:=`(site = site_lowercase)]
    fwrite(hist_bins, shareable_path(
      paste0("aim1_ttd_histogram_data_", site_lowercase, "_", today, ".csv")))

    message("  Saved aim1_ttd_histogram + data.")
  }
}


# =============================================================================
# Analysis 1.4: Table 1 — clinical characteristics by group
# =============================================================================

message("\n--- Analysis 1.4: Table 1 (shareable summary) ---")

# Build aggregated summary statistics by table1_group and subgroup.
# This produces a stacked long-format CSV: one row per variable × level × group.
# No patient-level data is shared.

build_table1 = function(dt, subgroup_label) {

  n_dt = nrow(dt)
  if (n_dt == 0) return(NULL)

  groups = unique(dt$table1_group)
  groups = groups[!is.na(groups)]

  result_rows = list()

  for (grp in groups) {
    g = dt[table1_group == grp]
    n_g = nrow(g)
    if (n_g == 0) next

    add_row = function(variable, level, n, pct = round(100 * n / n_g, 2),
                       mean_val = NA_real_, sd_val = NA_real_,
                       median_val = NA_real_, q25_val = NA_real_,
                       q75_val = NA_real_) {
      result_rows[[length(result_rows) + 1L]] <<- data.table(
        site        = site_lowercase,
        subgroup    = subgroup_label,
        table1_group = grp,
        n_group     = n_g,
        variable    = variable,
        level       = level,
        n           = n,
        pct         = pct,
        mean_val    = mean_val,
        sd_val      = sd_val,
        median_val  = median_val,
        q25_val     = q25_val,
        q75_val     = q75_val
      )
    }

    # --- N ---
    add_row("N", "total", n_g, pct = 100)

    # --- Age (continuous) ---
    add_row("Age", "continuous", n_g,
            mean_val   = round(mean(g$age_at_admission, na.rm = TRUE), 1),
            sd_val     = round(sd(g$age_at_admission, na.rm = TRUE), 1),
            median_val = round(median(g$age_at_admission, na.rm = TRUE), 1),
            q25_val    = round(quantile(g$age_at_admission, 0.25, na.rm = TRUE), 1),
            q75_val    = round(quantile(g$age_at_admission, 0.75, na.rm = TRUE), 1))

    # --- Age group ---
    for (ag in c("18-39", "40-64", "65-79", "80+")) {
      add_row("Age group", ag, sum(g$age_group == ag, na.rm = TRUE))
    }

    # --- Sex ---
    for (s in unique(na.omit(g$sex_category))) {
      add_row("Sex", s, sum(g$sex_category == s, na.rm = TRUE))
    }

    # --- Race ---
    for (r in unique(na.omit(g$race_category))) {
      add_row("Race", r, sum(g$race_category == r, na.rm = TRUE))
    }

    # --- Ethnicity ---
    for (e in unique(na.omit(g$ethnicity_category))) {
      add_row("Ethnicity", e, sum(g$ethnicity_category == e, na.rm = TRUE))
    }

    # --- LOS (continuous) ---
    add_row("LOS (days)", "continuous", n_g,
            mean_val   = round(mean(g$los_days, na.rm = TRUE), 1),
            sd_val     = round(sd(g$los_days, na.rm = TRUE), 1),
            median_val = round(median(g$los_days, na.rm = TRUE), 1),
            q25_val    = round(quantile(g$los_days, 0.25, na.rm = TRUE), 1),
            q75_val    = round(quantile(g$los_days, 0.75, na.rm = TRUE), 1))

    # --- ICU admission ---
    add_row("ICU admission", "yes", sum(g$icu_flag == 1L, na.rm = TRUE))

    # --- ICU within 48h ---
    add_row("ICU within 48h", "yes", sum(g$icu_48h_flag == 1L, na.rm = TRUE))

    # --- ED admission ---
    add_row("ED admission", "yes", sum(g$ed_flag == 1L, na.rm = TRUE))

    # --- Age >= 65 ---
    add_row("Age >= 65", "yes", sum(g$age_65_flag == 1L, na.rm = TRUE))

    # --- LST flags ---
    add_row("Any LST", "yes", sum(g$lst_any_flag == 1L, na.rm = TRUE))
    add_row("IMV", "yes", sum(g$imv_flag == 1L, na.rm = TRUE))
    add_row("Vasopressors", "yes", sum(g$vaso_flag == 1L, na.rm = TRUE))

    # --- Clinical inflection ---
    add_row("Clinical inflection (>=2/4)", "yes",
            sum(g$inflection_2of4 == 1L, na.rm = TRUE))
    add_row("Clinical inflection (>=1/4)", "yes",
            sum(g$inflection_1of4 == 1L, na.rm = TRUE))

    # --- Individual heuristic components ---
    add_row("Code status change", "yes",
            sum(g$flag_code_status_strict == 1L, na.rm = TRUE))
    add_row("Opioid initiation/escalation", "yes",
            sum(g$flag_opioid == 1L, na.rm = TRUE))
    add_row("Lab cessation", "yes",
            sum(g$flag_lab_cessation == 1L, na.rm = TRUE))
    add_row("LST withdrawal", "yes",
            sum(g$flag_withdrawal == 1L, na.rm = TRUE))

    # --- Time to death (for hospice deaths only) ---
    hosp_died = g[!is.na(hours_transition_to_death) & hours_transition_to_death >= 0]
    if (nrow(hosp_died) >= 5) {
      add_row("Hours: transition to death", "continuous", nrow(hosp_died),
              mean_val   = round(mean(hosp_died$hours_transition_to_death, na.rm = TRUE), 1),
              sd_val     = round(sd(hosp_died$hours_transition_to_death, na.rm = TRUE), 1),
              median_val = round(median(hosp_died$hours_transition_to_death, na.rm = TRUE), 1),
              q25_val    = round(quantile(hosp_died$hours_transition_to_death, 0.25, na.rm = TRUE), 1),
              q75_val    = round(quantile(hosp_died$hours_transition_to_death, 0.75, na.rm = TRUE), 1))
    }

    # --- Short survival (<48h) among hospice transitions ---
    if (any(g$any_hospice == 1L)) {
      add_row("Short survival (<48h)", "yes",
              sum(g$short_interval_death == 1L, na.rm = TRUE))
    }
  }

  rbindlist(result_rows)
}

# run for each subgroup
table1_all = rbindlist(lapply(subgroup_defs, \(sg) {
  sub = cohort[eval(sg$filter_expr)]
  build_table1(sub, sg$label)
}))

# print summary for full cohort
message("\n========== TABLE 1 SUMMARY (full cohort, key groups) ==========")
fc_t1 = table1_all[subgroup == "full_cohort"]
for (grp in c("Acute in-hospital death", "Hospice death (short-survival)",
              "Hospice death (standard)")) {
  gt = fc_t1[table1_group == grp]
  if (nrow(gt) == 0) next
  n_g = gt[variable == "N", n]
  message(sprintf("\n  %s (n = %s):", grp, format(n_g, big.mark = ",")))
  age_row = gt[variable == "Age" & level == "continuous"]
  if (nrow(age_row)) {
    message(sprintf("    Age: mean %.1f (SD %.1f), median %.1f",
                    age_row$mean_val, age_row$sd_val, age_row$median_val))
  }
  icu_row = gt[variable == "ICU admission"]
  if (nrow(icu_row)) {
    message(sprintf("    ICU admission: %s (%4.1f%%)",
                    format(icu_row$n, big.mark = ","), icu_row$pct))
  }
  lst_row = gt[variable == "Any LST"]
  if (nrow(lst_row)) {
    message(sprintf("    Any LST: %s (%4.1f%%)",
                    format(lst_row$n, big.mark = ","), lst_row$pct))
  }
  rd_row = gt[variable == "Recognized dying (>=2/4)"]
  if (nrow(rd_row)) {
    message(sprintf("    Recognized dying (>=2/4): %s (%4.1f%%)",
                    format(rd_row$n, big.mark = ","), rd_row$pct))
  }
}
message("================================================================\n")

# save
fwrite(
  table1_all,
  shareable_path(
            paste0("aim1_table1_", site_lowercase, "_", today, ".csv"))
)
message("  Saved aim1_table1.")


# =============================================================================
# Analysis 1.5: Sensitivity analysis at death-linked sites
# =============================================================================

# At death-linked sites (OHSU, UMN), the patient death_dttm captures deaths
# from state registries, including those occurring after external hospice
# discharge. This allows us to include external hospice deaths alongside
# same-facility hospice deaths and compare the two pathways.

is_death_linked = "death_linked_site" %in% names(cohort) &&
                  any(cohort$death_linked_site == 1L, na.rm = TRUE)

if (is_death_linked) {

  message("\n--- Analysis 1.5: Death-linked site sensitivity ---")

  # identify external hospice discharges with confirmed death within 30 days
  cohort[, external_hospice_death := fifelse(
    outcome_category == "external_hospice_discharge" &
      !is.na(days_discharge_to_death) &
      days_discharge_to_death <= 30,
    1L, 0L
  )]

  n_ext_dead = sum(cohort$external_hospice_death == 1L)
  n_ext_total = sum(cohort$outcome_category == "external_hospice_discharge")
  message(sprintf("  External hospice discharges: %s total, %s died within 30 days (%4.1f%%)",
                  format(n_ext_total, big.mark = ","),
                  format(n_ext_dead, big.mark = ","),
                  100 * n_ext_dead / max(n_ext_total, 1)))

  # combined hospice deaths: same-facility + external with confirmed death
  cohort[, hospice_death_combined := fifelse(
    outcome_category == "same_facility_hospice_death" |
      external_hospice_death == 1L,
    1L, 0L
  )]

  # --- 1.5a: Disposition with external hospice deaths reclassified ---

  disp_sens = cohort[, .N, by = outcome_category]
  disp_sens[, `:=`(site = site_lowercase, analysis = "primary")]

  # create a sensitivity version where external hospice deaths are separated
  cohort_sens = copy(cohort)
  cohort_sens[external_hospice_death == 1L,
              outcome_category := "external_hospice_death"]

  disp_sens2 = cohort_sens[, .N, by = outcome_category]
  disp_sens2[, `:=`(site = site_lowercase, analysis = "sensitivity_death_linked")]

  disp_sensitivity = rbind(disp_sens, disp_sens2)
  disp_sensitivity[, pct := round(100 * N / sum(N), 2), by = analysis]

  # --- 1.5b: Time-to-death comparison: same-facility vs external hospice ---

  # same-facility hospice deaths
  sf_deaths = cohort[outcome_category == "same_facility_hospice_death" &
                       !is.na(hours_transition_to_death)]

  # external hospice deaths: time = days_discharge_to_death * 24
  ext_deaths = cohort[external_hospice_death == 1L &
                        !is.na(days_discharge_to_death)]
  ext_deaths[, hours_discharge_to_death := days_discharge_to_death * 24]

  ttd_comparison = rbindlist(list(
    if (nrow(sf_deaths) >= 5) {
      data.table(
        pathway        = "Same-facility hospice",
        n              = nrow(sf_deaths),
        median_hours   = round(median(sf_deaths$hours_transition_to_death), 1),
        q25_hours      = round(quantile(sf_deaths$hours_transition_to_death, 0.25), 1),
        q75_hours      = round(quantile(sf_deaths$hours_transition_to_death, 0.75), 1),
        pct_died_48h   = round(100 * mean(sf_deaths$hours_transition_to_death <= 48), 2),
        pct_died_7d    = round(100 * mean(sf_deaths$hours_transition_to_death <= 168), 2)
      )
    },
    if (nrow(ext_deaths) >= 5) {
      data.table(
        pathway        = "External hospice",
        n              = nrow(ext_deaths),
        median_hours   = round(median(ext_deaths$hours_discharge_to_death), 1),
        q25_hours      = round(quantile(ext_deaths$hours_discharge_to_death, 0.25), 1),
        q75_hours      = round(quantile(ext_deaths$hours_discharge_to_death, 0.75), 1),
        pct_died_48h   = round(100 * mean(ext_deaths$hours_discharge_to_death <= 48), 2),
        pct_died_7d    = round(100 * mean(ext_deaths$hours_discharge_to_death <= 168), 2)
      )
    }
  ))

  if (nrow(ttd_comparison) > 0) {
    ttd_comparison[, site := site_lowercase]

    message("\n  Time-to-death comparison:")
    for (i in seq_len(nrow(ttd_comparison))) {
      row = ttd_comparison[i]
      message(sprintf("    %s (n = %s): median = %.1f hours | 48h: %4.1f%% | 7d: %4.1f%%",
                      row$pathway, format(row$n, big.mark = ","),
                      row$median_hours, row$pct_died_48h, row$pct_died_7d))
    }
  }

  # --- Save sensitivity outputs ---

  fwrite(
    disp_sensitivity,
    shareable_path(
              paste0("aim1_sensitivity_disposition_", site_lowercase, "_", today, ".csv"))
  )

  if (nrow(ttd_comparison) > 0) {
    fwrite(
      ttd_comparison,
      shareable_path(
                paste0("aim1_sensitivity_ttd_comparison_", site_lowercase, "_", today, ".csv"))
    )
  }

  message("  Saved aim1_sensitivity outputs.")

  rm(disp_sens, disp_sens2, disp_sensitivity, cohort_sens,
     sf_deaths, ext_deaths, ttd_comparison)

} else {
  message("\n--- Analysis 1.5: Skipped (not a death-linked site) ---")
}


# =============================================================================
# Analysis E.1: Temporal trends
# =============================================================================

message("\n--- Analysis E.1: Temporal trends ---")

# quarterly outcome counts, hospice rates, mortality rates, and ratios

cohort[, discharge_quarter := floor_date(discharge_dttm, "quarter")]

qtr = cohort[, .(
  n_total       = .N,
  n_death_acute = sum(outcome_category == "in_hospital_death_acute"),
  n_hospice_sf  = sum(outcome_category %in% c("same_facility_hospice_death",
                                                "same_facility_hospice_alive")),
  n_hospice_ext = sum(outcome_category == "external_hospice_discharge"),
  n_hospice_all = sum(outcome_category %in% c("same_facility_hospice_death",
                                                "same_facility_hospice_alive",
                                                "external_hospice_discharge")),
  n_transfer    = sum(outcome_category == "transfer_acute_facility"),
  n_short_surv  = sum(short_interval_death == 1L, na.rm = TRUE)
), by = discharge_quarter]

qtr[, `:=`(
  rate_death_acute = round(1000 * n_death_acute / n_total, 2),
  rate_hospice_sf  = round(1000 * n_hospice_sf / n_total, 2),
  rate_hospice_all = round(1000 * n_hospice_all / n_total, 2),
  rate_transfer    = round(1000 * n_transfer / n_total, 2),
  hospice_mortality_ratio = round(n_hospice_sf / pmax(n_death_acute, 1), 3),
  pct_short_survival = round(100 * n_short_surv / pmax(n_hospice_sf, 1), 1)
)]

qtr[, site := site_lowercase]
setorder(qtr, discharge_quarter)

# 4-quarter rolling means
if (nrow(qtr) >= 4) {
  qtr[, roll_rate_hospice := frollmean(rate_hospice_sf, n = 4, align = "right")]
  qtr[, roll_rate_death   := frollmean(rate_death_acute, n = 4, align = "right")]
  qtr[, roll_ratio        := frollmean(hospice_mortality_ratio, n = 4, align = "right")]
}

message(sprintf("  Quarterly data: %d quarters (%s to %s)",
                nrow(qtr),
                format(min(qtr$discharge_quarter), "%Y Q%q"),
                format(max(qtr$discharge_quarter), "%Y Q%q")))

# temporal trend plots
if (nrow(qtr) >= 4) {

  # stacked rates over time
  trend_long = melt(qtr,
                    id.vars = "discharge_quarter",
                    measure.vars = c("rate_hospice_sf", "rate_death_acute", "rate_transfer"),
                    variable.name = "outcome", value.name = "rate_per_1000")
  trend_long[, outcome := fcase(
    outcome == "rate_hospice_sf",  "Same-facility hospice",
    outcome == "rate_death_acute", "In-hospital death",
    outcome == "rate_transfer",    "Transfer"
  )]

  p_trends = ggplot(trend_long, aes(x = discharge_quarter, y = rate_per_1000,
                                     color = outcome)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    # pandemic shading
    annotate("rect",
             xmin = as.POSIXct("2020-03-01"), xmax = as.POSIXct("2023-05-11"),
             ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey50") +
    scale_color_manual(values = c(
      "Same-facility hospice" = pal_outcome[["Same-facility hospice"]],
      "In-hospital death"     = pal_outcome[["In-hospital death"]],
      "Transfer"              = pal_outcome[["Transfer"]]
    )) +
    labs(x = NULL, y = "Rate per 1,000 admissions",
         title = "Quarterly disposition rates",
         subtitle = paste0(toupper(site_lowercase), " | Shaded = pandemic period"),
         color = NULL) +
    theme_clif()

  # hospice-to-mortality ratio over time
  p_ratio = ggplot(qtr, aes(x = discharge_quarter, y = hospice_mortality_ratio)) +
    geom_line(linewidth = 0.7, color = pal_outcome[["Same-facility hospice"]]) +
    geom_point(size = 1.5, color = pal_outcome[["Same-facility hospice"]]) +
    {if ("roll_ratio" %in% names(qtr))
      geom_line(aes(y = roll_ratio), linewidth = 0.9, color = pal_lst[["LST"]], linetype = "solid")} +
    annotate("rect",
             xmin = as.POSIXct("2020-03-01"), xmax = as.POSIXct("2023-05-11"),
             ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey50") +
    geom_hline(yintercept = 0.27, linetype = "dashed", color = pal_outcome[["In-hospital death"]], linewidth = 0.4) +
    annotate("text", x = min(qtr$discharge_quarter), y = 0.27,
             label = "National median (0.27)", hjust = 0, vjust = -0.5,
             size = 3, color = pal_outcome[["In-hospital death"]]) +
    labs(x = NULL, y = "Hospice-to-mortality ratio",
         title = "Same-facility hospice to in-hospital mortality ratio",
         subtitle = paste0(toupper(site_lowercase),
                           " | Amber = 4-quarter rolling mean"),
         color = NULL) +
    theme_clif()

  p_combined = p_trends / p_ratio + plot_layout(heights = c(1, 1))

  ggsave(shareable_path(
                   paste0("aim1_temporal_trends_", site_lowercase, "_", today, ".pdf")),
         plot = p_combined, width = 10, height = 10, dpi = 300)
  message("  Saved aim1_temporal_trends.")

  suppressWarnings(rm(trend_long, p_trends, p_ratio, p_combined))
}

# save quarterly data for cross-site compilation
fwrite(qtr,
       shareable_path(
                 paste0("aim1_quarterly_trends_", site_lowercase, "_", today, ".csv")))
message("  Saved aim1_quarterly_trends.")


# =============================================================================
# Analysis E.2: Patients admitted from hospice
# =============================================================================

message("\n--- Analysis E.2: Patients admitted from hospice ---")

# Load the admitted-from-hospice cohort saved by 01_cohort.R
hospice_admit_file = list.files(
  private_dir,
  pattern = "cohort_admitted_from_hospice", full.names = TRUE
)

if (length(hospice_admit_file) > 0) {

  hospice_admits = fread(hospice_admit_file[order(hospice_admit_file, decreasing = TRUE)][1])
  hospice_admits[, hospitalization_id := as.character(hospitalization_id)]
  n_ha = nrow(hospice_admits)
  message(sprintf("  Admitted-from-hospice cohort: %s encounters", format(n_ha, big.mark = ",")))

  if (n_ha >= 20) {

    # disposition distribution
    ha_disp = hospice_admits[, .N, by = outcome_category]
    ha_disp[, pct := round(100 * N / n_ha, 2)]
    ha_disp[, `:=`(site = site_lowercase, group = "admitted_from_hospice", n_total = n_ha)]

    # compare to newly transitioned hospice patients in primary cohort
    new_hospice = cohort[outcome_category %in% c("same_facility_hospice_death",
                                                   "same_facility_hospice_alive")]
    n_nh = nrow(new_hospice)

    # clinical characteristics comparison
    compare_groups = list(
      list(dt = hospice_admits, label = "Admitted from hospice"),
      list(dt = new_hospice,    label = "Newly transitioned to hospice")
    )

    e2_comparison = rbindlist(lapply(compare_groups, \(g) {
      dt = g$dt; n = nrow(dt)
      if (n < 10) return(NULL)

      rows = list(
        data.table(variable = "N", level = "total", value = n),
        data.table(variable = "Age (median)", level = "continuous",
                   value = round(median(dt$age_at_admission, na.rm = TRUE), 1)),
        data.table(variable = "LOS days (median)", level = "continuous",
                   value = round(median(dt$los_days, na.rm = TRUE), 1)),
        data.table(variable = "ICU admission (%)", level = "pct",
                   value = round(100 * mean(dt$icu_flag == 1L, na.rm = TRUE), 1)),
        data.table(variable = "Any LST (%)", level = "pct",
                   value = round(100 * mean(dt$lst_any_flag == 1L, na.rm = TRUE), 1)),
        data.table(variable = "In-hospital death (%)", level = "pct",
                   value = round(100 * mean(dt$outcome_category == "in_hospital_death_acute",
                                             na.rm = TRUE), 1))
      )

      # time-to-death for those who died
      if ("hours_transition_to_death" %in% names(dt)) {
        died = dt[!is.na(hours_transition_to_death) & hours_transition_to_death >= 0]
        if (nrow(died) >= 5) {
          rows = c(rows, list(
            data.table(variable = "Hours to death (median, hospice deaths)", level = "continuous",
                       value = round(median(died$hours_transition_to_death), 1)),
            data.table(variable = "Died within 48h (%)", level = "pct",
                       value = round(100 * mean(died$hours_transition_to_death <= 48), 1))
          ))
        }
      }

      out = rbindlist(rows)
      out[, `:=`(group = g$label, site = site_lowercase, n_group = n)]
      out
    }))

    if (nrow(e2_comparison) > 0) {
      setcolorder(e2_comparison, c("site", "group", "n_group", "variable", "level", "value"))

      message("\n  Comparison:")
      for (grp in unique(e2_comparison$group)) {
        sub = e2_comparison[group == grp]
        message(sprintf("    %s (n = %s):", grp, format(sub$n_group[1], big.mark = ",")))
        for (i in seq_len(nrow(sub)))
          message(sprintf("      %-40s %s", sub$variable[i], sub$value[i]))
      }

      fwrite(e2_comparison,
             shareable_path(
                       paste0("aim1_e2_hospice_admit_comparison_", site_lowercase, "_", today, ".csv")))
      message("  Saved aim1_e2_hospice_admit_comparison.")
    }

    rm(ha_disp, new_hospice, e2_comparison, compare_groups)
  } else {
    message("  Too few admitted-from-hospice encounters for analysis.")
  }

  rm(hospice_admits)

} else {
  message("  No admitted-from-hospice cohort file found. Run 01_cohort.R first.")
}


# =============================================================================
# cleanup
# =============================================================================

rm(cohort, cohort_all, disposition_summary, ttd_results, ttd_summary,
   km_data_all, table1_all, subgroup_defs,
   heur_cols, qtr)

# remove plot objects if they exist
suppressWarnings(rm(p_km, p_hist, km_plot_dt, lst_km_dt, hist_dt, hospice_fc,
                    km_fc, lst_km, fc_t1,
                    is_death_linked))

gc()

message("\n===== 03_aim1.R complete. =====\n")
