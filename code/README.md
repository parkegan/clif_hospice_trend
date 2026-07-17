# CLIF Hospice Discharge Study — Analysis Pipeline (code reference)

Detailed per-script documentation for the federated CLIF hospice-discharge pipeline. For environment setup, running the pipeline, and version control, see the **top-level `README.md`**; this file documents what each script does.

**Protocol:** v7.0 (July 2026)
**Study period:** 2016-01-01 through 2024-06-30
**Cohort:** Adult (≥18) inpatient admissions across CLIF consortium sites (~40 hospitals)

Only de-identified, CCN-level aggregates leave a site. No patient-level data is transmitted.

---

## Pipeline overview

The pipeline is a linear sequence of numbered stages that communicate through files on disk — each stage reads the most recent output of the previous one, so stages can be run, inspected, and re-run independently. Run the whole thing with `Rscript code/run_pipeline.R` (see top-level README), which times each stage.

```
00_setup.R          (sourced by every stage; not run directly)
    │
    ▼
01_cohort.R                    cohort, outcomes, QI flags, hospice-pathway availability
    │
    ▼
02_clinical_inflection.R       code-status classification + inflection signals
    │
    ▼
03_aim1.R                      disposition patterns + competing-risks cumulative incidence
    │
    ▼
04_aim2a3_risk_adjusted_variation.R   risk-standardized hospital variation (Aim 2A)
    │
    ▼
05_aim3_clinical_trajectories.R       documented care-transition trajectories (Aim 3)

aim3_chart_validation.R        chart-validation substudy (run separately; manual review step)
```

**Death-linked sites:** OHSU and UMN carry state death-registry data in the CLIF `death_dttm` field. Stages 03–05 detect these sites and run additional sensitivity analyses that include external hospice deaths; other sites skip those blocks with no runtime cost.

**Federated design:** most of Aim 2A's hospital-variation quantities are not estimable at a single site. Each site emits per-hospital estimates standardized onto a common population; the coordinating center pools them across all ~40 hospitals (see the top-level README and `federation_estimand_design.md`).

---

## Script descriptions

### 00_setup.R

**Purpose:** Environment configuration, package loading, CLIF table access, and shared helpers. Sourced by every stage; never run directly.

- Loads packages (managed by renv) and detects cores/memory for parallel steps.
- Reads site config from `config/config.yaml` (site name, study dates, paths to the CLIF tables).
- Opens the CLIF parquet tables into `data_list` (named Arrow datasets) and validates required columns.
- Resolves the project root via `rprojroot::find_root(has_dir("config"))`.
- Loads reference files, each optional with a clear message when absent:
  - `reference/hospital_ccn_map.csv` — hospital_id → CMS CCN for CCN-level aggregation.
  - `reference/crosswalk_ccsr_to_acuity.csv` — CCSR → acuity diagnosis grouping (`load_ccsr_acuity_crosswalk()`).
  - `config/network_standard_population.csv` — the federated common standard population (`load_standard_population()`).
- Defines `save_output()` (routes to the private vs shareable directory), `shareable_path()`, `latest_output()` (loads the newest matching intermediate — the loader that keeps stages in sync), and `qi_flag()` (prints a data-quality banner without stopping the run).

**Outputs:** none (configures the session).

---

### 01_cohort.R

**Purpose:** Build the analytic cohort, classify outcomes, and emit site data-quality flags.

- Pulls adult inpatient admissions in the study window; tracks sequential exclusions with counts.
- **Earliest-timestamp death rule:** where `death_dttm` precedes `discharge_dttm` (an organ-donor admission recorded after death, or discharge paperwork logged after an in-hospital death), discharge is clamped to the death time, and any admission recorded entirely after death is dropped. This prevents negative death intervals rather than discarding the rows.
- Extracts ADT flags (ICU, ICU-within-48h, ED, hospice location, first location, hospital_id) and the LST subcohort (IMV, vasopressors, RRT, and the union `lst_any_flag`).
- Classifies each hospitalization into one of seven outcome categories (routine discharge, in-hospital death, same-facility hospice death, same-facility hospice alive, external hospice discharge, transfer, other/AMA), and derives `any_hospice_disposition` (same-facility **or** external), the network-wide outcome.
- **Hospice-pathway availability (§3.4.4):** classifies each hospital as `a_available_captured`, `b_unavailable`, or `c_uncertain` from its same-facility transition counts and yearly coverage. A pathway that opens or closes mid-study is forced to class (c). A site attestation file (`config/hospital_pathway_map.csv`) overrides the data-driven classification.
- **QI flags:** missing ADT `hospital_id` (with a by-admission-year and by-encounter-type breakdown to show whether the gap is random or systematic), residual negative death intervals after clamping, and uncertain hospice-pathway availability.

**Key outputs:** `cohort_all_hospitalizations_*` (private), `hospital_pathway_availability_*` (shareable), exclusion and subgroup summaries.

---

### 02_clinical_inflection.R

**Purpose:** Code-status classification and the documented clinical-inflection signals.

- Maps `code_status_category` to Full / Limited resuscitation / Comfort-focused / Other. Network decisions baked in: **AND ("Allow Natural Death") → Comfort-focused**, and **"Special/Partial" → Limited resuscitation** (captures the "LIMITED"-named orders, at the cost of a few VAD device rows the category can't separate).
- Windows code-status orders to each encounter and derives `first_comfort_dttm` and related timing.
- Builds the four core documented signals (code-status change, opioid initiation/escalation, lab cessation, LST withdrawal) and the composite `inflection_2of4`/`inflection_1of4`. **Note:** per v7.0 the composite is retained only as a diagnostic; the signals are reported separately downstream.

**Key outputs:** `clinical_inflection_flags_*` (private), signal-prevalence summaries.

---

### 03_aim1.R

**Purpose:** Aim 1 — hospice discharge patterns and timing.

- Disposition distribution overall and by age group, full cohort + LST subgroups.
- **Analysis 1.2 (competing risks, not Kaplan-Meier):** for same-facility hospice transitions, the **Aalen–Johansen cumulative incidence** of death in hospice versus live discharge, with the still-in-hospice state, reported at 24h/48h/72h/7d/14d. Live discharge is a competing event, so 1−KM would overstate death. **Transfer to another hospice setting** folds into live discharge where a site does not separately capture it (documented limitation).
- **Gray's test** replaces the log-rank for comparing the death CIF across pandemic eras.
- CIF curves plotted with the NHPCO benchmark reference; histogram of observed hospice duration.
- **Analysis 1.5 (death-linked sites only):** external hospice discharges with confirmed death ≤30 days; primary-vs-sensitivity disposition and timing comparisons.

**Key outputs:** `aim1_cif_summary_*`, `aim1_cif_data_*`, `aim1_cif_grays_test_*` (shareable), `aim1_cif_curve_*.pdf`.

---

### 04_aim2a3_risk_adjusted_variation.R

**Purpose:** Aim 2A.3 — risk-standardized hospital-level variation in discharge disposition. The most changed script under v7.0.

**Covariates (all fixed at or shortly after admission — the primary model excludes post-admission variables):**
- Age, sex, race/ethnicity.
- **`dx_category` = acuity grouping.** ICD-10 principal codes are mapped to CCSR via `function_ccsr_bodysystem.R`, then `ccsr_category` is joined to `reference/crosswalk_ccsr_to_acuity.csv` to produce the 12-level acuity grouping (cancer, cancer in remission, dementia/neurodegenerative/neurodevelopmental, chronic organ failure, frailty, acute infection/sepsis, …). Body system is retained as `dx_body_system`; if the crosswalk is absent, `dx_category` falls back to body system.
- Elixhauser comorbidity (van Walraven index).
- APACHE II severity from a **24-hour** early window (`function_apache_II.R`).
- **Preexisting code status** (the code status documented by admission +24h, including standing pre-admission orders) and **admission source** (schema-detected).

**Models:**
- **Primary:** per-disposition binary GLMMs (`lme4::glmer`, one-vs-routine-discharge) with admission-fixed covariates and hospital random intercepts. Dynamic random-effect structure adds a health-system term (where a site spans ≥3 systems) and a patient-level intercept (where estimable), with a prespecified **one-admission-per-patient** sensitivity as the fallback.
- **Secondary clinical-course model:** adds post-admission course (ICU, LST) back on — a **different estimand**, written to its own files.
- **Multinomial composition (Option C):** `mclogit::mblogit` supplies the risk-standardized disposition composition via marginal standardization onto a common population (probabilities sum to 1). The binary GLMMs supply the uncertainty suite, and the two are cross-checked for agreement.

**Variation profiling (replaces the ICC):** between-hospital variance with profile-likelihood CI, median odds ratio, risk-standardized probabilities with intervals, reliability, and rank uncertainty. ICC is retained only as a secondary column.

**Pathway restriction:** same-facility outcomes are modeled only on class (a) hospitals; `any_hospice_disposition` is modeled across all hospitals.

**Federation:** hospitals are standardized onto the distributed common population; each per-hospital estimate carries a standard error and a standardization-source/version tag. The site emits `aim2a3_standardization_margins_*` for assembling the pooled standard.

**Key outputs:** `aim2a3_model_coefficients_*`, `aim2a3_variation_summary_*`, `aim2a3_hospital_profile_*`, `aim2a3_rs_composition_*`, `aim2a3_model_consistency_*`, the `*_clinical_course_*` counterparts, and `aim2a3_standardization_margins_*` (all shareable).

---

### 05_aim3_clinical_trajectories.R

**Purpose:** Aim 3 (exploratory) — documented care-transition trajectories preceding hospice transition.

- Reports the four core documented signals **separately** (not as a composite). The `inflection_2of4`/`inflection_1of4` composites are carried only as diagnostics, tagged `signal_role = "diagnostic_composite"`.
- LST disposition flow, withdrawal-to-hospice and hospice-to-death timing by LST type, and code-status sequencing relative to withdrawal.
- Hospital-level LST-to-hospice transition rates among LST decedents; pathway-stratified comparisons.
- Death-linked sensitivity: trajectory comparison between same-facility and external hospice deaths.

**Key outputs:** `aim3_trajectory_summary_*`, `aim3_lst_*`, `aim3_code_status_sequence_*` (shareable).

---

### aim3_chart_validation.R

**Purpose:** Aim 3 chart-validation substudy (§3.6.3). Run separately — the middle step is a manual chart review.

- **Mode 1 (default):** draws a stratified random sample of algorithm-flagged hospice transitions and writes a blank review worksheet (kept private, since it holds encounter IDs).
- **Mode 2 (`AIM3_CV_SCORE=<worksheet.csv>`):** scores the positive predictive value of the hospice-transition algorithm and of the AND→comfort mapping, with Wilson confidence intervals. Refuses to run on a blank worksheet, so it cannot fabricate a gold standard.

---

### function_apache_II.R

**Purpose:** APACHE II severity scoring from CLIF vitals and labs. Sourced by `04`, not run standalone.

Extracts worst-scoring values for the acute physiology variables in the **first 24 hours** (v7.0 changed this from 48h; 24h is also APACHE II's canonical window), adds age and optional chronic-health points, and returns hospitalization-level scores with completeness flags.

---

### function_ccsr_bodysystem.R

**Purpose:** Map ICD-10 principal diagnosis codes to CCSR categories. Sourced by `04`, not run standalone. **Still required.**

- `load_ccsr_lookup()` loads the ICD-10 → CCSR lookup; `assign_ccsr_bodysystem()` returns both the 530-level `ccsr_category` and a 22-level body-system grouping.
- Its role changed under v7.0: `dx_category` is now the **acuity** grouping, but this function is load-bearing because `ccsr_category` is the **join key** for the acuity crosswalk — without it there is no CCSR code to map to acuity. It also provides the `dx_body_system` reference column and the fallback grouping used when the acuity crosswalk is not present.
- The name is now a mild misnomer (it does ICD-10→CCSR, not just body system). Renaming is optional cleanup; if you do, update the `list.files(pattern = ...)` lookup in `04`.

---

### run_pipeline.R / setup.R / _dependencies.R

Infrastructure, documented in the top-level README:
- `run_pipeline.R` — runs stages 01→05 as isolated processes with per-stage timing and a timing CSV; supports `--from`, `--only`, `--keep-going`, `--list`.
- `setup.R` — one-time renv bootstrap (restore from `renv.lock`, or initialise and snapshot).
- `_dependencies.R` — explicit package manifest for renv discovery; never sourced by the pipeline.

---

## Output routing

Each stage writes to two locations set in `config.yaml`:

- **Private** (`<site>_private_data/`) — intermediates that may contain identifiers; never leave the site.
- **Shareable** (`<site>_upload_to_box/`) — de-identified, CCN-level aggregates the coordinating center collects.

All shareable CSVs carry a `site` column for cross-site compilation. Watch the run log for **QI DATA-QUALITY FLAG** banners; several map directly onto the §3.7 site data-quality profile the protocol requires.

---

## Reference files

| File | Location | Shared or site-specific |
|---|---|---|
| `hospital_ccn_map.csv` | `reference/` | Site-specific (each site's hospitals) |
| `crosswalk_ccsr_to_acuity.csv` | `reference/` | **Shared** — all sites use the same version |
| `network_standard_population.csv` | `config/` | **Shared** — distributed by the coordinating center |
| `hospital_pathway_map.csv` | `config/` | Site-specific (attestation) |

Shared reference files must be identical across sites for cross-site estimates to be comparable; site-specific files are supplied locally (see `.gitignore`).

---

## Remaining gaps

| Item | Status |
|---|---|
| Federated estimand / pooling design | Designed (`federation_estimand_design.md`) and the site-side standardization is wired in; awaits the coordinating center finalizing the standard population |
| Comfort/AND windowing | Under investigation — few AND orders land inside encounter windows; likely pre-admission standing orders |
| ADT coverage gap (~30%) | Confirmed a genuine coverage issue, not a join bug; the by-year/by-type QI breakdown characterizes it per site |
| `SYM016` acuity mapping | Flagged for review but absent from the reviewed sheet; provisionally `Other/residual` |
| Early-window LST/organ-support for the clinical-course model | `icu_flag`/`lst_any_flag` are whole-stay; windowing to the early window is follow-up |
