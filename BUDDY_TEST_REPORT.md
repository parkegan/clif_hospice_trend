<!--
Buddy test report template. See guides/buddy-testing-guide.md for what each check means.
To use: copy this file to the project root as BUDDY_TEST_REPORT.md, fill it in, and commit it.
In the Result column, bold the outcome that applies (or delete the other), then add a note.
Pass with notes is only for the Overall verdict at the bottom.
-->

# Buddy Test Report

| | |
|---|---|
| **Buddy site / institution** | *UCMC* |
| **Tester** | *Kaveri Chhikara* |
| **Date** | *2026-10-19 (re-run 2026-10-20)* |

## Environment

| | |
|---|---|
| **OS** | *macOS 26.6.2* |
| **RAM** | *64 GB* |
| **R** | *4.5.2* |


## Checks

| # | Check | Result | Notes |
|:-:|-------|--------|-------|
| 1 | Environment reproduces (`uv sync` / `00_renv_restore.R`, nothing by hand) | Pass  | |
| 2 | Configuration works from `config/README.md` alone; no hardcoding | Pass | README now has config instructions. |
| 3 | Required tables/fields match what the code reads (mCIDE-valid) | Pass  | |
| 4 | Runs end to end with no manual edits between steps | Pass | Single-CCN crash in stage 04 fixed; all 5 stages run. |
| 5 | Outputs in `output/final_no_phi/` with right naming/type, no raw dumps | Pass  | Outputs go to `ucmc_upload_to_box/`|
| 6 | **Data security**: no PHI, every stat n ≥ 10, no raw data *(blocking)* | Pass with notes | Counts < 5 suppressed (n < 5 threshold accepted). |
| 7 | Clinical sanity: aggregates plausible for the cohort | Pass | |
| 8 | Documentation usable: could run from the README alone | Pass | See note 2. |

## Overall verdict

**Verdict:** Pass with notes


### Blocking issues (must fix before distribution)
None.

### Resolved (re-run 2026-10-20)
1. Stage 04 crashed for any single-CCN site, and stage 05 never ran. Fixed; all 5 stages run.
2. Counts < 5 in `aim1_disposition_summary` and `aim1_quarterly_trends`. Now suppressed.
3. `pct_short_survival` was 300–500%. Now NA at sites without same-facility hospice.
4. 2025 Q2 quarter with a single admission. Removed. 2025 Q1 (discharges of late-2024 admissions) stays in the output; the coordinating center drops quarters from 2025-01-01 on during aggregation.
5. `hospice_mortality_ratio` meant different things in the Aim 1 and Aim 2 files. Now `sf_hospice_to_mortality_ratio` and `all_hospice_to_mortality_ratio`.
6. Same-facility hospice is 0, which is correct: UCMC has no inpatient hospice unit (pipeline classified it b_unavailable).

### Non-blocking notes / suggestions

