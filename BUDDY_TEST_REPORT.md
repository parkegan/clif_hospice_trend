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
| **Date** | *2026-10-19* |

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
| 2 | Configuration works from `config/README.md` alone; no hardcoding | Pass | Need to add config update instructions to the readme. |
| 3 | Required tables/fields match what the code reads (mCIDE-valid) | Pass  | |
| 4 | Runs end to end with no manual edits between steps | Fail | Stage 04 crashes for any single-CCN site, and stage 05 never ran.|
| 5 | Outputs in `output/final_no_phi/` with right naming/type, no raw dumps | Pass  | Outputs go to `ucmc_upload_to_box/`|
| 6 | **Data security**: no PHI, every stat n ≥ 10, no raw data *(blocking)* | Pass with notes | Some output files have n < 5 |
| 7 | Clinical sanity: aggregates plausible for the cohort |  Fail | |
| 8 | Documentation usable: could run from the README alone | Fail | |

## Overall verdict

**Verdict:** Fail 

- **Pass** — runs cleanly, output sane and secure, docs followable. Ready to distribute.
- **Pass with notes** — works; address the notes below.
- **Fail** — at least one blocking issue (doesn't run, output wrong/insecure, docs unfollowable). A
  data-security failure (check 6) is always a Fail.

### Blocking issues (must fix before distribution)
1. Stage 04 crashes for any single-CCN site, and stage 05 never ran.

```
Hospitals with inpatient hospice unit: 0 of 1
  Variation summary (all hospitals):
Error in `[.data.table`(var_stats, version == "all_hospitals") : 
  i evaluates to a logical vector length 14 but there are 0 rows. Recycling of logical i is no longer allowed as it hides more bugs than is worth the rare convenience. Explicitly use rep(...,length=.N) if you really need to recycle.
Calls: nrow ... [ -> [.data.table -> stopf -> raise_condition -> signal
Execution halted
    -> FAILED (exit 1)  (37.6s)

==========================================================================
TIMING SUMMARY
--------------------------------------------------------------------------
  step stage                                               elapsed   status
  01   Cohort construction                                   26.9s   ok
  02   Clinical inflection & code status                     50.8s   ok
  03   Aim 1: discharge patterns & competing-risks CI         8.0s   ok
  04   Aim 2A: risk-standardized variation                   37.6s   FAILED
--------------------------------------------------------------------------
       TOTAL                                                2m 03s
==========================================================================
```

2. `aim1_disposition_summary`: N of 1, 2 and 5 in the 80+ age stratum, for example lst_imv / other_ama with N = 1. `aim1_quarterly_trends`: n_short_surv of 3–4, n_transfer of 1, and one quarter with n_total of 1.
3. `pct_short_survival` is 300–500% in aim1_quarterly_trends. 
4. There is a 2025 Q2 quarter with a single admission. Does the study window allow that?
5. `hospice_mortality_ratio` is 0 in the Aim 1 file and 0.88 in the Aim 2 file for the same quarter, so the column means different things in the two files.
6. Mortality (3.7%) and external hospice (1.9%) look plausible. Same-facility hospice is 0, which is correct: UCMC has no inpatient hospice unit (pipeline classified it b_unavailable).

### Non-blocking notes / suggestions
1. *…*

