# CLIF Hospice Discharge Study

A federated electronic-health-record analysis of hospice discharge patterns among hospitalized adults, run across the CLIF consortium using CLIF 2.1 data.

---

## Code Pipeline:

- **Aim 1 — discharge patterns and timing.** Where and when hospitalized adults are discharged to hospice, length of time between hospice discharge and death.
- **Aim 2A — risk-standardized disposition composition.** How hospitals vary in end-of-life discharge disposition after case-mix adjustment, reported as risk-standardized probabilities, between-hospital variance, the median odds ratio, reliability, and rank uncertainty.
- **Aim 3 — documented care-transition trajectories (exploratory).** Documented signals prior to hospice transition (code-status change, opioid initiation, lab cessation, life-sustaining-therapy withdrawal).

| Stage | Script | Produces |
|---|---|---|
| 00 | `00_setup.R` | Shared config, helpers, package loading (sourced by every stage; not run directly) |
| 01 | `01_cohort.R` | The analytic cohort, outcomes, and site data-quality (QI) flags |
| 02 | `02_clinical_inflection.R` | Code-status classification and clinical-inflection signals |
| 03 | `03_aim1.R` | Aim 1 disposition tables and the competing-risks cumulative incidence |
| 04 | `04_aim2a3_risk_adjusted_variation.R` | Aim 2A models and the risk-standardized variation profile |
| 05 | `05_aim3_clinical_trajectories.R` | Aim 3 documented-signal trajectories |
| — | `aim3_chart_validation.R` | Chart-validation substudy (will be completed at a later date) |

Utilities are `function_ccsr_bodysystem.R` and `function_apache_II.R`

---

## Configuration

1. Navigate to the `config/` directory.
2. Rename `config_template.yaml` to `config.yaml`.
3. Update the `config.yaml` with site-specific settings.

## Setup, by operating system

Do this once per machine, from the project root.

### macOS

From a terminal in the project root:
   ```bash
   Rscript code/setup.R
   ```

### Windows

From PowerShell or Command Prompt in the project root:
   ```powershell
   Rscript code/setup.R
   ```

### Linux (Debian/Ubuntu)

From the project root:
   ```bash
   Rscript code/setup.R
   ```

---

## Running the pipeline

From the project root:

```bash
Rscript code/run_pipeline.R
```

The runner executes stages 01 → 05 in order and prints a timing summary.

```
==========================================================================
TIMING SUMMARY
--------------------------------------------------------------------------
  step stage                                               elapsed   status
  01   Cohort construction                                 4m 12s    ok
  02   Clinical inflection & code status                   1m 58s    ok
  03   Aim 1: discharge patterns & competing-risks CIF       47s     ok
  04   Aim 2A: risk-standardized variation                 6m 03s    ok
  05   Aim 3: documented trajectories                      2m 20s    ok
--------------------------------------------------------------------------
       TOTAL                                              15m 20s
==========================================================================
```


Useful flags:

| Command | Effect |
|---|---|
| `Rscript code/run_pipeline.R --list` | List the stages and exit |
| `Rscript code/run_pipeline.R --from 03` | Resume from stage 03 onward |
| `Rscript code/run_pipeline.R --only 01,04` | Run only the named stages |
| `Rscript code/run_pipeline.R --keep-going` | Continue past a failing stage instead of stopping |

The runner exits non-zero if any stage fails, so it drops cleanly into a scheduler or CI job. If a stage fails, fix it and resume with `--from <stage>` rather than rerunning from the top.

---

## The chart-validation substudy - ignore for now (will have one or two sites do this later)

Aim 3's chart validation is not part of the automated run because the middle step is a manual chart review.

---

## Outputs

Each stage writes to two folders:

- a **private** directory (`<site>_private_data/`) for intermediates which never leave the site, and
- a **shareable** directory (`<site>_upload_to_box/`) for the de-identified, CCN-level aggregates.

---

## Reproducibility notes

- Package versions are pinned by `renv.lock`; do not `install.packages()` outside renv.
- Randomized steps (the standardization sample, rank-uncertainty draws, the chart-validation sample) are seeded.
