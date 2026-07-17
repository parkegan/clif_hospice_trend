# CLIF Hospice Discharge Study

A federated electronic-health-record analysis of hospice discharge patterns among hospitalized adults, run across the CLIF consortium (~40 hospitals, 2016–2024) using the CLIF common data model.

---

## Code Pipeline:

- **Aim 1 — discharge patterns and timing.** Where and when hospitalized adults are discharged to hospice, and, for same-facility hospice transitions, the competing-risks cumulative incidence of death in hospice versus live discharge (Aalen–Johansen, with Gray's test).
- **Aim 2A — risk-standardized disposition composition.** How hospitals vary in end-of-life discharge disposition after case-mix adjustment, reported as risk-standardized probabilities, between-hospital variance, the median odds ratio, reliability, and rank uncertainty. Composition comes from a mixed multinomial model; the variance and uncertainty suite comes from complementary binary mixed models.
- **Aim 3 — documented care-transition trajectories (exploratory).** The documented signals around a hospice transition (code-status change, opioid initiation, lab cessation, life-sustaining-therapy withdrawal), reported separately, plus a chart-validation substudy estimating the positive predictive value of the hospice-transition algorithm.

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

## Directory layout

The pipeline locates the project root by looking for a `config/` directory, so the layout matters:

```
<project_root>/
├── .here                         # anchor so here::here() resolves to the root (created by setup.R)
├── setup.R                       # one-time environment setup (renv)
├── renv.lock                     # locked package versions (created on first setup; commit it)
├── README.md
├── config/
│   ├── config.yaml               # REQUIRED: site name, paths to the CLIF tables, options
│   ├── network_standard_population.csv   # (federation) common standard population, if distributed
│   └── hospital_pathway_map.csv          # (optional) site attestation of hospice-pathway availability
├── reference/
│   ├── hospital_ccn_map.csv      # hospital_id -> CMS CCN, for CCN-level aggregation
│   └── crosswalk_ccsr_to_acuity.csv      # CCSR -> acuity diagnosis-category crosswalk
└── code/
    ├── 00_setup.R
    ├── 01_cohort.R … 05_aim3_clinical_trajectories.R
    ├── function_ccsr_bodysystem.R
    ├── function_apache_II.R
    ├── aim3_chart_validation.R
    ├── _dependencies.R           # package manifest for renv
    └── run_pipeline.R            # the pipeline runner
```

Place the finalized reference CSVs where the layout shows: the acuity crosswalk at `reference/crosswalk_ccsr_to_acuity.csv`, and — for the federated Aim 2A — the coordinating center's standard population at `config/network_standard_population.csv`. Both are shared reference files: every site must use the same version for cross-site estimates to be comparable.

---

## Prerequisites

- **R 4.2 or newer.**
- Roughly 16 GB of RAM is comfortable for a large site; the stages release memory between each other because each runs as its own process.
- Package management is handled by **renv**, so you do not install packages by hand — `setup.R` does it and pins the versions.

---

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

The same table is written to `pipeline_timing_<timestamp>.csv` in the project root, so you can track run time across sites and over time.

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

Watch the run log for **QI DATA-QUALITY FLAG** banners. They surface issues a site should resolve or attest to before its results are pooled — missing ADT `hospital_id`, uncertain hospice-pathway availability, residual data-quality problems — and several map directly onto the site data-quality profile the protocol requires.

---

## Reproducibility notes

- Package versions are pinned by `renv.lock`; do not `install.packages()` outside renv.
- Randomized steps (the standardization sample, rank-uncertainty draws, the chart-validation sample) are seeded, so a given input and package set reproduce the same numbers.
- Most of Aim 2A is not estimable at a single site — the hospital-variation quantities become estimable only when the coordinating center pools per-hospital estimates across the network. A single-site run standardizes locally and reports what it can; that is expected.
