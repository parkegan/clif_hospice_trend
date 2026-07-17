# =============================================================================
# _dependencies.R  —  explicit package manifest for renv.
# -----------------------------------------------------------------------------
# renv discovers dependencies by scanning the project for library()/:: calls.
# Listing every package here in one place guarantees a complete, reviewable
# lockfile even for packages referenced only inside tryCatch()/requireNamespace()
# (e.g. cmprsk, mclogit). This file is NEVER sourced by the pipeline; it exists
# solely so `renv::snapshot()` captures the full set.
#
# Grouped by role for readability.

# --- data I/O and manipulation ---
library(arrow)        # reads the CLIF parquet tables
library(data.table)   # primary data engine across all scripts
library(dplyr)        # arrow/dplyr verbs at the ingestion layer
library(collapse)     # fast grouped aggregation
library(fst)          # fast on-disk intermediate storage
library(stringr)      # string handling
library(lubridate)    # datetime handling

# --- project plumbing ---
library(here)         # path resolution relative to project root
library(rprojroot)    # project-root detection (also used by run_pipeline.R)
library(yaml)         # config.yaml
library(rlang)        # tidy evaluation helpers
library(ps)           # process / system introspection

# --- clinical measures ---
library(comorbidity)  # Elixhauser / van Walraven comorbidity scoring

# --- statistics ---
library(lme4)         # binary mixed models (Laplace) for Aim 2A profiling
library(mclogit)      # multinomial mixed model (mblogit) for disposition composition
library(survival)     # Aalen-Johansen competing-risks CIF (Aim 1.2)
library(cmprsk)       # Gray's test for competing risks (Aim 1.2)

# --- figures ---
library(ggplot2)      # all plots
library(scales)       # axis scales/formatting
library(patchwork)    # multi-panel composition

# --- parallelism ---
# 'parallel' ships with base R and needs no installation; listed for the record.
# library(parallel)
