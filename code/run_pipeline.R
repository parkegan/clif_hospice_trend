#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R  —  run the CLIF hospice pipeline end to end, with timing.
# -----------------------------------------------------------------------------
# Each stage is run as its own Rscript process (the stages communicate through
# files on disk, so process isolation is both correct and safer: a failure in
# one stage cannot corrupt the R session of the next, and memory is released
# between stages). Wall-clock time is recorded per stage and summarized at the
# end, and also written to a timing CSV.
#
# Usage (from the PROJECT ROOT, the folder that contains config/ and code/):
#   Rscript code/run_pipeline.R                # run every stage in order
#   Rscript code/run_pipeline.R --from 03      # resume from stage 03 onward
#   Rscript code/run_pipeline.R --only 01,04   # run only these stages
#   Rscript code/run_pipeline.R --keep-going   # continue past a failing stage
#   Rscript code/run_pipeline.R --list         # list stages and exit
# =============================================================================

if (!requireNamespace("rprojroot", quietly = TRUE))
  stop("Package 'rprojroot' is required. Run  Rscript setup.R  first.", call. = FALSE)

project_root <- rprojroot::find_root(rprojroot::has_dir("config"))
code_dir     <- file.path(project_root, "code")
rscript      <- file.path(R.home("bin"),
                          if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

# ---- pipeline definition (ordered) ------------------------------------------
steps <- list(
  list(id = "01", script = "01_cohort.R",                         label = "Cohort construction"),
  list(id = "02", script = "02_clinical_inflection.R",            label = "Clinical inflection & code status"),
  list(id = "03", script = "03_aim1.R",                           label = "Aim 1: discharge patterns & competing-risks CIF"),
  list(id = "04", script = "04_aim2a3_risk_adjusted_variation.R", label = "Aim 2A: risk-standardized variation"),
  list(id = "05", script = "05_aim3_clinical_trajectories.R",     label = "Aim 3: documented trajectories")
)

# ---- argument parsing -------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
flagval  <- function(flag) { i <- which(args == flag); if (length(i) && i < length(args)) args[i + 1] else NA_character_ }
keepgoing <- "--keep-going" %in% args
ids       <- vapply(steps, `[[`, "", "id")

if ("--list" %in% args) {
  cat("\nPipeline stages:\n")
  for (s in steps) cat(sprintf("  %s  %-40s %s\n", s$id, s$script, s$label))
  quit(save = "no")
}

only <- flagval("--only"); from <- flagval("--from")
keep <- rep(TRUE, length(steps))
if (!is.na(only)) keep <- ids %in% trimws(strsplit(only, ",")[[1]])
if (!is.na(from)) {
  m <- match(from, ids)
  if (is.na(m)) stop("--from stage '", from, "' is not a known stage id.", call. = FALSE)
  keep <- seq_along(steps) >= m
}
steps <- steps[keep]
if (!length(steps)) stop("No stages selected.", call. = FALSE)

# ---- elapsed-time formatter -------------------------------------------------
fmt <- function(sec) {
  if (is.na(sec)) return("       -")
  h <- sec %/% 3600; m <- (sec %% 3600) %/% 60; s <- sec %% 60
  if (h > 0)      sprintf("%dh %02dm %02ds", h, m, as.integer(round(s)))
  else if (m > 0) sprintf("%dm %02ds", m, as.integer(round(s)))
  else            sprintf("%.1fs", s)
}
rule <- function(ch = "-", n = 74) cat(strrep(ch, n), "\n", sep = "")

# ---- run --------------------------------------------------------------------
cat("\n"); rule("=")
cat(sprintf("CLIF hospice pipeline\n  project : %s\n  R %s | %s\n  started : %s\n",
            project_root, getRversion(), R.version$platform,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
rule("=")

log       <- data.frame(step = character(), stage = character(),
                        seconds = numeric(), status = character(),
                        stringsAsFactors = FALSE)
t_all     <- Sys.time()
failed    <- FALSE

for (s in steps) {
  path <- file.path(code_dir, s$script)
  if (!file.exists(path)) {
    cat(sprintf("\n==> [%s] %s\n    SKIP: %s not found\n", s$id, s$label, s$script))
    log <- rbind(log, data.frame(step = s$id, stage = s$label, seconds = NA_real_, status = "missing"))
    failed <- TRUE
    if (!keepgoing) break else next
  }
  cat(sprintf("\n==> [%s] %s\n", s$id, s$label)); flush.console()
  t0   <- Sys.time()
  code <- system2(rscript, args = shQuote(path), stdout = "", stderr = "")
  el   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok   <- isTRUE(code == 0)
  cat(sprintf("    -> %s  (%s)\n", if (ok) "OK" else sprintf("FAILED (exit %s)", code), fmt(el)))
  log  <- rbind(log, data.frame(step = s$id, stage = s$label, seconds = el,
                                status = if (ok) "ok" else "FAILED"))
  if (!ok) { failed <- TRUE; if (!keepgoing) break }
}
total <- as.numeric(difftime(Sys.time(), t_all, units = "secs"))

# ---- summary ----------------------------------------------------------------
cat("\n"); rule("=")
cat("TIMING SUMMARY\n"); rule("-")
cat(sprintf("  %-4s %-46s %12s   %s\n", "step", "stage", "elapsed", "status"))
for (i in seq_len(nrow(log)))
  cat(sprintf("  %-4s %-46s %12s   %s\n",
              log$step[i], substr(log$stage[i], 1, 46), fmt(log$seconds[i]), log$status[i]))
rule("-")
cat(sprintf("  %-4s %-46s %12s\n", "", "TOTAL", fmt(total)))
rule("=")

# ---- persist timing ---------------------------------------------------------
out <- rbind(
  transform(log, seconds = round(seconds, 1)),
  data.frame(step = "", stage = "TOTAL", seconds = round(total, 1),
             status = if (failed) "FAILED" else "ok")
)
tl <- file.path(project_root,
                sprintf("pipeline_timing_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))
ok_write <- tryCatch({ write.csv(out, tl, row.names = FALSE); TRUE }, error = function(e) FALSE)
if (ok_write) cat(sprintf("\nTiming log: %s\n", tl))

if (failed) {
  cat("\nPipeline FAILED. Fix the failing stage and re-run with --from <stage>.\n")
  quit(save = "no", status = 1)
}
cat("\nPipeline completed successfully.\n")
quit(save = "no", status = 0)
