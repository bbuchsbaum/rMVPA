# Collate benchmark receipts into tools/bench/SCOREBOARD.md.
# Run from the package root:  Rscript tools/bench/scoreboard.R
#
# Uses the latest receipt per (tool, rMVPA branch, scenario, method). Ratio is
# rMVPA time / competitor time, so < 1 means rMVPA is faster. A ratio counts
# only when the output digests agree. Otherwise the row says the two sides
# computed different things.

`%||%` <- function(a, b) if (is.null(a)) b else a

receipt_dir <- file.path("tools", "bench", "receipts")
files <- list.files(receipt_dir, pattern = "[.]jsonl$", full.names = TRUE)
rows <- unlist(lapply(files, function(f) lapply(readLines(f), jsonlite::fromJSON)), recursive = FALSE)

flat <- do.call(rbind, lapply(rows, function(r) {
  data.frame(
    timestamp = r$timestamp,
    tool = r$tool,
    side = if (identical(r$tool, "rMVPA")) paste0(r$rmvpa_branch %||% "?", " @", r$rmvpa_git_sha %||% "?") else "competitor",
    scenario = r$scenario,
    method = r$method,
    engine = if (is.null(r$engine) || is.na(r$engine)) "" else as.character(r$engine),
    blas = as.character(r$blas),
    ms = r$ms_per_unit,
    iqr_ms = 1000 * r$iqr_s / r$n_units,
    unit = r$unit,
    digest = as.numeric(unlist(r$digest))[1],
    stringsAsFactors = FALSE
  )
}))
flat <- flat[order(flat$timestamp), ]
latest <- flat[!duplicated(flat[, c("side", "scenario", "method")], fromLast = TRUE), ]

comp <- latest[latest$side == "competitor", ]
ours <- latest[latest$side != "competitor", ]

fmt_ms <- function(x) ifelse(is.na(x), "—", ifelse(x < 1, sprintf("%.3f", x), sprintf("%.1f", x)))
lines <- c(
  "# Benchmark scoreboard",
  "",
  sprintf("Generated %s from `tools/bench/receipts/` by `tools/bench/scoreboard.R`.", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Single-threaded medians. **Ratio = rMVPA / competitor (< 1 means rMVPA is faster).**",
  "A ratio is shown only when the two output digests agree to 1e-4 relative;",
  "otherwise the methods are not computing the same thing and the row is marked.",
  "Losses are recorded as losses.",
  ""
)
for (sc in unique(latest$scenario)) {
  unit <- latest$unit[latest$scenario == sc][1]
  lines <- c(lines, sprintf("## %s (ms per %s)", sc, unit), "",
             "| method | rMVPA side | engine | rMVPA | competitor | tool | ratio | digests |",
             "|---|---|---|---|---|---|---|---|")
  for (i in which(ours$scenario == sc)) {
    o <- ours[i, ]
    c_ <- comp[comp$scenario == sc & comp$method == o$method, ]
    if (nrow(c_)) {
      agree <- isTRUE(abs(o$digest - c_$digest) <= 1e-4 * max(1, abs(c_$digest)))
      ratio <- if (agree) sprintf("**%.2f**", o$ms / c_$ms) else "n/a"
      dg <- if (agree) "match" else sprintf("differ (%.6g vs %.6g)", o$digest, c_$digest)
      lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |", o$method, o$side, o$engine,
                                fmt_ms(o$ms), fmt_ms(c_$ms), c_$tool, ratio, dg))
    } else {
      lines <- c(lines, sprintf("| %s | %s | %s | %s | — | — | — | rMVPA only |", o$method, o$side, o$engine, fmt_ms(o$ms)))
    }
  }
  lines <- c(lines, "")
}
blas <- unique(paste(latest$tool, latest$blas, sep = ": "))
lines <- c(lines, "## Environment", "", paste0("- ", blas),
           "- R uses the reference BLAS and numpy uses Apple Accelerate, so BLAS-heavy rows favour the competitor.",
           "- PyMVPA 2.6.5 cannot be installed with numpy >= 1.26 / Python 3.12 (needs numpy.distutils); not benchmarked.",
           "- CoSMoMVPA needs MATLAB/Octave; not available.")
writeLines(lines, file.path("tools", "bench", "SCOREBOARD.md"))
cat(lines, sep = "\n")
