# rmvpa DSL: the single `rmvpa` command.

.dsl_cli_usage <- function() {
  c(
    "usage: rmvpa <command> [arguments]",
    "",
    "  new <bids_root> [--dir DIR] [--target COLUMN]   create rmvpa.yml and a first analysis",
    "  add <verb> <name> [--scope searchlight|regions|whole_brain]   add analyses/<name>.yml",
    "  explain <analysis>                              print the resolved analysis with origins",
    "  check [analysis...] [--subjects 01,02]          validate without fitting; estimate cost",
    "  draft [analysis...]                             one subject, randomized searchlight, no permutations",
    "  run [analysis...] [--subjects 01,02] [--force] [--yes]   check, then run changed cells",
    "  status [analysis...]                            done / stale / missing per subject",
    "  lock [analysis...] [--note TEXT]                append a preregistration record",
    "  migrate <old_config.yaml> [--dir DIR]           convert an rmvpa-searchlight/regional config",
    "  studio [--port N] [--no-browser] [--background IMG] [--surface FILE]   open the studio in the browser",
    "  export <dir> [analysis...] [--background IMG] [--surface FILE] [--force]   write a static studio folder",
    "  verbs | classifiers | version | help",
    "",
    "Run inside a project directory, or pass --project DIR."
  )
}

.dsl_cli_parse <- function(argv) {
  flags <- list(); pos <- character()
  value_flags <- c("dir", "target", "subjects", "scope", "note", "project", "port", "background", "surface")
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (grepl("=", key, fixed = TRUE)) {
        flags[[sub("=.*$", "", key)]] <- sub("^[^=]*=", "", key)
      } else if (key %in% value_flags) {
        if (i == length(argv)) stop("--", key, " needs a value", call. = FALSE)
        flags[[key]] <- argv[i + 1L]; i <- i + 1L
      } else {
        flags[[key]] <- TRUE
      }
    } else {
      pos <- c(pos, a)
    }
    i <- i + 1L
  }
  unknown <- setdiff(names(flags), c(value_flags, "force", "yes", "help", "no-browser"))
  if (length(unknown)) stop("unknown option --", unknown[1], call. = FALSE)
  list(cmd = if (length(pos)) pos[1] else "help", args = pos[-1], flags = flags)
}

.dsl_cli_split <- function(x) if (is.null(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1]]

.dsl_cli_confirm <- function(passes, yes) {
  if (isTRUE(yes) || passes <= 50) return(TRUE)
  msg <- sprintf("This run needs about %s searchlight/ROI passes. Continue? [y/N] ", format(passes, big.mark = ","))
  if (!isatty(stdin())) {
    message(msg, "(no terminal; continuing. Pass --yes to silence this.)")
    return(TRUE)
  }
  cat(msg)
  ans <- readLines(file("stdin"), n = 1L)
  length(ans) == 1L && tolower(trimws(ans)) %in% c("y", "yes")
}

#' Command-line entry point for the rmvpa DSL
#'
#' Implements the `rmvpa` command (see `rmvpa help`). The launcher installed
#' by [install_rmvpa_command()] calls this function; it can also be called from
#' R with a character vector.
#'
#' @param argv Command-line arguments, without the program name.
#' @return An exit status (0 on success), invisibly.
#' @export
rmvpa_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
  status <- tryCatch({
    p <- .dsl_cli_parse(argv)
    f <- p$flags
    proj <- f$project %||% "."
    analyses <- if (length(p$args)) p$args else NULL
    switch(p$cmd,
      help = , `--help` = cat(.dsl_cli_usage(), sep = "\n"),
      version = cat(sprintf("rmvpa DSL %s (rmvpadsl %s, rMVPA %s)\n", paste(.dsl_versions, collapse = ", "),
                            utils::packageVersion("rmvpadsl"), utils::packageVersion("rMVPA"))),
      verbs = {
        v <- .dsl_verbs()
        for (n in names(v)) cat(sprintf("%-14s %-9s scopes: %s\n", n, v[[n]]$family,
                                        paste(.dsl_scope_support(n), collapse = ", ")))
      },
      classifiers = cat(rmvpa_classifiers(), sep = "\n"),
      new = {
        if (!length(p$args)) stop("usage: rmvpa new <bids_root> [--dir DIR] [--target COLUMN]", call. = FALSE)
        rmvpa_new(p$args[1], dir = f$dir %||% ".", target = f$target)
      },
      add = {
        if (length(p$args) != 2L) stop("usage: rmvpa add <verb> <name>", call. = FALSE)
        rmvpa_add(proj, p$args[1], p$args[2], scope = f$scope)
      },
      explain = {
        if (length(p$args) != 1L) stop("usage: rmvpa explain <analysis>", call. = FALSE)
        rmvpa_explain(proj, p$args[1])
      },
      check = {
        chk <- rmvpa_check(proj, analyses, subjects = .dsl_cli_split(f$subjects))
        if (!isTRUE(attr(chk, "ok"))) stop("check failed", call. = FALSE)
      },
      draft = rmvpa_run(proj, analyses, draft = TRUE, force = TRUE),
      run = {
        chk <- rmvpa_check(proj, analyses, subjects = .dsl_cli_split(f$subjects), quiet = TRUE)
        if (!isTRUE(attr(chk, "ok"))) { print(chk); stop("check failed; nothing was run", call. = FALSE) }
        if (!.dsl_cli_confirm(attr(chk, "passes"), f$yes)) stop("cancelled", call. = FALSE)
        res <- rmvpa_run(proj, analyses, subjects = .dsl_cli_split(f$subjects), force = isTRUE(f$force),
                         check = FALSE)
        print(res, row.names = FALSE)
      },
      status = print(rmvpa_status(proj, analyses), row.names = FALSE),
      lock = rmvpa_lock(proj, analyses, note = f$note),
      studio = rmvpa_studio(proj, port = if (!is.null(f$port)) as.integer(f$port), launch = !isTRUE(f[["no-browser"]]),
                            background = f$background, surfaces = .dsl_cli_split(f$surface), block = TRUE),
      export = {
        if (!length(p$args)) stop("usage: rmvpa export <dir> [analysis...]", call. = FALSE)
        rmvpa_export_studio(proj, p$args[1], analyses = if (length(p$args) > 1) p$args[-1],
                            background = f$background, surfaces = .dsl_cli_split(f$surface), overwrite = isTRUE(f$force))
      },
      migrate = {
        if (length(p$args) != 1L) stop("usage: rmvpa migrate <old_config.yaml> [--dir DIR]", call. = FALSE)
        rmvpa_migrate(p$args[1], dir = f$dir %||% ".")
      },
      stop(sprintf("unknown command '%s'%s. Run `rmvpa help`.", p$cmd,
                   .dsl_did_you_mean(p$cmd, c("new", "add", "explain", "check", "draft", "run", "status",
                                              "lock", "migrate", "studio", "export", "verbs", "classifiers", "version", "help"))),
           call. = FALSE)
    )
    0L
  }, error = function(e) {
    message("rmvpa: ", conditionMessage(e))
    1L
  })
  invisible(status)
}
