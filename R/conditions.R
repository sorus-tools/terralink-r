# Condition and logging helpers for TerraLink

terralink_abort <- function(message, class = "terralink_error", details = NULL, fix = NULL, call = NULL) {
  detail_items <- NULL
  fix_items <- NULL
  if (!is.null(details) && length(details) > 0) {
    detail_items <- stats::setNames(as.character(details), rep("x", length(details)))
  }
  if (!is.null(fix) && length(fix) > 0) {
    fix_items <- stats::setNames(as.character(fix), rep("i", length(fix)))
  }
  bullets <- c(message, detail_items, fix_items)
  cli::cli_abort(bullets, class = unique(c(class, "terralink_error")), call = call)
}

terralink_warn <- function(message, ctx = NULL, call = NULL) {
  if (!is.null(ctx) && is.environment(ctx)) {
    ctx$warnings <- c(ctx$warnings, as.character(message))
  }
  cli::cli_warn(message, call = call)
  invisible(message)
}

terralink_inform <- function(message, ctx = NULL, call = NULL, level = 1) {
  if (is.null(ctx) || !is.environment(ctx)) {
    cli::cli_inform(message, call = call)
    return(invisible(message))
  }
  verbose <- ctx$verbose %||% 0L
  if (verbose >= level) {
    cli::cli_inform(message, call = call)
  }
  invisible(message)
}

terralink_new_run_context <- function(verbose = 0, progress = FALSE) {
  ctx <- new.env(parent = emptyenv())
  ctx$warnings <- character(0)
  ctx$diagnostics <- list()
  ctx$verbose <- as.integer(verbose)
  ctx$progress <- isTRUE(progress)
  ctx$progress_id <- NULL
  ctx
}

terralink_progress_start <- function(ctx, total = 100, message = NULL) {
  if (is.null(ctx) || !is.environment(ctx) || !isTRUE(ctx$progress)) return(invisible(NULL))
  ctx$progress_id <- tryCatch(
    cli::cli_progress_bar(total = total, format = "{cli::pb_bar} {cli::pb_percent} {cli::pb_message}"),
    error = function(e) NULL
  )
  if (is.null(ctx$progress_id)) {
    ctx$progress <- FALSE
    return(invisible(NULL))
  }
  if (!is.null(message)) {
    terralink_progress_update(ctx, 0, message = message)
  }
  invisible(ctx$progress_id)
}

terralink_progress_update <- function(ctx, value, message = NULL) {
  if (is.null(ctx) || !is.environment(ctx) || !isTRUE(ctx$progress) || is.null(ctx$progress_id)) {
    return(invisible(NULL))
  }
  ok <- tryCatch({
    update_formals <- names(formals(cli::cli_progress_update))
    if (!is.null(message) && "message" %in% update_formals) {
      cli::cli_progress_update(ctx$progress_id, set = value, message = message)
    } else {
      cli::cli_progress_update(ctx$progress_id, set = value)
    }
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    ctx$progress_id <- NULL
    ctx$progress <- FALSE
    return(invisible(NULL))
  }
  invisible(value)
}

terralink_progress_done <- function(ctx) {
  if (is.null(ctx) || !is.environment(ctx) || !isTRUE(ctx$progress) || is.null(ctx$progress_id)) {
    return(invisible(NULL))
  }
  tryCatch(
    cli::cli_progress_done(ctx$progress_id),
    error = function(e) NULL
  )
  ctx$progress_id <- NULL
  invisible(TRUE)
}

terralink_format_bytes <- function(bytes) {
  if (!is.finite(bytes)) return("unknown")
  units <- c("B", "KB", "MB", "GB", "TB")
  pow <- if (bytes <= 0) 0 else floor(log(bytes, base = 1024))
  pow <- min(pow, length(units) - 1)
  value <- bytes / (1024 ^ pow)
  sprintf("%.1f %s", value, units[pow + 1])
}
