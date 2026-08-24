# Messaging helpers -----------------------------------------------------------
#
# scatools emits user-facing progress messages with cli, which raises them as R
# conditions so they are captured by evaluate/knitr and appear inside rendered
# rmarkdown and quarto documents (see #9).
#
# cli has no notion of log levels. This helper restores the one level scatools
# actually relied on when it used logger: DEBUG sat below logger's default INFO
# threshold, so debug output was silent unless explicitly enabled.

#' Emit a debug-level message
#'
#' Internal counterpart to [cli::cli_alert_info()] for chatty diagnostics that
#' should stay hidden during normal use. Silent unless the user opts in with
#' `options(scatools.debug = TRUE)`.
#'
#' @param ... Passed on to [cli::cli_alert_info()]. Supports cli inline markup
#'   and `{}` interpolation.
#' @param .envir Environment in which to evaluate `{}` expressions.
#'
#' @return `NULL`, invisibly. Called for its side effect.
#' @noRd
log_debug <- function(..., .envir = parent.frame()) {
  if (isTRUE(getOption("scatools.debug", FALSE))) {
    cli::cli_alert_info(..., .envir = .envir)
  }
  invisible(NULL)
}
