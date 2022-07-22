.onLoad <- function(libname, pkgname){
  if (interactive()) {
    logger::log_layout(logger::layout_glue_colors)
  }
}
