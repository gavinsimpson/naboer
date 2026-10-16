check_coverage <- function(coverage) {
  if (covr::percent_coverage(coverage) < 90) stop("Overall coverage must reach 90%.")
  data <- covr::tally_coverage(coverage, by = "line")
  selected <- grepl("/(folds|indices|nei|utils)\\.R$", data$filename)
  if (!any(selected)) stop("Core files were not included in the coverage report.")
  core <- 100 * mean(data$value[selected] > 0)
  if (core < 95) stop(sprintf("Core coverage %.1f%% must reach 95%%.", core))
  invisible(TRUE)
}
