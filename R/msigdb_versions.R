#' Returns a table of available msigdb versions within this package.
#'
#' @export
#' @noRd
msigdb_versions <- function(dat.dir = NULL,
                               dir.pattern = "msigdb-collection_v.*rds$",
                               version.regex = "(\\d+)\\.(\\d+)\\.(\\d+)?") {
  if (is.null(dat.dir)) {
    dat.dir <- system.file("extdata", package = "GeneSetDb.MSigDB")
  }
  assert_directory_exists(dat.dir, "r")

  coll.fns <- dir(dat.dir, pattern = dir.pattern)
  if (length(coll.fns) == 0) {
    stop("No msigdb-collection*.rds file found in: ", dat.dir)
  }

  parsed <- stringr::str_match(coll.fns, version.regex)
  out <- tibble(
    path = file.path(dat.dir, coll.fns),
    version = parsed[, 1],
    major = as.integer(parsed[, 2]),
    minor = as.integer(parsed[, 3]),
    patch = as.integer(parsed[, 4]))
  out$patch <- ifelse(is.na(out$patch), 0, out$patch)
  dplyr::arrange(out, desc(major), desc(minor), desc(patch))
}
