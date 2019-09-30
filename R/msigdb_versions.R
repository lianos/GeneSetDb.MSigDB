#' Returns a table of available msigdb versions within this package.
#'
#' The parameters in this function really shouldn't be used. By default
#' it checks the internal package directory to see what msigdb definitions
#' are included.
#'
#' @export
#' @param dat.dir the directory that holds the version data. By default it
#'   looks at the package's internal `"extdata"` pacakge for the collection
#'   data.
#' @param dir.pattern what the misgdb files look like
#' @param version.regex how the versions are encoded
#' @return a tibble of msigdb version info
msigdb_versions <- function(dat.dir = NULL,
                            dir.pattern = "msigdb-collection_v.*rds$",
                            version.regex = "(\\d+)\\.(\\d+)\\.(\\d+)?") {
  if (is.null(dat.dir)) {
    dat.dir <- system.file("extdata", package = "msigdb.data")
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
