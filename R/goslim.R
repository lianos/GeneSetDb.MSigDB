#' @noRd
.slim.dir <- function() {
  system.file("extdata", "go_subsets", package = "msigdb.data")
}

#' Returns the available subsets/slims
#'
#' @export
slim_available <- function() {
  slims <- dir(.slim.dir(), "json$")
  tibble(
    name = sub(".*?_", "", sub(".json", "", slims)),
    path = file.path(.slim.dir(), slims))
}

#' Returns the ID's that map to a specific GO subset / slim
#'
#' @export
#' @importFrom rjson fromJSON
#' @return a vector of GO ids
slim_ids <- function(name = "generic") {
  all.slims <- slim_available()
  info <- filter(all.slims, .data[["name"]] == .env[["name"]])
  if (nrow(info) != 1L) {
    stop("name must be one of: ", all.slimes[["name"]])
  }
  slim.json <- fromJSON(file = info[["path"]])
  slim.nodes <- slim.json[[1]][[1]][["nodes"]]
  go.ids <- sapply(slim.nodes, function(x) {
    basename(x[["id"]])
  })

  # Every GO id in the edges is also found in the nodes
  if (FALSE) {
    slim.edges <- slim.json[[1]][[1]][["edges"]]
    ee = basename(unlist(slim.edges, recursive = TRUE))
    ee <- ee[grepl("^GO", ee)]
    e.xtra <- setdiff(ee, go.ids) # length() == 0
  }

  sub("_", ":", go.ids)
}
