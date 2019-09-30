#' Species information table and lookup
#'
#' `species_info()` returns a table of ncbi species id's, names, common name,
#' etc. `species_lookup()` performs fuzzy matches against the table to return
#' the species-level info for a query. The fuzzy match must start left-to-right.
#'
#' @export
#' @param hcop_only Only retrun species information for organisms that are
#'   in the [hcop()] database.
#' @return A single-row tibble with (numeric) species id, name, etc.
#' @examples
#' species_lookup("human")
#' species_lookup("scer")
#' species_lookup("hsapiens")
species_info <- function(hcop_only = FALSE) {
  out <- tibble::tribble(
    ~species_id,    ~species_name,              ~common_name,  ~bioc_name,
    9606L,          "Homo sapiens",             "human",       "Hsapiens",
    4932L,          "Saccharomyces cerevisiae", "yeast",       "Scerevisiae",
    6239L,          "Caenorhabditis elegans",   "celegans",    "Celegans",
    7227L,          "Drosophila melanogaster",  "fly",         "Dmelanogaster",
    7955L,          "Danio rerio",              "zebrafish",   "Drerio",
    9031L,          "Gallus gallus",            "chicken",     "Ggallus",
    9615L,          "Canis lupus familiaris",   "dog",         "Cfamiliaris",
    9823L,          "Sus scrofa",               "boar",        "Sscrofa",
    9913L,          "Bos taurus",               "cow",         "Btaurus",
    10090L,         "Mus musculus",             "mouse",       "Mmusculus",
    10116L,         "Rattus norvegicus",        "rat",         "Rnorvegicus",
    9541L,          "Macaca fascicularis",      "cyno",        "Mfascicularis")
  out <- mutate(out, species_name_ = sub(" ", "_", species_name))
  if (isTRUE(hcop_only)) {
    sids.fn <- system.file("extdata", "hcop_species.csv",
                           package = "msigdb.data",
                           mustWork = TRUE)
    sids <- read.csv(sids.fn, colClasses = c("integer", "character"))
    out <- filter(out, species_id %in% sids[["species_id"]])
  }
  out
}

#' @export
#' @rdname species_info
species_lookup <- function(query, stable = species_info(FALSE),
                           ignore.case = TRUE) {
  stopifnot(is.atomic(query), length(query) == 1L)
  row.idx <- NA_integer_
  for (cname in colnames(stable)) {
    vals <- stable[[cname]]
    row.idx <- grep(sprintf("^%s", query), vals, ignore.case = ignore.case)
    if (length(row.idx) > 1L) {
      stop(sprintf("ambiguous match for query `%s` -> `%s`",
                   query, paste(vals[row.idx], collapse = ",")))
    }
    if (length(row.idx) == 1L) break
  }
  if (is.na(row.idx)) {
    stop("couldn't match species against: ", query)
  }
  stable[row.idx,,drop=FALSE]
}
