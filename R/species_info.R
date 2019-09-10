#' Species information
#'
#' @export
species_info <- function(hcop_only = TRUE) {
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
  if (isTRUE(hcop_only)) {
    sids.fn <- system.file("extdata", "hcop_species.csv",
                           package = "GeneSetDb.MSigDB",
                           mustWork = TRUE)
    sids <- read.csv(sids.fn, colClasses = c("integer", "character"))
    out <- filter(out, species_id %in% sids[["species_id"]])
  }

  out
}
