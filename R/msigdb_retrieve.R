#' Retrieves a msgidb collection data.frame for the given species
#'
#' @export
#' @param species the name of the species (human, mouse, etc.) see entries
#'   in GeneSetDb.MSigDB:::.species_tbl()[["common_name]] (+ "human")
#' @param collections character of MSigDb collections ("h", "c1", ..., "c7")
#' @param id_type "ensembl", "entrez", or "symbol"
msigdb_retrieve <- function(species, collections = NULL,
                            id_type = c("ensembl", "entrez", "symbol")) {

}
