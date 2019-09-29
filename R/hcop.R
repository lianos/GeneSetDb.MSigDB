#' Retrieves the hcop ortholog map data.frame
#'
#' This is the data.frame that is used to map the collection data from its
#' source/raw input (human, entrez), to ensembl or other species.
#'
#' HCOPS is the HGNC Comparison of Orthology Predictions.
#'
#' @export
#' @seealso https://www.genenames.org/help/hcop/
#'
#' @param cache If `TRUE` (default), will store the data.frame internally so
#'   that it doesn't have to be re-read the next time the function is called
#' @return a data.frame with human to species ortholog mapping info
hcop <- function(cache = TRUE) {
  out <- .CACHE$hcop
  if (is.null(out)) {
    hfn <- system.file("extdata", "hcop.rds", package = "GeneSetDb.MSigDB")
    out <- readRDS(hfn)
    if (cache) {
      assign("hcop", out, envir = .CACHE)
    }
  }
  out
}


#' Retrieves and parses the latest HCOP data
#'
#' HCOPS is the HGNC Comparison of Orthology Predictions.
#' This code was taken and modifed from the msigdbr::msigdbr-prepare.R script
#' We would have just stuck with that, however they dump the ensembl identifiers
#' and I want to keep them
#'
#' We use the human hcop file because MSigDB's genesets (as far as I understand)
#' are referenced to human identifiers first.
#'
#' This script is not exported on purpose
#'
#' @seealso https://www.genenames.org/help/hcop/
#' @noRd
generate_hcop_orthologs <- function(hcop_txt_url = NULL) {
  if (!test_file_exists(hcop_txt_url)) {
    message("Downloading hcop file ================================")
    hcop_txt_url = file.path(
      "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop",
      "human_all_hcop_sixteen_column.txt.gz")
  }
  hcop = readr::read_tsv(hcop_txt_url)

  # Keep mapping only genes found in multiple ortholog/homolog databases
  msigdbr_orthologs =
    hcop %>%
    dplyr::select(
      human_entrez_gene,
      human_ensembl_gene,
      human_gene_symbol = human_symbol,
      species_id = ortholog_species,
      entrez_gene = ortholog_species_entrez_gene,
      ensembl_gene = ortholog_species_ensembl_gene,
      gene_symbol = ortholog_species_symbol,
      sources = support
    ) %>%
    dplyr::filter(
      human_entrez_gene != "-",
      entrez_gene != "-",
      gene_symbol != "-"
    ) %>%
    dplyr::mutate(
      human_entrez_gene = as.integer(human_entrez_gene),
      entrez_gene = as.integer(entrez_gene),
      num_sources = stringr::str_count(sources, ",") + 1
    ) %>%
    dplyr::filter(
      # human_entrez_gene %in% msigdb_entrez_genes,
      # num_sources > 1
      num_sources >= 0
    )

  # Names and IDs of common species
  species_tbl <- species_info(hcop_only = FALSE) %>%
    filter(common_name != "human") %>%
    select(species_id, species_name)

  # Add species names
  msigdbr_orthologs = dplyr::inner_join(
    species_tbl,
    msigdbr_orthologs, by = "species_id")

  # Keep all of the mappings from the "top tier" mapping, ie. the ortholog
  # mapping that is supported by the most databases
  orthologs =
    msigdbr_orthologs %>%
    dplyr::group_by(human_entrez_gene, species_name) %>%
    dplyr::filter(num_sources == max(num_sources)) %>%
    # dplyr::top_n(1, num_sources) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(human_entrez_gene = as.character(human_entrez_gene),
           entrez_gene = as.character(entrez_gene),
           species_id = as.integer(species_id))
  attr(orthologs, "creation_date") <- format(Sys.time())
  orthologs
}
