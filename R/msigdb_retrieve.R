.MSIG <- list()

#' Loads the reference msigdb collection data.frame
#'
#' These are all-human identifiers.
#'
#' Info for v7: http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_v7.0_Release_Notes
#'
#' @export
#' @param version Specifies the version of the database to load. If `NULL`
#'   (the default), the latest version is retrieved. To see what other versions
#'   are available, refer to the `version` column in the [msigdb_versions()]
#'   data.frame
#' @return The full data.frame of gene set collections for the give version.
#'   The version is added as a `msigdb_version` attribute to the returned
#'   result.
msigdb_load <- function(version = NULL, cache = TRUE) {
  msig.info <- head(msigdb_versions(), 1)

  if (!is.null(version)) {
    warning("Skipping version for now, always pulling in latest")
  }
  # TODO: hack logic in here to support version selection
  version. <- msig.info$version
  mdb <- .MSIG[[version.]]
  if (is.null(mdb)) {
    mdb <- readRDS(msig.info$path)
    attr(mdb, "msigdb_version") <- version.
    if (cache) {
      .MSIG[[version.]] <<- mdb
    }
  }
  mdb
}

#' Retrieves a msgidb collection data.frame for the given species
#'
#' @export
#' @param species the name of the species (human, mouse, etc.) see entries
#'   in GeneSetDb.MSigDB:::.species_tbl()[["common_name]] (+ "human")
#' @param collections character of MSigDb collections ("h", "c1", ..., "c7")
#' @param id_type "ensembl", "entrez", or "symbol"
#' @param version if `NULL` (default), the latest will be returned.
#' @param min_ortho_sources Filter the hcop table on `num_sources` column to
#'   ensure a minimum number of databases that support the ortholog map
#' @return a geneset data.frame with the msigdb collecitons mapped to the given
#'   species. This result can be passed into `multiGSEA::GeneSetDb()` to get
#'   gene set mojo started.
msigdb_retrieve <- function(species, collections = NULL,
                            id_type = c("ensembl", "entrez", "symbol"),
                            version = NULL, slim = TRUE,
                            allow_multimap = TRUE, min_ortho_sources = 2,
                            cache = TRUE, ...) {
  id_type <- match.arg(id_type)
  sinfo <- species_lookup(species)
  all.colls <- c("H", paste0("C", 1:7))
  if (is.null(collections)) {
    collections <- all.colls
  }
  collections <- unique(assert_subset(toupper(collections), all.colls))
  if (sinfo$common_name != "human" && "c1" %in% collections) {
    warning("The C1 collection doesn't make sense for anything but human, ...",
            immediate. = TRUE)
  }

  db.all <- msigdb_load(version = version, cache = cache)
  db.all <- rename(db.all, name = "gs_name", collection = "gs_cat",
                   subcategory = "gs_subcat", msigdb_id = "gs_id")
  db.all <- filter(db.all, collection %in% collections)

  gene.cols <- c("human_entrez_symbol", "human_entrez_id",
                 "human_ensembl_id", "human_ensembl_symbol")

  # human is a special case, since there is no homology mapping required
  if (sinfo[["common_name"]] == "human") {
    if (id_type == "ensembl") {
      fid.col <- "human_ensembl_id"
      sym.col <- "human_ensembl_symbol"
    } else if (id_type == "entrez") {
      fid.col <- "human_entrez_id"
      sym.col <- "human_entrez_symbol"
    } else {
      fid.col <- "human_entrez_symbol"
      sym.col <- "human_entrez_id"
    }
    out <- db.all %>%
      rename(featureId = fid.col, symbol = sym.col) %>%
      select(collection, name, featureId, symbol, subcategory, everything()) %>%
      filter(nchar(featureId) > 0, !is.na(featureId), featureId != "-") %>%
      distinct(collection, name, featureId, .keep_all = TRUE)
    if (slim) {
      out <- select(out, collection:subcategory)
    }
  } else {
    out <- .msigdb_map_ortho(db.all, id_type, sinfo, slim,
                             allow_multimap, min_ortho_sources)
  }

  out <- out %>%
    mutate(subcategory = ifelse(nchar(subcategory) == 0, NA, subcategory))
  attr(out, "species_info") <- sinfo
  out
}

#' This function will return the db.all data.frame with a featureId
#' and symbol column.
#'
#' We assume all parameter checking is already done.
#'
#' Joining on entrez id for now. Seems more primary to the internal Msigdb
#' mappings.
#'
#' @noRd
#' @param allow_multimap Orthologos mapping of a gene from human to whatever
#'   can result in multiple genes. If you want to allow 1:many mapping, leave
#'   this unchaged (`TRUE`), otherwise set to `FALSE` to only include one
#'   mapping. The map was already filtered to only include those with the
#'   highest number of database support, but this can still be more than one.
#' @param min_ortho_sources Filter the hcop table on `num_sources` column to
#'   ensure a minimum number of databases that support the ortholog map
.msigdb_map_ortho <- function(db.all, id_type, sinfo, slim,
                              allow_multimap = TRUE,
                              min_ortho_sources = 2) {
  hcop.all <- filter(hcop(), num_sources >= min_ortho_sources)
  hcop. <- filter(hcop.all, species_id == sinfo[["species_id"]])
  if (!allow_multimap) {
    hcop. <- distinct(hcop., human_entrez_gene, .keep_all = TRUE)
  }
  out <- inner_join(db.all, hcop.,
                    by = c("human_entrez_id" = "human_entrez_gene"))

  if (id_type == "ensembl") {
    fid.col <- "ensembl_gene"
    sym.col <- "gene_symbol"
  } else if (id_type == "entrez") {
    fid.col <- "entrez_gene"
    sym.col <- "gene_symbol"
  } else {
    fid.col <- "gene_symbol"
    sym.col <- "entrez_gene"
  }

  out <- out %>%
    rename(featureId = fid.col, symbol = sym.col) %>%
    filter(nchar(featureId) > 0, !is.na(featureId), featureId != "-") %>%
    transmute(collection, name, featureId, symbol, subcategory,
              num_sources, human_ensembl_id,
              human_symbol = human_entrez_symbol)

  if (slim) {
    out <- select(out, collection:num_sources)
  }

  out <- arrange(out, collection, name, featureId, desc(num_sources))
  distinct(out, collection, name, featureId, .keep_all = TRUE)
}
