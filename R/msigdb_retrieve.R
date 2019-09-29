.MSIG <- list()

#' Retrieves a msgidb collection data.frame for the given species
#'
#' This function retrieves arbitrary subset of the MSigDB collections, which are
#' determined by the values passed into `collections`. This parameter is
#' a character vector that can be any of the collections themselves (ie.
#' "H", C1:C7).
#'
#' The C2 collection includes subets of curated databases in it, such as
#' biocarta, kegg, pid, and reactome. You can also use these identifiers in the
#' `collections` parameter if you don't want to retrieve all of C2. When the
#' `promote_subcategory_to_collection` is `TRUE`, these databases will be pulled
#' out of the C2 collection, and promoted to their own collections, themselves.
#' The C5/GO annotations will also be split up into three collecionts, GO_MF,
#' GO_CC, GO_BP.
#'
#' @export
#' @param collections character of MSigDb collections ("h", "c1", ..., "c7")
#' @param species the name of the species (human, mouse, etc.) see entries
#'   in GeneSetDb.MSigDB:::.species_tbl()[["common_name]] (+ "human")
#' @param id_type "ensembl", "entrez", or "symbol"
#' @param version if `NULL` (default), the latest will be returned.
#' @param allow_multimap If `TRUE` (default), multiple entries are returned
#'   for a single gene in the original dataset, in the event that a single
#'   human,entrez identifier maps to more than one ensembl identifier.
#' @param min_ortho_sources The minimum number of databases that need to support
#'   the ortholog mapping. This number is provided by the [hcop()] database.
#' @param min_ortho_sources Filter the hcop table on `num_sources` column to
#'   ensure a minimum number of databases that support the ortholog map
#' @param promote_subcategory_to_collection When `TRUE` (default is `FALSE`),
#'   the database collections in C2 (like reactome, biocarta) will be pulled out
#'   of the c2 collection and promoted to their own (ie. there will be
#'   "reactome_c2" collection).
#' @return a geneset data.frame with the msigdb collecitons mapped to the given
#'   species. This result can be passed into `multiGSEA::GeneSetDb()` to get
#'   gene set mojo started.
msigdb_retrieve <- function(collections = "H", species = "human",
                            id_type = c("ensembl", "entrez", "symbol"),
                            version = NULL, slim = TRUE,
                            allow_multimap = TRUE, min_ortho_sources = 2,
                            promote_subcategory_to_collection = FALSE,
                            cache = TRUE, ...) {
  if (is.null(collections)) {
    collections <- c("H", paste0("C", 1:7))
  }
  sinfo <- species_lookup(species)
  id_type <- match.arg(id_type)

  all.colls <- c("H", paste0("C", 1:7))
  all.dbs <- c("BIOCARTA", "KEGG", "PID", "REACTOME")
  collections <- assert_subset(toupper(collections), c(all.colls, all.dbs))
  collections <- unique(collections)

  if (sinfo$common_name != "human" && "c1" %in% collections) {
    warning("The C1 collection doesn't make sense for anything but human, ...",
            immediate. = TRUE)
  }

  db.all <- msigdb_load(version = version, cache = cache)
  mversion <- attr(db.all, "msigdb_version")
  db.all <- rename(db.all, name = "gs_name", collection = "gs_cat",
                   subcategory = "gs_subcat", msigdb_id = "gs_id")

  colls <- intersect(collections, all.colls)
  cdbs <- intersect(collections, all.dbs)

  db.all <- db.all %>%
    filter(collection %in% colls | subcategory %in% paste0("CP:", cdbs))

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

  if (promote_subcategory_to_collection) {
    # we aren't doing all of them!
    cp <- paste0("CP:", cdbs)
    go <- c("BP", "CC", "MF")
    cp.df <- filter(out, subcategory %in% cp)
    go.df <- filter(out, subcategory %in% go)

    rest <- out
    if (nrow(cp.df)) {
      rest <- anti_join(rest, cp.df, by = c("collection", "name"))
      cp.df <- cp.df %>%
        mutate(collection = sub("CP:", "", subcategory),
               name = sub("^.*?_" ,"", name),
               subcategory = "C2")
    }
    if (nrow(go.df)) {
      rest <- anti_join(rest, go.df, by = c("collection", "name"))
      go.df <- go.df %>%
        mutate(collection = paste0("GO_", subcategory),
               name = sub("^GO_" , "", name),
               subcategory = "C5")
    }
    out <- bind_rows(rest, cp.df, go.df)
  }

  attr(out, "species_info") <- sinfo
  attr(out, "msigdb_version") <- mversion
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
  mdb <- .CACHE[["msigdb"]][[version.]]
  if (is.null(mdb)) {
    mdb <- readRDS(msig.info$path)
    attr(mdb, "msigdb_version") <- version.
    if (cache) {
      mcache <- .CACHE[["msigdb"]]
      mcache[[version.]] <- mdb
      assign("msigdb", mcache, envir = .CACHE)
    }
  }
  mdb
}
