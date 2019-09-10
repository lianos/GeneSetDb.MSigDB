#' Generates a species specific MSigDB GeneSetDb
#'
#' @noRd
.generate_msigdb_gdb <- function(gdb.dir = "~/workspace/data/GeneSetDb",
                                 species = c("human", "mouse", "fly"),
                                 id.type = c("ensembl", "entrez"),
                                 debug = FALSE) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr package required to build MSigDB GeneSetDb objects")
  }
  assert_directory_exists(gdb.dir, "w")
  hcop.fn <- file.path(gdb.dir, "hcop.rds")
  if (!test_file_exists(hcop.fn, "r")) {
    stop("Need to generate `hcop.rds` in gdb.dir, please run: ",
         sprintf('`multiGSEA:::generate_hcop_orthologs("%s")` first', gdb.dir))
  }
  species <- match.arg(species)
  species. <- sub("_", " ", resolve.species(species))

  hdf <- msigdbr::msigdbr(species = "Homo sapiens")

  xdf <- dplyr::transmute(
    hdf,
    collection = gs_cat, name = gs_name,
    featureId = as.character(entrez_gene),
    symbol = human_gene_symbol,
    symbol_source = gene_symbol,
    subcategory = ifelse(nchar(gs_subcat) == 0, NA, gs_subcat),
    systematic_name = gs_id)

  if (species. == "Homo sapiens") {
    # source and symbol_source are the same for human
    if (id.type == "entrez") {
      out <- select(xdf, -symbol_source)
    } else if (id.type == "ensembl") {
      # convert IDs with biomaRt
      mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      exref <- biomaRt::getBM(
        attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol"),
        filters = "entrezgene_id",
        values = unique(out$featureId),
        mart = mart) %>%
        transmute(entrez_gene = as.character(entrezgene_id),
                  ensembl_gene = ensembl_gene_id, symbol = hgnc_symbol)
      exref.u <- exref %>%
        dplyr::filter(!is.na(symbol), nchar(symbol) > 0) %>%
        dplyr::distinct(entrez_gene, symbol, .keep_all = TRUE)
      out <- inner_join(
        select(xdf, -symbol),
        exref.u, by = c("featureId" = "entrez_gene"))
      out <- transmute(
        out,
        collection, name, featureId = ensembl_gene, symbol, symbol_source,
        subcategory, systematic_name)
    }
  } else {
    hcop.all <- readRDS(hcop.fn)
    hcop. <- hcop.all %>%
      dplyr::filter(species_name == species., num_sources > 1) %>%
      dplyr::select(featureId = human_entrez_gene, entrez_gene,
                    ensembl_gene, symbol = gene_symbol)
    ortho <- dplyr::inner_join(
      dplyr::select(xdf, -symbol),
      hcop., by = c("featureId"))
    ortho <- select(ortho, -featureId)

    if (id.type == "ensembl") {
      ortho <- dplyr::rename(ortho, featureId = "ensembl_gene")
    } else {
      ortho <- dplyr::rename(ortho, featureId = "entrez_gene")
    }

    out <- dplyr::transmute(
      ortho,
      collection, name, featureId, symbol, symbol_source,
      subcategory, systematic_name)
  }
  # check what genesets are missing
  gmissing <- setdiff(xdf$name, out$name)
  if (length(gmissing)) {
    warning(length(gmissing), " gene sets missing after converion, set ",
            "`debug = TRUE` to see which ones")
    if (debug) browser()
  }
  # Check that there are no duplicate collection,name,featureId pairs
  gstats <- out %>%
    dplyr::group_by(collection, name, featureId) %>%
    dplyr::summarize(n = n())
  multimap <- dplyr::filter(gstats, n > 1)
  if (any(gstats$n > 1)) {
    multi <- dplyr::distinct(multimap, name, .keep_all = TRUE)
    warning(nrow(multi), " gene sets have multimap genes id's, set ",
            "`debug = TRUE` to see which ones")
    if (debug) browser()
    out <- dplyr::distinct(out, collection, name, featureId, .keep_all = TRUE)
  }

  gdb <- GeneSetDb(out)
  msigdb.version <- packageVersion("msigdbr")
  if (id.type == "ensembl") {
    id.fn <- GSEABase::ENSEMBLIdentifier()
  } else {
    id.fn <- GSEABase::EntrezIdentifier()
  }

  url.fn <- function(collection, name) {
    url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
    sprintf(url, name)
  }
  for (col in unique(geneSets(gdb)$collection)) {
    geneSetCollectionURLfunction(gdb, col) <- url.fn
    featureIdType(gdb, col) <- id.fn
    gdb <- addCollectionMetadata(
      gdb, col, 'source',
      sprintf("msigdbr_v%s", msigdb.version))
  }

  sname.parts <- unlist(strsplit(species., " "))
  sname <- paste0(substr(sname.parts[[1]], 1, 1), sname.parts[2])
  org(gdb) <- resolve.species(species)

  out.fn <- file.path(
    gdb.dir,
    sprintf("GeneSetDb.MSigDB.%s.v%s.rds", sname, msigdb.version))
  saveRDS(gdb, out.fn)
  invisible(gdb)
}
