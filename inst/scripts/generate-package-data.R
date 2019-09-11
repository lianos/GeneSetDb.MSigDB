# This script makes four files inst/extdata:
#
# 1. inst/extdata/hcop.rds, which is parsed from `hcop.txt.fn`. This table
#    holds ortholg mapping for species in species_info()
# 2. hcop_sepcies.csv -- lists the species_id and species_name of the organism
#    in the hcop ortholog file. This is used to match organism requests against
#    when asking for an msigdb object
# 3. msigdb-collection_v<VERSION>.rds
#
# To retrieve msigdb collections for other organims except human, we take the
# colelction from (3), merge against the orthologs for that organism in (1),
# and make sure that there are no duplicate collecetion,name,featureId entries.
library(msigdbr)
library(multiGSEA)
library(dplyr)
library(data.table)
devtools::load_all(".")


# hcop.txt.fn <- file.path(dat.dir, "human_all_hcop_sixteen_column.txt.gz")
hcop.txt.fn <- file.path("~/tmp", "human_all_hcop_sixteen_column.txt.gz")
if (!file.exists(hcop.txt.fn)) {
  curl::curl_download(
    url = file.path(
      "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop",
      "human_all_hcop_sixteen_column.txt.gz"),
    destfile = hcop.txt.fn)
}

msigdb.version <- packageVersion("msigdbr")

dat.dir <- system.file("inst", "extdata", package = "GeneSetDb.MSigDB")
msigdb.fn <- file.path(
  dat.dir,
  sprintf("msigdb-collection_v%s.rds", msigdb.version))
sids.fn <- file.path(dat.dir, "hcop_species.csv")
# Retrieve and save updated HCOPS file for ortholog mapping

hcop.all <- GeneSetDb.MSigDB:::generate_hcop_orthologs(hcop.txt.fn)
hcop.rds.fn <- file.path(dat.dir, "hcop.rds")

# 4mb as of 2019-09-10
# takes ~ 0.45s to load via readRDS
saveRDS(hcop.all, hcop.rds.fn)
sids <- distinct(hcop.all, species_id, species_name)

# Let's save the species_id's that are in hcop. This is used in the
# .spcecies_tbl(hcop_only = TRUE)
write.csv(sids, sids.fn, row.names = FALSE)

# Generate base human msigdb data.frame with both entrez and ensembl info
hdf <- local({
  tmp <- msigdbr::msigdbr(species = "Homo sapiens")
  stopifnot(all.equal(tmp$human_gene_symbol, tmp$gene_symbol))
  stopifnot(all(tmp$species_name == "Homo sapiens"))
  tmp %>%
    rename(human_entrez_symbol = human_gene_symbol,
           human_entrez_id = entrez_gene) %>%
    select(-species_name, -gene_symbol, -sources) %>%
    mutate(human_entrez_id = as.character(human_entrez_id))
})

# Add ensembl identifiers ... this will create duplicate gs_name,gene rows
# which will be filtered out in the end
hdf.ens <- local({
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mart.result <- biomaRt::getBM(
    attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
    filters = "entrezgene_id",
    values = unique(hdf$human_entrez_id),
    mart = mart)

  xx <- mart.result %>%
    transmute(human_entrez_id = as.character(entrezgene_id),
              human_ensembl_id = ensembl_gene_id,
              human_ensembl_symbol = hgnc_symbol)
  left_join(hdf, xx, by = "human_entrez_id")
})

saveRDS(hdf.ens, msigdb.fn)

