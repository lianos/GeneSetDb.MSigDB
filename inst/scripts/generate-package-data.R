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
devtools::load_all(".")

# Parse and HCOP ===============================================================
# HCOP is the HGNC Comparison of Orthology Predictions:
# https://www.genenames.org/help/hcop/
#
# This file contains ortholog maps from human to many other species, which can
# be queried by entreaz or ensembl id. We use this instead of biomaRt since the
# table we parse and serialize isn't too big (~ 4mb), and we don't have to
# bang on the biomaRt server each time the user requests an ortholog map.

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

dat.dir <- system.file("inst", "extdata", package = "msigdb.data")
msigdb.fn <- file.path(
  dat.dir,
  sprintf("msigdb-collection_v%s.rds", msigdb.version))
sids.fn <- file.path(dat.dir, "hcop_species.csv")
# Retrieve and save updated HCOPS file for ortholog mapping

hcop.all <- msigdb.data:::generate_hcop_orthologs(hcop.txt.fn)
hcop.rds.fn <- file.path(dat.dir, "hcop.rds")

# 4mb as of 2019-09-10
# takes ~ 0.45s to load via readRDS
saveRDS(hcop.all, hcop.rds.fn)
sids <- distinct(hcop.all, species_id, species_name)

# Let's save the species_id's that are in hcop. This is used in the
# .spcecies_tbl(hcop_only = TRUE)
write.csv(sids, sids.fn, row.names = FALSE)

# Parse and save the human reference MSigDb collection =========================
# The data.frame of MSigDb collections from msigdbr provides gene symbols and
# entrez identifiers. We will append the human ensembl id's to this.
#
# We do not yet enforce a 1:1 map here, and the serialized MSigDb catalog will
# have duplicated entrez entries for some genesets wince one entrez id can map
# to multiple ensembl IDs.
#
# When there are multiple entrez -> ensembl maps, however, we will remove
# entries that map to alternate chromosomes.

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
ens.xref <- local({
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mart.result <- biomaRt::getBM(
    attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
    filters = "entrezgene_id",
    values = unique(hdf$human_entrez_id),
    mart = mart)

  mart.result %>%
    transmute(human_entrez_id = as.character(entrezgene_id),
              human_ensembl_id = ensembl_gene_id,
              human_ensembl_symbol = hgnc_symbol,
              chromosome = chromosome_name)
})

# annotate alt chrommosomes
ex <- ens.xref %>%
  mutate(alt_chromosome = grepl("_|\\.", chromosome)) %>%
  group_by(human_entrez_id) %>%
  mutate(n = n(), n_std = sum(!alt_chromosome), n_alt = sum(alt_chromosome)) %>%
  ungroup()

ex <- filter(ex, n == 1 | (n > 1 & !alt_chromosome))
oops <- setdiff(ex$human_entrez_id, ens.xref$human_entrez_id)
stopifnot(length(oops) == 0)

# now annotate multimatches where one match doens't have an ensembl symbol
#
# 1. when an entrez id maps to multiple ensemble's, and some of those ensembl
#    id's do not have symbols, we'll drop the ones with out symbols
# 2. if an entrez id only maps to ensembl id's that do not have symbols, then
#    we keep them all.
exs <- select(ex, human_entrez_id:chromosome) %>%
  mutate(nosymbol = nchar(human_ensembl_symbol) == 0) %>%
  # let's remove fusions/read through genes by categorizing them as nosymbol
  mutate(nosymbol = nosymbol | grepl("-", human_ensembl_symbol, fixed = TRUE)) %>%
  group_by(human_entrez_id) %>%
  mutate(n = n(), n_symbol = sum(!nosymbol), n_nosymbol = sum(nosymbol)) %>%
  ungroup()
exs.keep <- with(exs, {
  n == 1 | (n > 1 & (!nosymbol & n_symbol > 0) | (n_symbol == 0))
})
exs.cut <- exs[exs.keep,]
oops <- setdiff(exs.cut$human_entrez_id, ens.xref$human_entrez_id)
stopifnot(length(oops) == 0)

hdf.ens <- hdf %>%
  left_join(select(exs.cut, human_entrez_id:chromosome), by = "human_entrez_id")

oops <- setdiff(hdf.ens$gs_name, hdf$gs_name)
stopifnot(length(oops) == 0)

# got smaller
orig.count <- count(hdf, gs_name)
new.count <- count(hdf.ens, gs_name)
cmp.count <- inner_join(orig.count, new.count, by = "gs_name")
stopifnot(nrow(cmp.count) == nrow(orig.count))

smaller <- filter(cmp.count, n.x > n.y)
stopifnot(nrow(smaller) == 0)
bigger <- filter(cmp.count, n.x < n.y) %>% arrange(n.x)
head(bigger)

subset(hdf, gs_name == bigger$gs_name[1])
subset(hdf.ens, gs_name == bigger$gs_name[1])

# Map Gene Ontology identifiers to C5 genesets, this needs the original
# msigdb*.xml file
msigdb_xml <- "~/workspace/data/msigdb/msigdb_v7.0.xml"
library(xml2)
library(dplyr)

msigdb_doc <- read_xml(msigdb_xml)
geneset_records <- xml_find_all(msigdb_doc, xpath = ".//GENESET")
msigdbr_genesets <-
  tibble(
    gs_name             = xml_attr(geneset_records, attr = "STANDARD_NAME"),
    gs_id               = xml_attr(geneset_records, attr = "SYSTEMATIC_NAME"),
    gs_cat              = xml_attr(geneset_records, attr = "CATEGORY_CODE"),
    geo_id              = xml_attr(geneset_records, attr = "EXACT_SOURCE"))

go.info <- msigdbr_genesets %>%
  filter(gs_cat == "C5") %>%
  transmute(gs_name, gs_id, geo_id)

# replace gs_id in C5 category with GO:XXXXXXX id
c5.tmp <- hdf.ens %>%
  inner_join(go.info, by = c("gs_name", "gs_id")) %>%
  mutate(gs_id = geo_id) %>%
  select(-geo_id)

hdf.out <- hdf.ens %>%
  filter(gs_cat != "C5") %>%
  bind_rows(c5.tmp)

saveRDS(hdf.out, msigdb.fn)

