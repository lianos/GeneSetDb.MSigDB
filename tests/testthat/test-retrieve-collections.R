context("Exercise MSigDB Retrieval")

# The human entrez reference is the version of the data that is closest to its
# original form, so we'll use that as a standard to compare conversions to other
# identifiers and species.

test_that("human entrez and ensembl id conversion are roughly equal", {
  df.hentrez <- msigdb_retrieve("human", collections = c("h", "c2"),
                                id_type = "entrez")
  df.hensembl <- msigdb_retrieve("human", collections = c("h", "c2"),
                                 id_type = "ensembl", rm_meta = FALSE)
  # these should be ensembl genes
  expect_true(all(grepl("ENSG\\d+$", df.hensembl$featureId)))

  # we should have the same number of genesets among both collections
  og.stats <- count(df.hentrez, collection, name)
  hens.stats <- count(df.hensembl, collection, name)
  expect_equal(hens.stats$name, og.stats$name)

  # the average difference in gene set size should be minimal
  diff.size <- abs(og.stats$n - hens.stats$n)
  expect_true(all(diff.size) < 3)
})

test_that("allow_multimap = FALSE ensures (arbitrary) 1:1 ortholog mapping", {
  df.multi <- msigdb_retrieve("mouse", collections = c("h", "c2"),
                              id_type = "ensembl", rm_meta = FALSE)
  df.uniq <- msigdb_retrieve("mouse", collections = c("h", "c2"),
                             id_type = "ensembl", allow_multimap = FALSE,
                             rm_meta = FALSE)
  expect_true(nrow(df.uniq) < nrow(df.multi))

  # df.multi should have some map duplication
  multi.count <- df.multi %>%
    distinct(human_symbol, featureId) %>%
    count(human_symbol) %>%
    arrange(desc(n))
  expect_true(any(multi.count$n > 1))

  # df.uniq should only have one mapped featureId per human_symbol
  uniq.count <- df.uniq %>%
    distinct(human_symbol, featureId) %>%
    count(human_symbol) %>%
    arrange(desc(n))
  expect_true(all(uniq.count$n == 1))
})

test_that("ortholog mapping seems approximately correct", {
  df.og <- msigdb_retrieve("human", collections = c("h", "c2"),
                           id_type = "ensembl")
  df.mouse <- msigdb_retrieve("mouse", collections = c("h", "c2"),
                              id_type = "ensembl", rm_meta = FALSE)
  df.rat <- msigdb_retrieve("rat", collections = c("h", "c2"),
                            id_type = "ensembl", rm_meta = FALSE)

  # ensure we have the right type of identifiers.
  expect_true(all(grepl("ENSMUSG\\d+$", df.mouse$featureId)))
  expect_true(all(grepl("ENSRNOG\\d+$", df.rat$featureId)))

  # we should have the same number of gensets among these
  og.stats <- count(df.og, collection, name)
  mouse.stats <- count(df.mouse, collection, name)
  rat.stats <- count(df.rat, collection, name)
  expect_equal(mouse.stats$name, og.stats$name)
  expect_equal(rat.stats$name, og.stats$name)

  # the genesets should have approximately the same number of genes
  mouse.diff <- abs(og.stats$n - mouse.stats$n) / og.stats$n
  expect_true(mean(mouse.diff <= 0.15) > 0.95)

  rat.diff <- abs(og.stats$n - rat.stats$n) / og.stats$n
  expect_true(mean(rat.diff <= 0.18) > 0.95)

  # the gene names should largely be case-insenstive matches of each other
  unique.m <- distinct(df.mouse, human_symbol, .keep_all = TRUE)
  mouse.match <- tolower(unique.m$symbol) == tolower(unique.m$human_symbol)
  expect_true(mean(mouse.match) > 0.85)

  unique.r <- distinct(df.rat, human_symbol, .keep_all = TRUE)
  rat.match <- tolower(df.rat$symbol) == tolower(df.rat$human_symbol)
  expect_true(mean(rat.match) > 0.92)
})

test_that("ortholog maps with mismatched symbol names are legit", {
  # The randomness in our sample_n calls makes it hard to set a cutoff for the
  # ortholog matching, ie. sometimes
  #   mean(has.match$matched) > 0.96; and
  #   mean(cmp$human_symbol == cmp$HGNC.symbol) > 0.95
  # and other times they are in the low 90s
  # Let's set a seed here to be deterministic
  set.seed(0xDECADE)
  dfm <- msigdb_retrieve("mouse", collections = c("h", "c2"),
                         id_type = "ensembl", rm_meta = FALSE) %>%
    distinct(human_ensembl_id, featureId, .keep_all = TRUE) %>%
    select(human_ensembl_id, human_symbol, symbol, featureId)

  mismatched <- filter(dfm, tolower(dfm$symbol) != tolower(dfm$human_symbol))
  mm <- mismatched %>%
    distinct(human_ensembl_id, .keep_all = TRUE) %>%
    sample_n(100)

  bm <- loadNamespace("biomaRt")
  human <- bm$useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- bm$useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  xref <- bm$getLDS(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                    filters = "ensembl_gene_id",
                    values = mm$human_ensembl_id,
                    mart = human,
                    attributesL = c("ensembl_gene_id", "mgi_symbol"),
                    martL = mouse)

  # Since we are allowing for multimapping, we need to identify the rows in
  # dfm where the symbols don't match, then use the total mappings in dfm
  # to see if any of the orthologs mapped
  check.me <- dfm %>%
    filter(human_ensembl_id %in% mm$human_ensembl_id) %>%
    arrange(human_ensembl_id)

  cmp <- inner_join(xref, check.me, by = c("Gene.stable.ID" = "human_ensembl_id"))
  # human symbols are the same from query
  expect_true(mean(cmp$human_symbol == cmp$HGNC.symbol) > 0.97)

  # for the mouse map, we'll have to group this table by human ensembl id,
  # then look within those rows to see if the biomaRt retrieved ensembl IDs
  # match our mapped ids
  has.match <- cmp %>%
    group_by(Gene.stable.ID) %>%
    summarize(n = n(),
              matched = length(intersect(Gene.stable.ID.1, featureId)) > 0)
  expect_true(mean(has.match$matched) > 0.93)
})
