context("GO slim")

test_that("go_slim = 'generic' returns slim ontology", {
  c5.all <- msigdb_retrieve("C5", go_slim = FALSE)
  c5.slim <- msigdb_retrieve("C5", go_slim = TRUE)

  sum.all <- c5.all %>%
    group_by(name) %>%
    summarize(gs_id = gs_id[1], subcategory = subcategory[1], n = n()) %>%
    ungroup() %>%
    arrange(subcategory, name)
  sum.slim <- c5.slim %>%
    group_by(name) %>%
    summarize(gs_id = gs_id[1], subcategory = subcategory[1], n = n()) %>%
    ungroup() %>%
    arrange(subcategory, name)

  expect_lt(nrow(sum.slim), nrow(sum.all))
  expect_setequal(sum.slim[["subcategory"]], sum.all[["subcategory"]])
})
