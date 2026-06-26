test_that("assemble_matrices builds aligned core objects from a long table", {
  d <- utils::read.csv(
    system.file("extdata", "site_env_spp_simulated.csv.gz", package = "invasimapr")
  )
  names(d)[names(d) == "site_id"] <- "site"
  chr <- vapply(d, is.character, logical(1))
  d[chr] <- lapply(d[chr], as.factor)

  core <- assemble_matrices(long_df = d, make_plots = FALSE)

  expect_type(core, "list")
  expect_true(all(c("site_df", "env_df", "comm_res", "traits_res") %in% names(core)))
  # community columns (residents) must match trait rows
  expect_setequal(colnames(core$comm_res), rownames(core$traits_res))
  # site rows must be aligned across the core matrices
  expect_identical(rownames(core$site_df), rownames(core$comm_res))
})
