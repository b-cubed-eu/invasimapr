test_that("prepare_inputs returns an invasimapr_fit with meta counts", {
  d <- utils::read.csv(
    system.file("extdata", "site_env_spp_simulated.csv.gz", package = "invasimapr")
  )
  names(d)[names(d) == "site_id"] <- "site"
  chr <- vapply(d, is.character, logical(1))
  d[chr] <- lapply(d[chr], as.factor)

  fit <- prepare_inputs(
    long_df      = d,
    site_col     = "site",
    env_prefix   = "^env",
    trait_prefix = "^trait",
    make_plots   = FALSE
  )

  expect_s3_class(fit, "invasimapr_fit")
  expect_true(!is.null(fit$inputs))
  expect_gt(fit$meta$n_sites, 0)
  expect_gt(fit$meta$n_residents, 0)
})
