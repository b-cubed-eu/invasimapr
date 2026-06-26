test_that("simulate_invaders returns the requested number of invaders", {
  traits <- utils::read.csv(
    system.file("extdata", "species_traits.csv", package = "invasimapr")
  )
  set.seed(1)
  inv <- simulate_invaders(
    resident_traits = traits,
    n_inv           = 5,
    species_col     = "species",
    mode            = "columnwise"
  )
  expect_s3_class(inv, "data.frame")
  expect_equal(nrow(inv), 5)
})

test_that("simulate_invaders errors on an empty resident table", {
  empty <- data.frame(species = character(0), trait_cont1 = numeric(0))
  expect_error(simulate_invaders(empty, n_inv = 3))
})
