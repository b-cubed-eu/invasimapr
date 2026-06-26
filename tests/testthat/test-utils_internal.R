test_that("new_invasimapr_fit creates a classed container", {
  fit <- invasimapr:::new_invasimapr_fit(list(
    inputs = list(a = 1),
    meta   = list(n_sites = 3L, n_residents = 5L, n_invaders = 2L)
  ))
  expect_s3_class(fit, "invasimapr_fit")
  expect_true(is.list(fit))
})

test_that("print.invasimapr_fit prints a compact summary and returns invisibly", {
  fit <- invasimapr:::new_invasimapr_fit(list(
    inputs = list(),
    meta   = list(n_sites = 3L, n_residents = 5L, n_invaders = 2L)
  ))
  expect_output(print(fit), "<invasimapr_fit>")
  expect_output(print(fit), "sites: 3")
  expect_invisible(print(fit))
})
