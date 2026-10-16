test_that("grouped_regression_data is shared by HLRM and CHLRM", {
  frame <- data.frame(
    response = c(1, 2, 3, 4, 5),
    x1 = c(-1, 0, 1, 2, 3),
    x2 = c(2, 1, 0, -1, -2),
    farm = c("b", "a", "b", "c", "a"))
  grouped <- grouped_regression_data(
    frame, "response", c("x1", "x2"), "farm")

  expect_s3_class(grouped, "vb_grouped_regression_data")
  expect_equal(grouped$mapping$group_label, c("b", "a", "c"))
  expect_equal(lengths(grouped$y), c(2, 2, 1))
  expect_equal(colnames(grouped$X[[1L]]),
               c("(Intercept)", "x1", "x2"))
  expect_equal(grouped$X[[1L]][, "x1"], c(-1, 1))
  expect_silent(.validate_fit_data(hlrm_model(3), grouped))
  expect_silent(.validate_fit_data(chlrm_model(2, 3), grouped))
})

test_that("grouped_regression_data validates roles and numeric content", {
  frame <- data.frame(y = 1:4, x = 2:5, g = c("a", "a", "b", "b"))
  no_intercept <- grouped_regression_data(frame, 1, 2, 3, intercept = FALSE)
  expect_equal(ncol(no_intercept$X[[1L]]), 1)
  expect_error(grouped_regression_data(frame, "y", "y", "g"), "distinct")
  bad <- frame; bad$x[1] <- NA_real_
  expect_error(grouped_regression_data(bad, "y", "x", "g"), "finite numeric")
  bad_group <- frame; bad_group$g[1] <- NA_character_
  expect_error(grouped_regression_data(
    bad_group, "y", "x", "g"), "missing")
})
