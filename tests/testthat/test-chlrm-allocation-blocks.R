test_that("blocked allocation marginals equal full-history summaries exactly", {
  for (storage in c("integer", "double")) for (m in c(1L, 7L))
    for (S in c(1L, 17L, 1031L)) {
      z <- matrix(rep(c(1L, 3L, 2L, 1L, 4L), length.out = m * S), m, S,
        dimnames = list(paste0("group", seq_len(m)), paste0("draw", seq_len(S))))
      storage.mode(z) <- storage
      expected <- matrix(0, m, 5L)
      for (k in 1:5) expected[, k] <- rowMeans(z == k)
      occupied <- apply(z, 2L, function(x) length(unique(x)))
      for (size in c(1L, 16L, 1024L, 2048L)) {
        actual <- .chlrm_allocation_marginals(z, 5L, size)
        expect_identical(actual$probabilities, expected)
        expect_identical(actual$occupied, occupied)
      }
    }
})

test_that("bounded label validation checks boundaries and keeps valid input unchanged", {
  z <- matrix(rep(1:3, length.out = 70000L), 10L, 7000L)
  before <- z
  expect_true(.valid_chlrm_allocation_labels(z, ncol(z), 3L))
  expect_identical(z, before)
  for (index in c(1L, 65536L, 65537L, length(z)))
    for (value in c(NA_real_, NaN, Inf, -Inf, 0, 4, 1.5)) {
      bad <- z; bad[index] <- value
      expect_false(.valid_chlrm_allocation_labels(bad, ncol(z), 3L))
    }
  expect_false(.valid_chlrm_allocation_labels(z, ncol(z) + 1L, 3L))
  expect_false(.valid_chlrm_allocation_labels(as.vector(z), ncol(z), 3L))
})
