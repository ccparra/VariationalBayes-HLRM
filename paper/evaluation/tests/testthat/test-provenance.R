test_that("configuration fingerprints are stable and sensitive to content", {
  first <- list(beta = data.frame(x = 1:2, y = c("a", "b")), alpha = 3)
  reordered <- list(alpha = 3, beta = data.frame(x = 1:2, y = c("a", "b")))
  changed <- list(alpha = 4, beta = data.frame(x = 1:2, y = c("a", "b")))
  expect_identical(.paper_fingerprint(first), .paper_fingerprint(reordered))
  expect_false(identical(.paper_fingerprint(first), .paper_fingerprint(changed)))
})

test_that("run provenance fingerprints package, suite and configuration", {
  configuration <- list(profile = "smoke", seed = 17L)
  first <- paper_run_provenance("TEST", "smoke", configuration, root)
  second <- paper_run_provenance("TEST", "smoke", configuration, root)
  expect_s3_class(first, "paper_run_provenance")
  expect_identical(first$schema_version, "paper-run-provenance-1")
  expect_identical(first$configuration_fingerprint,
                   second$configuration_fingerprint)
  expect_identical(first$source_fingerprint, second$source_fingerprint)
  expect_identical(first$suite_source_fingerprint,
                   second$suite_source_fingerprint)
  expect_identical(first$run_fingerprint, second$run_fingerprint)
  expect_identical(first$output_fingerprint, second$output_fingerprint)
  expect_match(first$run_fingerprint, "^md5:[0-9a-f]{32}$")
})

test_that("provenance is written and controls cache compatibility", {
  provenance <- paper_run_provenance(
    "TEST", "smoke", list(profile = "smoke"), root)
  output <- tempfile("paper-provenance-output-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  written <- write_paper_provenance(provenance, output)
  expect_true(file.exists(file.path(output, "provenance.csv")))
  expect_true(file.exists(file.path(output, "provenance.rds")))
  expect_identical(written$run_fingerprint, provenance$run_fingerprint)

  result <- structure(list(
    schema_version = "paper-test-result-1",
    task = data.frame(task_id = "task-1"), provenance = provenance),
    class = c("paper_test_result", "list"))
  expect_true(.paper_result_compatible(
    result, "paper_test_result", "paper-test-result-1", "task-1",
    provenance))
  changed <- paper_run_provenance(
    "TEST", "smoke", list(profile = "smoke", changed = TRUE), root)
  expect_false(.paper_result_compatible(
    result, "paper_test_result", "paper-test-result-1", "task-1", changed))
})
