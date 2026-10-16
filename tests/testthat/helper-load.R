root <- Sys.getenv("VB_HLRM_ROOT", unset = "")
if (nzchar(root) && !"package:VariationalBayesHLRM" %in% search()) {
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Package 'pkgload' is required for source-tree validation")
  pkgload::load_all(root, quiet = TRUE, export_all = TRUE,
                    helpers = FALSE, attach_testthat = FALSE)
}
