#' Create a compact, non-resumable CHLRM MCMC reporting view
#'
#' Remove stored group-coefficient and group-variance histories from a new set
#' of backends and orchestration lists. CHLRM public reports use component
#' draws and allocations, which are retained in full. The original fit is not
#' modified. Plain numeric arrays may share memory until modified under R
#' copy-on-modify semantics; mutable backend environments are independent.
#'
#' @param fit A CHLRM MCMC `vb_model_fit` with plain serializable state.
#' @return A `vb_mcmc_reporting_view` inheriting from `vb_model_fit`. Public
#'   posterior summaries, draws, predictions, classifications, likelihoods and
#'   retained orchestration diagnostics remain available. `resume_mcmc()` and
#'   `save_mcmc_checkpoint()` reject the view. Keep the complete original fit
#'   or checkpoint for continuation and full native histories. Replacing the
#'   caller's fit can release memory; retaining it keeps its original histories
#'   alive. Construction itself is not a bound on whole-process RSS.
#' @export
mcmc_reporting_view <- function(fit) {
  if (!inherits(fit, "vb_model_fit") || !identical(fit$model, "CHLRM") ||
      !identical(fit$method, "MCMC"))
    stop("mcmc_reporting_view() requires a CHLRM MCMC vb_model_fit")
  if (inherits(fit, "vb_mcmc_reporting_view")) return(fit)
  checkpoint <- .typed_mcmc_checkpoint(fit, copy_state = identity)
  .checkpoint_plain_state(checkpoint)
  prune <- function(draws) {
    draws$Betaj <- NULL
    draws$sigma2j <- NULL
    draws
  }
  compact_backend <- function(state) {
    state$chains <- prune(state$chains)
    .restore_mcmc_backend_checkpoint(state, copy_state = identity)
  }
  out <- checkpoint$fit
  attr(out, "backend") <- compact_backend(checkpoint$summary_backend)
  backends <- lapply(checkpoint$chain_backends, compact_backend)
  if (length(backends) > 1L) {
    attr(out, "backends") <- backends
    orchestration <- checkpoint$orchestration
    if (!inherits(orchestration, "mcmc_multichain_fit"))
      stop("Multichain fit does not retain valid orchestration state")
    orchestration$chains <- lapply(orchestration$chains, prune)
    if (!is.null(orchestration$aligned_chains))
      orchestration$aligned_chains <- lapply(orchestration$aligned_chains, prune)
    attr(orchestration, "models") <- backends
    attr(out, "orchestration") <- orchestration
  }
  attr(out, "reporting_view") <- list(schema_version = "chlrm-mcmc-reporting-view-1",
    removed_native_fields = c("Betaj", "sigma2j"), resumable = FALSE)
  class(out) <- c("vb_mcmc_reporting_view", class(out))
  out
}

.require_complete_mcmc_fit <- function(fit) {
  if (inherits(fit, "vb_mcmc_reporting_view"))
    stop("A reporting view cannot be resumed or saved as a checkpoint; use the complete original fit or checkpoint")
  invisible(NULL)
}
