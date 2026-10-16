root <- normalizePath(Sys.getenv("VB_HLRM_ROOT", unset = "."))
old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(root)
Sys.setenv(VB_HLRM_ROOT = root)

scripts <- c("thesis/final/code/runs/01-farms-hlrm/run.R",
             "thesis/final/code/runs/02-farms-chlrm/run.R",
             "thesis/final/code/runs/03-simulation/run.R")

status <- vector("list", length(scripts))
for (i in seq_along(scripts)) {
  message("Running ", scripts[i])
  started <- proc.time()[["elapsed"]]
  env <- new.env(parent = globalenv())
  err <- tryCatch({ sys.source(scripts[i], envir = env); NULL }, error = identity)
  status[[i]] <- data.frame(script = scripts[i],
                            seconds = proc.time()[["elapsed"]] - started,
                            ok = is.null(err),
                            error = if (is.null(err)) "" else conditionMessage(err))
  if (!is.null(err)) stop("Pipeline failed in ", scripts[i], ": ", conditionMessage(err))
  rm(env); gc(FALSE)
}

out_root <- file.path(root, "thesis", "final", "code")
run_output <- file.path(out_root, "run-output")
dir.create(run_output, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, status), file.path(run_output, "run-status.csv"), row.names = FALSE)
files <- list.files(file.path(out_root, "runs"), pattern = "^output$",
                    recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
files <- unlist(lapply(files, list.files, recursive = TRUE, full.names = TRUE),
                use.names = FALSE)
info <- file.info(files)
manifest <- data.frame(file = sub(paste0("^", out_root, "/?"), "", files),
                       bytes = info$size, modified = as.character(info$mtime),
                       stringsAsFactors = FALSE)
write.csv(manifest, file.path(run_output, "manifest.csv"), row.names = FALSE)
print(do.call(rbind, status))
