.history_test_fit <- local({
  cache <- new.env(parent = emptyenv())
  function(model = 'CHLRM', chains = 2L, p = 2L, K = 2L) {
    key <- paste(model, chains, p, K)
    if (exists(key, cache, inherits = FALSE)) return(cache[[key]])
    X <- if (p == 1L) matrix(1, 8, 1) else cbind(1, seq(-1, 1, length.out = 8))
    data <- list(X = rep(list(X), 4L), y = lapply(1:4, function(j) seq(-1, 1, length.out = 8) + j / 3))
    spec <- switch(model, LRM = lrm_model(p), HLRM = hlrm_model(p), CHLRM = chlrm_model(K, p))
    prior <- switch(model, LRM = lrm_prior(rep(0, p), diag(p), p+3, 1),
      HLRM = hlrm_prior(rep(0,p),diag(p),p+3,diag(p),4,2,2),
      CHLRM = chlrm_prior(rep(0,p),diag(p),p+3,diag(p),4,2,2,rep(1,K)))
    if (model == 'LRM') data <- list(X = X, y = data$y[[1L]])
    cache[[key]] <- fit_model(spec, data, prior, 'mcmc', initialization_from_prior(),
      mcmc_control(warmup=2, post_warmup_draws=8, chains=chains, evaluation_draws=2,
        stopping=if(chains>1) mcmc_stopping(min_post_warmup_draws=8, check_every=8,
          rhat_max=100,ess_bulk_min=.1,ess_tail_min=.1,
          relative_mcse_mean_max=1e10,relative_mcse_sd_max=1e10) else NULL,
        max_memory_mb=50),seed=72031,verbose=FALSE)
    cache[[key]]
  }
})
.history_test_create <- function(fit, directory) {
  checkpoint <- file.path(directory,'baseline.rds')
  save_mcmc_checkpoint(fit,checkpoint,file.path(directory,'baseline-store'),chunk_size=3L)
  create_mcmc_history(checkpoint,file.path(directory,'archive'),chunk_size=2L)
}
.history_test_resume <- function(fit,n) withCallingHandlers(resume_mcmc(fit,n),warning=function(w) {
  if(identical(trimws(conditionMessage(w)),'The ESS has been capped to avoid unstable estimates.')) invokeRestart('muffleWarning')
})
.history_test_subset <- function(draws, model, indices, fields=names(draws)) {
  out <- lapply(fields,function(name) .checkpoint_axis_subset(draws[[name]],indices,
    if(model=='CHLRM' && name=='z') 2L else 1L)); names(out)<-fields; out
}

test_that('immutable native histories equal complete continuation for all families and chains',{
  for(model in c('LRM','HLRM','CHLRM')) for(chains in c(1L,2L,4L)) {
    directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE),add=TRUE)
    fit<-.history_test_fit(model,chains); history<-.history_test_create(fit,directory)
    original<-history; state<-mcmc_history_restart_state(history)
    expect_identical(state,mcmc_restart_state(fit))
    base_files<-list.files(file.path(directory,'baseline-store'),recursive=TRUE,full.names=TRUE)
    base_hashes<-unname(tools::md5sum(base_files)); rng<-rng_snapshot()
    for(i in rev(seq_len(chains))) {
      a<-advance_mcmc_state(state,3L,chain=i);history<-append_mcmc_history(history,a);state<-a$state
      b<-advance_mcmc_state(state,2L,chain=i);history<-append_mcmc_history(history,b);state<-b$state
    }
    expect_identical(mcmc_history_restart_state(history),state)
    complete<-.history_test_resume(fit,5L)
    backends<-attr(complete,'backends',exact=TRUE);if(is.null(backends))backends<-list(.fit_backend(complete))
    base_backends<-attr(fit,'backends',exact=TRUE);if(is.null(base_backends))base_backends<-list(.fit_backend(fit))
    for(i in seq_len(chains)) {
      actual<-mcmc_history_draws(history,i,max_block_mb=10)
      expect_identical(actual$draws,backends[[i]]$chains)
      expect_identical(actual$draw_index$iteration,1:13)
      expect_identical(actual$draw_index$chain,rep.int(as.integer(i),13))
      selected<-c(1L,3L,7L,8L,9L,11L,12L,13L)
      fields<-rev(names(backends[[i]]$chains))
      expect_identical(mcmc_history_draws(history,i,selected,fields,max_block_mb=10)$draws,
        .history_test_subset(backends[[i]]$chains,model,selected,fields))
      expect_identical(mcmc_history_draws(original,i,max_block_mb=10)$draws,base_backends[[i]]$chains)
      expect_identical(state$chains[[i]]$rng_state,backends[[i]]$rng_state)
    }
    expect_identical(unname(tools::md5sum(base_files)),base_hashes)
    expect_identical(rng_snapshot(),rng)
    expect_identical(read_mcmc_history(history$directory,0L),original)
  }
})

test_that('derived CHLRM encoding reconstructs every group field without storing duplicate arrays',{
  directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE))
  history<-.history_test_create(.history_test_fit(),directory)
  block<-advance_mcmc_state(mcmc_history_restart_state(history),5L,chain=1L)
  history<-append_mcmc_history(history,block)
  x<-.history_manifest(history);record<-tail(x$blocks,1)[[1]]
  expect_identical(record$representation,'chlrm-derived-groups-1')
  store<-.history_store(history$directory);cache<-new.env(parent=emptyenv());decoded<-new.env(parent=emptyenv())
  expect_error(.history_child(record$root,'Betaj',store,cache,decoded),'Missing')
  expect_error(.history_child(record$root,'sigma2j',store,cache,decoded),'Missing')
  fields<-c('Betaj','sigma2j','z','Betak')
  expect_identical(mcmc_history_draws(history,1L,9:13,fields)$draws,block$draws[fields])
  for(model in c('LRM','HLRM','CHLRM')) {
    small<-file.path(directory,model);dir.create(small)
    fit<-.history_test_fit(model,1L,p=1L,K=1L);h<-.history_test_create(fit,small)
    block<-advance_mcmc_state(mcmc_history_restart_state(h),1L);h<-append_mcmc_history(h,block)
    expect_identical(mcmc_history_draws(h,1L,9L)$draws,block$draws)
  }
})

test_that('parent identities, chain frontiers and native consistency reject malformed appends',{
  directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE))
  history<-.history_test_create(.history_test_fit(),directory)
  state<-mcmc_history_restart_state(history); block<-advance_mcmc_state(state,3L,chain=1L)
  bad<-block;bad$parent_state_id<-'wrong';expect_error(append_mcmc_history(history,bad),'parent state')
  bad<-block;bad$draw_index$iteration[1]<-8L;expect_error(append_mcmc_history(history,bad),'gaps')
  bad<-block;bad$state$chains[[2]]$run_settings$samples<-9L;expect_error(append_mcmc_history(history,bad),'another chain')
  bad<-block;bad$draws$Beta[3,1]<-bad$draws$Beta[3,1]+1;expect_error(append_mcmc_history(history,bad),'Terminal native')
  bad<-block;bad$draws$Betaj[1,1,1]<-bad$draws$Betaj[1,1,1]+1;expect_error(append_mcmc_history(history,bad),'Derived group coefficient')
  bad<-block;bad$draws$z[1,1]<-999L;expect_error(append_mcmc_history(history,bad),'allocations')
  next_history<-append_mcmc_history(history,block)
  expect_error(append_mcmc_history(history,block),'Stale')
  expect_error(append_mcmc_history(next_history,block),'parent state')
  expect_false(dir.exists(file.path(history$directory,'.write-lock')))
  expect_identical(mcmc_history_restart_state(history),state)
  expect_equal(nrow(mcmc_history_draws(next_history,2L)$draw_index),8L)
  expect_error(mcmc_history_draws(history,1L,c(1,1)),'strictly increasing')
  expect_error(mcmc_history_draws(history,1L,9L),'in range')
  expect_error(mcmc_history_draws(history,3L),'Unknown')
  expect_error(mcmc_history_draws(history,1L,fields='unknown'),'unique native')
  expect_error(mcmc_history_draws(history,1L,max_output_mb=.00001),'max_output_mb')
  expect_error(mcmc_history_draws(history,1L,fields='Betaj',max_block_mb=.00001),'max_block_mb')
})

test_that('failed publication retains old snapshots and allows an explicit retry without resampling',{
  directory<-tempfile();dir.create(directory);on.exit(unlink(directory,recursive=TRUE))
  history<-.history_test_create(.history_test_fit(),directory)
  block<-advance_mcmc_state(mcmc_history_restart_state(history),3L,chain=1L)
  before<-mcmc_history_draws(history,1L)$draws
  original_publisher<-.publish_mcmc_history_snapshot
  local_mocked_bindings(.publish_mcmc_history_snapshot=function(...)stop('deliberate history publication failure'))
  expect_error(append_mcmc_history(history,block),'deliberate history')
  expect_identical(read_mcmc_history(history$directory),history)
  expect_identical(mcmc_history_draws(history,1L)$draws,before)
  expect_false(dir.exists(file.path(history$directory,'.write-lock')))
  local_mocked_bindings(.publish_mcmc_history_snapshot=original_publisher)
  next_history<-append_mcmc_history(history,block)
  expect_equal(next_history$generation,1L)
  expect_identical(mcmc_history_draws(next_history,1L,9:11)$draws,block$draws)
  lock<-file.path(history$directory,'.write-lock');dir.create(lock)
  expect_error(append_mcmc_history(next_history,advance_mcmc_state(block$state,1L,chain=1L)),'locked')
  expect_true(dir.exists(lock))
})

test_that('relocation and selected-chunk corruption preserve dependency and read boundaries',{
  directory<-tempfile();dir.create(directory);on.exit(unlink(c(directory,paste0(directory,'-moved')),recursive=TRUE))
  history<-.history_test_create(.history_test_fit(),directory)
  block<-advance_mcmc_state(mcmc_history_restart_state(history),3L,chain=1L);history<-append_mcmc_history(history,block)
  expected<-mcmc_history_draws(history,1L)$draws
  moved<-paste0(directory,'-moved');expect_true(file.rename(directory,moved))
  history<-read_mcmc_history(file.path(moved,'archive'))
  expect_identical(mcmc_history_draws(history,1L)$draws,expected)
  store<-.history_store(history$directory);x<-.history_manifest(history)
  cache<-new.env(parent=emptyenv());decoded<-new.env(parent=emptyenv())
  root<-x$blocks[[1]]$root
  key<-.history_child(root,'Beta',store,cache,decoded);node<-.checkpoint_store_get(key,store,cache)
  file<-file.path(store,'objects',paste0(node$chunks[1],'.rds'));writeBin(as.raw(1:4),file)
  expect_error(mcmc_history_draws(history,1L,9L,'Beta'),'checksum mismatch')
  expect_identical(mcmc_history_draws(history,1L,1L,'Beta')$draws,
    .history_test_subset(expected,'CHLRM',1L,'Beta'))
  # An unrelated new-block field still reads; no whole-history validation claim.
  expect_equal(nrow(mcmc_history_draws(history,1L,9L,'omega')$draw_index),1L)
  baseline<-file.path(moved,'baseline.rds');writeBin(as.raw(1:4),baseline)
  expect_error(mcmc_history_draws(history,1L,1L,'omega'),'Baseline checkpoint')
})
