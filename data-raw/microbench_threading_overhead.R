# Microbenchmark: the cost of compiling with threading, run single-threaded.
# ------------------------------------------------------------------------
# Answers two questions with numbers rather than theory:
#
#   (A) TLS penalty  - does compiling a model WITH stan_threads but running it
#       on ONE thread make it slower than the same model compiled WITHOUT
#       threading? (The "unthreaded model in a threaded package" cost, entirely
#       from the thread-local autodiff stack.)
#   (B) reduce_sum overhead - how much extra does the *threaded* model cost when
#       run on one thread (TLS + reduce_sum scheduling), and where does adding
#       threads turn it back into a win?
#
# Uses cmdstanr because it can toggle stan_threads per compile. The same .stan
# source (hmm.stan) is compiled two ways to isolate (A) cleanly; runs are paired
# by seed so each variant does identical sampling work.
#
# Run from the package root:  Rscript data-raw/microbench_threading_overhead.R

suppressMessages({
  library(cmdstanr)
  devtools::load_all(".", quiet = TRUE)
})

## ---- Tunables -----------------------------------------------------
n_hh          <- 10     # households (HMM fits are slow; keep modest)
iter_warmup   <- 150
iter_sampling <- 150
n_reps        <- 10     # replicates per config (paired by seed)
ctx_threads   <- min(4L, parallel::detectCores())  # for the "does it pay off" row

## ---- Data + init (shared by every fit) ---------------------------
dat <- sir[sir$hh_id <= n_hh, ]
inf_mod <- make_infection_model(
  transmit(from = "S", to = "I"),
  progress(from = "I", to = "R", gamma = NA)
)
obs_mod <- make_observation_model(
  pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
  igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
)
dat_stan <- make_stan_data(
  inf_mod, obs_mod, dat, init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10)
)
dat_stan <- lapply(dat_stan, function(x) if (is.data.frame(x)) as.matrix(x) else x)

init1 <- list(
  logit_params = array(rep(logit(0.5), dat_stan$n_params)),
  beta_eh      = logit(0.02),
  beta_ih      = array(rep(logit(0.02), dat_stan$n_inf_prob))
)
if (dat_stan$n_mult_params > 0) {
  init1$logit_mult_params <- array(rep(logit(0.5), dat_stan$n_mult_params))
}

## ---- Compile the three variants (each to its own dir) ------------
d <- function(x) { p <- file.path(tempdir(), x); dir.create(p, showWarnings = FALSE); p }
message("Compiling variants ...")
mod_plain <- cmdstan_model("inst/stan/hmm.stan", dir = d("plain"))
mod_tls   <- cmdstan_model("inst/stan/hmm.stan",
                           cpp_options = list(stan_threads = TRUE), dir = d("tls"))
mod_thr   <- cmdstan_model("inst/stan/hmm_threaded.stan",
                           cpp_options = list(stan_threads = TRUE), dir = d("thr"))

## ---- Timing helper (CmdStan-reported wall time, excludes compile) -
time_fit <- function(mod, seed, threads = NULL) {
  args <- list(
    data = dat_stan, init = list(init1), chains = 1,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed, refresh = 0, show_messages = FALSE
  )
  if (!is.null(threads)) args$threads_per_chain <- threads
  fit <- do.call(mod$sample, args)
  fit$time()$total
}

## ---- Run all configs, paired by seed -----------------------------
configs <- list(
  "hmm  (no stan_threads)"        = list(mod = mod_plain, threads = NULL),
  "hmm  (stan_threads, 1 thread)" = list(mod = mod_tls,   threads = 1L),
  "hmm_threaded (1 thread)"       = list(mod = mod_thr,   threads = 1L),
  "hmm_threaded (N threads)"      = list(mod = mod_thr,   threads = ctx_threads)
)
names(configs)[4] <- paste0("hmm_threaded (", ctx_threads, " threads)")

times <- matrix(NA_real_, nrow = length(configs), ncol = n_reps,
                dimnames = list(names(configs), paste0("rep", seq_len(n_reps))))
for (r in seq_len(n_reps)) {
  seed <- 100L + r
  message("rep ", r, " (seed ", seed, ") ...")
  for (i in seq_along(configs)) {
    times[i, r] <- time_fit(configs[[i]]$mod, seed, configs[[i]]$threads)
  }
}

## ---- Report ------------------------------------------------------
mn   <- apply(times, 1, min)                     # least-contended run (cleanest)
med  <- apply(times, 1, median)
pct_min <- 100 * (mn  - mn[1])  / mn[1]          # overhead from the min
pct_med <- 100 * (med - med[1]) / med[1]         # overhead from the median

res <- data.frame(
  config      = names(configs),
  min_s       = round(mn, 2),
  median_s    = round(med, 2),
  vs_base_min = sprintf("%+.1f%%", pct_min),
  vs_base_med = sprintf("%+.1f%%", pct_med),
  row.names   = NULL
)

cat("\n==================  THREADING OVERHEAD  ==================\n")
cat(sprintf("Data: %d households, %d obs | %d warmup + %d sampling | %d reps\n\n",
            n_hh, dat_stan$n_obs, iter_warmup, iter_sampling, n_reps))
cat("Per-rep CmdStan wall time (s):\n")
print(round(times, 2))
cat("\nOverhead vs plain hmm (min = cleanest, median = robust):\n")
print(res, row.names = FALSE)
cat(sprintf(
  "\n(A) TLS penalty  [hmm threaded-compile, 1 thread vs plain]: %+.1f%% (min) / %+.1f%% (median)\n",
  pct_min[2], pct_med[2]))
cat(sprintf(
  "(B) reduce_sum @1 thread [hmm_threaded 1 thread vs plain]:  %+.1f%% (min) / %+.1f%% (median)\n",
  pct_min[3], pct_med[3]))
cat(sprintf(
  "    threaded @%d threads vs plain:                          %+.1f%% (min) / %+.1f%% (median)\n",
  ctx_threads, pct_min[4], pct_med[4]))
cat("==========================================================\n")
