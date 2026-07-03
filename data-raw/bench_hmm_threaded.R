# Benchmark + equivalence check: hmm.stan vs hmm_threaded.stan
# ------------------------------------------------------------------
# Fits BOTH models on the same small SIR subset via cmdstanr, times each
# (wall clock, sampling only), and compares the posterior means of the shared
# parameters. The two are the same statistical model, so their posteriors
# should agree within Monte Carlo error; the threaded one just parallelises the
# forward algorithm over households with reduce_sum.
#
# Run from the package root:  Rscript data-raw/bench_hmm_threaded.R
# Requires: cmdstanr + a working CmdStan toolchain, and devtools.

suppressMessages({
  library(cmdstanr)
  library(posterior)
  devtools::load_all(".", quiet = TRUE)   # exposes make_stan_data(), logit(), sir
})

## ---- Tunables -----------------------------------------------------
n_hh          <- 20     # households (keep small; increase to see threading pay off)
n_chains      <- 2
n_threads     <- 4      # threads_per_chain for the threaded model
iter_warmup   <- 500
iter_sampling <- 500
seed          <- 123

## ---- 1. Small dummy dataset --------------------------------------
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
  inf_mod, obs_mod, dat,
  init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10)
)
# cmdstanr needs matrices, not data frames.
dat_stan <- lapply(dat_stan, function(x) if (is.data.frame(x)) as.matrix(x) else x)

## ---- 2. Init (mirrors run_model's defaults) ----------------------
make_init <- function() {
  out <- list(
    logit_params = array(rep(logit(0.5), dat_stan$n_params)),
    beta_eh      = logit(0.02),
    beta_ih      = array(rep(logit(0.02), dat_stan$n_inf_prob))
  )
  if (dat_stan$n_mult_params > 0) {
    out$logit_mult_params <- array(rep(logit(0.5), dat_stan$n_mult_params))
  }
  out
}
init <- rep(list(make_init()), n_chains)

## ---- 3. Compile both (outside the timed region) ------------------
message("Compiling hmm.stan ...")
mod_serial <- cmdstan_model("inst/stan/hmm.stan", dir = tempdir())
message("Compiling hmm_threaded.stan (stan_threads = TRUE) ...")
mod_thread <- cmdstan_model(
  "inst/stan/hmm_threaded.stan",
  cpp_options = list(stan_threads = TRUE),
  dir = tempdir()
)

## ---- 4. Fit both, identical settings -----------------------------
common <- list(
  data = dat_stan, init = init, chains = n_chains,
  parallel_chains = n_chains, iter_warmup = iter_warmup,
  iter_sampling = iter_sampling, seed = seed, refresh = 0,
  show_messages = FALSE
)

message("Fitting serial (hmm.stan) ...")
t_serial <- system.time(
  fit_serial <- do.call(mod_serial$sample, common)
)[["elapsed"]]

message("Fitting threaded (hmm_threaded.stan, ", n_threads, " threads/chain) ...")
t_thread <- system.time(
  fit_thread <- do.call(mod_thread$sample,
                        c(common, list(threads_per_chain = n_threads)))
)[["elapsed"]]

## ---- 5. Compare posteriors ---------------------------------------
pars <- c("eh_prob", "ih_prob", "params")
if (dat_stan$n_mult_params > 0) pars <- c(pars, "mult_params")

s1 <- summarise_draws(fit_serial$draws(pars), mean, sd, mcse_mean)
s2 <- summarise_draws(fit_thread$draws(pars), mean, sd, mcse_mean)

cmp <- data.frame(
  variable    = s1$variable,
  mean_serial = round(s1$mean, 5),
  mean_thread = round(s2$mean, 5),
  abs_diff    = abs(s1$mean - s2$mean),
  mcse_comb   = sqrt(s1$mcse_mean^2 + s2$mcse_mean^2)
)
cmp$within_3mcse <- cmp$abs_diff < 3 * cmp$mcse_comb

cmp$abs_diff /cmp$mcse_comb

## ---- 6. Report ---------------------------------------------------
cat("\n====================  RESULTS  ====================\n")
cat(sprintf("Data: %d households, %d observations\n", n_hh, dat_stan$n_obs))
cat(sprintf("Sampling: %d chains x (%d warmup + %d sampling)\n\n",
            n_chains, iter_warmup, iter_sampling))

cat("Wall clock (sampling only):\n")
cat(sprintf("  hmm.stan       (serial)      : %6.1f s\n", t_serial))
cat(sprintf("  hmm_threaded   (%d thr/chain) : %6.1f s   (%.2fx)\n\n",
            n_threads, t_thread, t_serial / t_thread))

cat("Posterior means (should agree within Monte Carlo error):\n")
print(cmp, row.names = FALSE, digits = 4)

ok <- all(cmp$within_3mcse)
cat(sprintf(
  "\nEquivalence: %s  (all means within 3x combined MCSE: %s)\n",
  if (ok) "PASS" else "CHECK", ok
))
cat(sprintf("Max |difference in means|: %.3g\n", max(cmp$abs_diff)))
cat("===================================================\n")
