test_that("return one set of chains per parameter", {
  # Subset to first ten households
  siir_sub <- sir[sir$hh_id <= 10, ]

  # Infection process model
  inf_mod <- make_infection_model(
    transmit(from = "S", to = c("Is", "Ia"), split = "phi"),
    progress(from = "Is", to = "R", gamma = NA),
    progress(from = "Ia", to = "R", gamma = NA)
  )

  obs_mod <- make_observation_model(
    pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.05, "R" = 0.05),
    igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8)
  )

  suppressWarnings(suppressMessages(
    model_out <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = siir_sub,
      init_probs = c(1 - 3 * 1e-10, 1e-10, 1e-10, 1e-10),
      stan_opts = stan_options(iter = 10)
    )
  ))

  # Variable is the third dimension in a draws array_object
  # Should be 4 variables: eh_prob, ih_prob, gamma, phi
  expect_equal(length(posterior::variables(model_out)), 4)
})

test_that("multiple infections probabilities supported", {
  # Subset to first ten households
  siir_sub <- sir[sir$hh_id <= 10, ]

  # Infection process model
  inf_mod <- make_infection_model(
    transmit(from = "S", to = c("Is", "Ia"), split = "phi"),
    progress(from = "Is", to = "R", gamma = NA),
    progress(from = "Ia", to = "R", gamma = NA),
    mult_inf_probs = TRUE
  )

  obs_mod <- make_observation_model(
    pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.05, "R" = 0.05),
    igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8)
  )

  suppressWarnings(suppressMessages(
    model_out <- run_model(
      inf_model = inf_mod,
      obs_model = obs_mod,
      data = siir_sub,
      init_probs = c(1 - 3 * 1e-10, 1e-10, 1e-10, 1e-10),
      stan_opts = stan_options(iter = 10)
    )
  ))

  # Variable is the third dimension in a draws array_object
  # Should be 5 variables: eh_prob, ih_prob_Is, ih_prob_Ia, gamma, phi
  expect_equal(length(posterior::variables(model_out)), 5)

  # Should be two variables that match ih_prob
  expect_equal(length(grep("ih_prob", posterior::variables(model_out))), 2)
})

test_that("entry validation rejects bad run_model inputs before sampling", {
  inf_mod <- make_infection_model(
    transmit(from = "S", to = "I"),
    progress(from = "I", to = "R", gamma = NA)
  )
  obs_mod <- make_observation_model(
    pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
    igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
  )
  good_data <- sir[sir$hh_id <= 10, ]
  good_init <- c(1 - 2 * 1e-10, 1e-10, 1e-10)

  # init_probs that do not sum to 1
  expect_error(
    run_model(inf_mod, obs_mod, good_data, init_probs = c(0.2, 0.2, 0.2)),
    "sum to 1"
  )

  # data missing a required structural column
  expect_error(
    run_model(
      inf_mod,
      obs_mod,
      good_data[, setdiff(names(good_data), "hh_size")],
      init_probs = good_init
    ),
    "hh_size"
  )

  # data missing an outcome column named in the observation model
  expect_error(
    run_model(
      inf_mod,
      obs_mod,
      good_data[, setdiff(names(good_data), "igg")],
      init_probs = good_init
    ),
    "igg"
  )

  # obs_model that is not a named list
  expect_error(
    run_model(inf_mod, list(), good_data, init_probs = good_init),
    "obs_model"
  )
})

# Integration test
test_that("State probabilities match expectation for toy examples", {
  # Some utility functions

  softmax <- function(x) {
    exp(x) / sum(exp(x))
  }

  get_diagonal_element <- function(m, i) {
    out = 1
    for (j in 1:nrow(m)) {
      if (j != i) {
        out = out - m[j, i]
      }
    }
    return(out)
  }

  replace_zeroes <- function(m, epsilon) {
    out <- matrix(NA, nrow = nrow(m), ncol = ncol(m))
    out <- m
    for (i in 1:nrow(m)) {
      for (j in 1:ncol(m)) {
        if (m[i, j] == 0) {
          out[i, j] <- epsilon
        }
      }
    }
    return(out)
  }

  normalize_cols <- function(m) {
    out <- matrix(NA, nrow = nrow(m), ncol = ncol(m))

    for (i in 1:ncol(m)) {
      out[, i] <- m[, i] / sum(m[, i])
    }
    return(out)
  }

  # R function to generate logalpha
  gen_la <- function(dat_stan, ih_prob, eh_prob, params, mult_params) {
    logalpha <- list()
    hh_size <- dat_stan$hh_size
    n_states <- dat_stan$n_states
    hh_tmax <- dat_stan$hh_tmax
    y <- dat_stan$y
    t_day <- dat_stan$t_day
    hh_start_ind <- dat_stan$hh_start_ind
    hh_end_ind <- dat_stan$hh_end_ind
    part_id <- dat_stan$part_id
    n_obs_type <- dat_stan$n_obs_type
    inf_states <- dat_stan$inf_states
    trans <- dat_stan$trans
    trans_index <- dat_stan$trans_index
    param_index <- dat_stan$param_index
    n_trans_fit <- dat_stan$n_trans_fit
    source_states <- dat_stan$source_states
    n_mult_fit <- dat_stan$n_mult_fit
    mult_index <- dat_stan$mult_index
    transition_multiplier <- dat_stan$transition_multiplier
    mult_param_index <- dat_stan$mult_param_index
    obs_prob <- dat_stan$obs_prob
    n_hh <- dat_stan$n_hh
    epsilon <- dat_stan$epsilon

    for (h in 1:n_hh) {
      alpha <- matrix(NA, nrow = hh_size[h] * n_states, ncol = max(hh_tmax))
      logalpha[[h]] <- matrix(
        NA,
        nrow = hh_size[h] * n_states,
        ncol = max(hh_tmax)
      )
      i_rows <- matrix(NA, nrow = hh_size[h], ncol = n_states) # rows in alpha corresponding to infectious states

      # subset to data only for the given HH
      y_hh = y[(hh_start_ind[h]):(hh_end_ind[h]), ]
      t_day_hh = t_day[(hh_start_ind[h]):(hh_end_ind[h])]
      part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])]

      index = 1

      {
        # START FORWARD ALGORITHM

        # fill first column of alpha using starting probabilities
        for (i in 1:hh_size[h]) {
          obs <- matrix(NA, nrow = n_obs_type, ncol = n_states) # observation component for enrolled memebrs, set to 1 if no observation for this time step

          ref = (n_states * (i - 1) + 1):(n_states * (i - 1) + n_states)

          obs_switch = 0

          if (t_day_hh[index] == 1) {
            if (part_id_hh[index] == i) {
              obs_switch = 1
            }
          }

          if (obs_switch == 1) {
            for (k in 1:n_obs_type) {
              if (y_hh[index, k] != -1) {
                obs[k, ] = obs_prob[k, y_hh[index, k], ]
              } else {
                obs[k, ] = 1
              }
            }
          } else {
            obs = matrix(1, nrow = n_obs_type, ncol = n_states)
          }

          if (obs_switch == 1) {
            index = min(index + 1, hh_end_ind[h] - hh_start_ind[h] + 1)
          }

          # Fill in starting probability for SIR states
          logalpha[[h]][ref, 1] = log(init_probs)
          for (k in 1:n_obs_type) {
            logalpha[[h]][ref, 1] = logalpha[[h]][ref, 1] + log(obs[k, ])
          }
          for (s in inf_states) {
            i_rows[i, s] = n_states * (i - 1) + s
          }

          # normalize and convert to the probability scale
          alpha[
            (n_states * (i - 1) + 1):(n_states * (i - 1) + n_states),
            1
          ] = softmax(logalpha[[h]][ref, 1])
        } # end participant loop - t=1, update logalpha with observation probability
        for (tt in 2:(hh_tmax[h])) {
          for (p in 1:hh_size[h]) {
            obs <- matrix(NA, nrow = n_obs_type, ncol = n_states)
            no_inf_prob <- numeric(n_states) # probability of avoiding all infections
            no_hh_inf_prob <- matrix(NA, nrow = hh_size[h], ncol = n_states) # probability of avoiding infection from each HH member

            ref = (n_states * (p - 1) + 1):(n_states * (p - 1) + n_states)

            logalpha_temp = logalpha[[h]][ref, tt - 1]

            obs_switch = 0

            if (t_day_hh[index] == tt) {
              if (part_id_hh[index] == p) {
                obs_switch = 1
              }
            }

            if (obs_switch == 1) {
              for (k in 1:n_obs_type) {
                if (y_hh[index, k] != -1) {
                  obs[k, ] = obs_prob[k, y_hh[index, k], ]
                } else {
                  obs[k, ] = 1
                }
              }
            } else {
              obs = matrix(1, nrow = n_obs_type, ncol = n_states)
            }

            if (obs_switch == 1) {
              index = min(index + 1, hh_end_ind[h] - hh_start_ind[h] + 1)
            }

            ct = 1
            for (s in 1:n_states) {
              if (s %in% inf_states) {
                no_hh_inf_prob[, s] = alpha[i_rows[, s], tt - 1] *
                  (1 - ih_prob[ct]) +
                  (1 - alpha[i_rows[, s], tt - 1]) # Pr of avoiding infection from each household member
                ct = ct + 1
                no_hh_inf_prob[p, s] = 1 # Particpant can't infect themselves
              } else {
                no_hh_inf_prob[, s] = 1
              }
              no_inf_prob[s] = prod(no_hh_inf_prob[, s]) # Probability of avoiding infection from all household members
            }

            # fill in tranistions that are being fit
            trans_temp = trans # rebuild from the base each step
            for (m in 1:n_trans_fit) {
              if (sum(source_states[m, ]) == 0) {
                trans_temp[
                  trans_index[m, 1],
                  trans_index[m, 2]
                ] = params[param_index[m]]
              } else {
                no_inf = 1
                for (s in 1:n_states) {
                  if (source_states[m, s] == 1) {
                    no_inf = no_inf * no_inf_prob[s]
                  }
                }
                trans_temp[trans_index[m, 1], trans_index[m, 2]] = 1 -
                  (no_inf * (1 - eh_prob))
              }
            }

            # fill in multipliers that are being fit
            mult_temp = transition_multiplier # need to reset mult_temp since it is self-referential
            if (n_mult_fit > 0) {
              for (m in 1:n_mult_fit) {
                if (mult_param_index[m] > 0) {
                  mult_temp[
                    mult_index[m, 1],
                    mult_index[m, 2]
                  ] = mult_params[mult_param_index[m]]
                } else {
                  mult_temp[mult_index[m, 1], mult_index[m, 2]] = mult_temp[
                    mult_index[m, 1],
                    mult_index[m, 2]
                  ] -
                    mult_params[-mult_param_index[m]]
                }
              }
            }

            # transition splits
            trans_temp <- trans_temp * mult_temp

            # fill in diagonals (columns must sum to one)
            for (i in 1:ncol(trans_temp)) {
              trans_temp[i, i] = get_diagonal_element(trans_temp, i)
            }

            # replace zeroes with epsilon and normalize
            trans_temp = replace_zeroes(trans_temp, epsilon)
            trans_temp = normalize_cols(trans_temp)

            # Compute the probability of each epidemiological state
            logalpha[[h]][ref, tt] = log(trans_temp %*% exp(logalpha_temp))
            for (k in 1:n_obs_type) {
              logalpha[[h]][ref, tt] = logalpha[[h]][ref, tt] + log(obs[k, ])
            }

            # normalize and convert to probability scale
            alpha[
              (n_states * (p - 1) + 1):(n_states * (p - 1) + n_states),
              tt
            ] = softmax(logalpha[[h]][ref, tt])
          } # end participant loop - update logalpha with observation probability
        } # end time loop
      } # END FORWARD ALGORITHM
    }

    return(do.call(rbind, logalpha))
  }

  # Infection process model
  inf_process <- make_infection_model(
    transmit(from = "S", to = "I"),
    progress(from = "I", to = "R", gamma = NA)
  )

  # Observation process model
  obs_process <- make_observation_model(
    pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
    igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
  )

  init_probs <- c(1 - 2 * 1e-10, 1e-10, 1e-10)

  ######## 2 HH, all observations for first 10 days ###################
  sir_sub <- sir %>%
    filter(t <= 10, hh_id <= 2)

  np <- nrow(unique(sir_sub %>% dplyr::select(hh_id, part_id)))

  dat_stan <- hestia:::make_stan_data(
    inf_model = inf_process,
    obs_model = obs_process,
    data = sir_sub,
    init_probs = init_probs,
    epsilon = 1e-10,
    ih_cov = NULL,
    eh_cov = NULL
  )

  # Run hestia
  test_res <- run_model(
    inf_model = inf_process,
    obs_model = obs_process,
    data = sir_sub,
    init_probs = init_probs,
    stan_opts = stan_options(iter = 50),
    save_states = TRUE
  )

  # Select an iteration to check
  ch <- 3
  it <- 13

  # extract all logalpha variables at iteration i, chain c
  vals <- test_res[it, ch, grepl("^logalpha\\[", variables(test_res))]
  # reshape to 24 x 10
  logalpha_mat <- matrix(vals, nrow = np * 3, ncol = 10, byrow = FALSE)

  ih_prob <- as.numeric(test_res[it, ch, "ih_prob"])
  eh_prob <- as.numeric(test_res[it, ch, "eh_prob"])
  params <- c(test_res[it, ch, "gamma"])
  mult_params <- numeric(0)

  test_la <- gen_la(dat_stan, ih_prob, eh_prob, params, mult_params)

  # Test that values are equal
  expect_equal(logalpha_mat, test_la)

  ######## 2 HH, observations for days 6-10 ###################
  sir_sub <- sir %>%
    filter(t <= 10, t > 5, hh_id <= 2)

  np <- nrow(unique(sir_sub %>% dplyr::select(hh_id, part_id)))

  dat_stan <- hestia:::make_stan_data(
    inf_model = inf_process,
    obs_model = obs_process,
    data = sir_sub,
    init_probs = init_probs,
    epsilon = 1e-10,
    ih_cov = NULL,
    eh_cov = NULL
  )

  # Run hestia
  test_res <- run_model(
    inf_model = inf_process,
    obs_model = obs_process,
    data = sir_sub,
    init_probs = init_probs,
    stan_opts = stan_options(iter = 50),
    save_states = TRUE
  )

  # Select an iteration to check
  ch <- 3
  it <- 13

  # extract all logalpha variables at iteration i, chain c
  vals <- test_res[it, ch, grepl("^logalpha\\[", variables(test_res))]
  # reshape to 24 x 10
  logalpha_mat <- matrix(vals, nrow = np * 3, ncol = 10, byrow = FALSE)

  ih_prob <- as.numeric(test_res[it, ch, "ih_prob"])
  eh_prob <- as.numeric(test_res[it, ch, "eh_prob"])
  params <- c(test_res[it, ch, "gamma"])
  mult_params <- numeric(0)

  test_la <- gen_la(dat_stan, ih_prob, eh_prob, params, mult_params)

  # Test that values are equal
  expect_equal(logalpha_mat, test_la)

  ######## 2 HH, observations for first 10 days, SIIR ###################

  inf_process <- make_infection_model(
    transmit(
      from = "S",
      to = c("Is", "Ia"),
      source = c("Is", "Ia"),
      split = "phi"
    ),
    progress(from = "Is", to = "R", gamma_s = NA),
    progress(from = "Ia", to = "R", gamma_a = NA),
    mult_inf_probs = TRUE
  )

  obs_process <- make_observation_model(
    pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.95, "R" = 0.05),
    igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8),
    symp = c("S" = 0.03, "Is" = 1 - 1e-10, "Ia" = 0.03, "R" = 0.03)
  )

  siir_sub <- siir %>%
    filter(hh_id <= 2, t <= 10)

  init_probs <- c(1 - 3 * 1e-10, 1e-10, 1e-10, 1e-10)

  np <- nrow(unique(siir_sub %>% dplyr::select(hh_id, part_id)))

  dat_stan <- hestia:::make_stan_data(
    inf_model = inf_process,
    obs_model = obs_process,
    data = siir_sub,
    init_probs = init_probs,
    epsilon = 1e-10,
    ih_cov = NULL,
    eh_cov = NULL
  )

  test_res <- run_model(
    inf_model = inf_process,
    obs_model = obs_process,
    data = siir_sub,
    init_probs = init_probs,
    stan_opts = stan_options(iter = 50),
    save_states = TRUE
  )

  # Select an iteration to check
  ch <- 3
  it <- 13

  # extract all logalpha variables at iteration i, chain c
  vals <- test_res[it, ch, grepl("^logalpha\\[", variables(test_res))]
  # reshape to 24 x 10
  logalpha_mat <- matrix(vals, nrow = np * 4, ncol = 10, byrow = FALSE)

  ih_prob <- c(
    as.numeric(test_res[it, ch, "ih_prob_Is"]),
    as.numeric(test_res[it, ch, "ih_prob_Ia"])
  )
  eh_prob <- as.numeric(test_res[it, ch, "eh_prob"])
  params <- c(
    as.numeric(test_res[it, ch, "gamma_s"]),
    as.numeric(test_res[it, ch, "gamma_a"])
  )
  mult_params <- c(test_res[it, ch, "phi"])

  test_la <- gen_la(dat_stan, ih_prob, eh_prob, params, mult_params)

  # Test that values are equal
  expect_equal(logalpha_mat, test_la)
})
