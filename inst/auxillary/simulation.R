make_new_obs <- function(d, i, hh_size, part_ids, enroll_per_hh, state_info) {
  if (d == 1) {
    data.frame(
      t = rep(d, hh_size[i]),
      part_id = part_ids[[i]],
      enroll = c(
        rep(1, enroll_per_hh[i]),
        rep(0, hh_size[i] - enroll_per_hh[i])
      ),
      state = sample(
        seq_along(state_info),
        hh_size[i],
        replace = TRUE,
        prob = state_info
      ),
      hh_size = rep(hh_size[i], hh_size[i]),
      hh_id = rep(i, hh_size[i])
    )
  } else {
    data.frame(
      t = rep(d, hh_size[i]),
      part_id = part_ids[[i]],
      enroll = c(
        rep(1, enroll_per_hh[i]),
        rep(0, hh_size[i] - enroll_per_hh[i])
      ),
      state = state_info,
      hh_size = rep(hh_size[i], hh_size[i]),
      hh_id = rep(i, hh_size[i])
    )
  }
}

#' @title Simulate SIR data
#' @description Simulate household data for a pathogen following an SIR model
#'
sim_sir <- function(
  eh_prob = 0.01,
  ih_prob = 0.05,
  n_hh = 100,
  hh_size = 1:5,
  tmax = 100,
  gamma = 1 / 5,
  covs_eh = c(0, 0),
  covs_ih = c(0, 0),
  obs_prob = list(c(0.05, 0.95, 0.05), c(0.01, 0.01, 0.8)),
  start_prob = c(1, 0, 0),
  complete_enroll = TRUE
) {
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH

  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))

  for (i in seq_along(covs_ih)) {
    x[, i] <- rbinom(sum(hh_size), 1, 0.4)
  }

  # Create participant IDs
  part_ids <- list()
  for (i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    if (hh_size[i] == 1) {
      enroll_per_hh[i] <- 1
    } else {
      enroll_per_hh[i] <- sample(1:hh_size[i], 1)
    }
  }

  # Infection state for all household members on all time steps
  complete_obs <- data.frame(
    t = numeric(),
    part_id = numeric(),
    enroll = numeric(),
    state = numeric(),
    hh_size = numeric(),
    hh_id = numeric()
  )

  last_x <- 0

  for (i in seq_along(hh_size)) {
    # Move HH members through SIR states
    for (d in 1:tmax) {
      if (d == 1) {
        new_obs <- make_new_obs(
          d,
          i,
          hh_size,
          part_ids,
          enroll_per_hh,
          start_prob
        )

        complete_obs <- complete_obs |>
          bind_rows(new_obs)

        prior <- new_obs$state
      } else {
        prior_inf <- sum(prior == 2)
        new_states <- rep(0, hh_size[i])
        for (part in 1:hh_size[i]) {
          if (prior[part] == 1) {
            eh_prob_x <- inv_logit(
              logit(eh_prob) +
                sum(x[last_x + part, ] * covs_eh)
            )
            ih_prob_x <- inv_logit(
              logit(ih_prob) +
                sum(x[last_x + part, ] * covs_ih)
            )
            no_inf_prob <- (1 - eh_prob_x) * (1 - ih_prob_x)^prior_inf
            new_states[part] <- sample(
              x = c(1, 2, 3),
              size = 1,
              prob = c(no_inf_prob, (1 - no_inf_prob), 0)
            )
          } else if (prior[part] == 2) {
            new_states[part] = sample(
              x = c(1, 2, 3),
              size = 1,
              prob = c(0, 1 - gamma, gamma)
            )
          } else {
            new_states[part] = sample(
              x = c(1, 2, 3),
              size = 1,
              prob = c(0, 0, 1)
            )
          }
        }

        new_obs <- make_new_obs(
          d,
          i,
          hh_size,
          part_ids,
          enroll_per_hh,
          new_states
        )

        complete_obs <- complete_obs |>
          bind_rows(new_obs)

        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }

  complete_obs <- complete_obs |>
    arrange(hh_id, t, part_id)

  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  outcome_names <- paste0("y", seq_along(obs_prob))
  colnames(outcome) <- outcome_names
  for (i in seq_len(nrow(complete_obs))) {
    for (j in seq_along(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i, j] <- sample(c(0, 1), 1, prob = c(1 - p1, p1))
    }
  }

  complete_obs <- complete_obs |>
    bind_cols(as.data.frame(outcome))

  if (!all(c(covs_eh, covs_ih) == 0)) {
    x <- as.data.frame(x)
    names(x) <- paste0("x", seq_len(ncol(x)))
  }

  if (!complete_enroll) {
    obs <- complete_obs[complete_obs$enroll == 1, ]
  } else {
    obs <- complete_obs
  }

  out <- list(obs = obs, complete_obs = complete_obs)

  if (!all(c(covs_eh, covs_ih) == 0)) {
    out <- append(out, list(x = x))
  }

  out
}


sim_siir <- function(
  eh_prob = 0.01,
  ih_prob = 0.05,
  n_hh = 100,
  hh_size = 1:5,
  tmax = 100,
  gamma = c(1 / 5, 1 / 30),
  split = c(0.7, 0.3),
  covs_eh = c(0, 0),
  covs_ih = c(0, 0),
  obs_prob = list(c(0.05, 0.95, 0.95, 0.05), c(0.01, 0.01, 0.01, 0.8)),
  start_prob = c(1, 0, 0, 0),
  complete_enroll = TRUE
) {
  if (length(ih_prob) == 1) {
    ih_prob <- rep(ih_prob, 2)
  }

  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH

  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))

  for (i in seq_along(covs_ih)) {
    x[, i] <- rbinom(sum(hh_size), 1, 0.4)
  }

  # Create participant IDs
  part_ids <- list()
  for (i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    if (hh_size[i] == 1) {
      enroll_per_hh[i] <- 1
    } else {
      enroll_per_hh[i] <- sample(1:hh_size[i], 1)
    }
  }

  # Infection state for all household members on all time steps
  complete_obs <- data.frame(
    t = numeric(),
    part_id = numeric(),
    enroll = numeric(),
    state = numeric(),
    hh_size = numeric(),
    hh_id = numeric()
  )

  last_x <- 0

  for (i in seq_along(hh_size)) {
    # Move HH members through SIR states
    for (d in 1:tmax) {
      if (d == 1) {
        new_obs <- make_new_obs(
          d,
          i,
          hh_size,
          part_ids,
          enroll_per_hh,
          start_prob
        )

        complete_obs <- complete_obs |>
          bind_rows(new_obs)

        prior <- new_obs$state
      } else {
        prior_inf <- c(sum(prior == 2), sum(prior == 3))
        new_states <- rep(0, hh_size[i])
        for (part in 1:hh_size[i]) {
          if (prior[part] == 1) {
            eh_prob_x <- inv_logit(
              logit(eh_prob) +
                sum(x[last_x + part, ] * covs_eh)
            )
            ih_prob_x <- inv_logit(
              logit(ih_prob) +
                sum(x[last_x + part, ] * covs_ih)
            )
            no_inf_prob <- (1 - eh_prob_x) * prod((1 - ih_prob_x)^prior_inf)
            new_states[part] <- sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(
                no_inf_prob,
                (1 - no_inf_prob) * split[1],
                (1 - no_inf_prob) * split[2],
                0
              )
            )
          } else if (prior[part] == 2) {
            new_states[part] = sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(0, 1 - gamma[1], 0, gamma[1])
            )
          } else if (prior[part] == 3) {
            new_states[part] = sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(0, 0, 1 - gamma[2], gamma[2])
            )
          } else {
            new_states[part] = sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(0, 0, 0, 1)
            )
          }
        }

        new_obs <- make_new_obs(
          d,
          i,
          hh_size,
          part_ids,
          enroll_per_hh,
          new_states
        )

        complete_obs <- complete_obs |>
          bind_rows(new_obs)

        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }

  complete_obs <- complete_obs |>
    arrange(hh_id, t, part_id)

  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  outcome_names <- paste0("y", seq_along(obs_prob))
  colnames(outcome) <- outcome_names
  for (i in seq_len(nrow(complete_obs))) {
    for (j in seq_along(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i, j] <- sample(c(0, 1), 1, prob = c(1 - p1, p1))
    }
  }

  complete_obs <- complete_obs |>
    bind_cols(as.data.frame(outcome))

  if (!complete_enroll) {
    obs <- complete_obs[complete_obs$enroll == 1, ]
  } else {
    obs <- complete_obs
  }

  list(obs = obs, complete_obs = complete_obs)
}

sim_seir <- function(
  eh_prob = 0.01,
  ih_prob = 0.05,
  n_hh = 100,
  hh_size = 1:5,
  tmax = 100,
  sigma = 1 / 2,
  gamma = 1 / 5,
  covs_eh = c(0, 0),
  covs_ih = c(0, 0),
  obs_prob = list(c(0.05, 0.05, 0.95, 0.05), c(0.01, 0.01, 0.1, 0.8)),
  start_prob = c(1, 0, 0, 0),
  complete_enroll = TRUE
) {
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH

  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))

  for (i in seq_along(covs_ih)) {
    x[, i] <- rbinom(sum(hh_size), 1, 0.4)
  }

  # Create participant IDs
  part_ids <- list()
  for (i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    if (hh_size[i] == 1) {
      enroll_per_hh[i] <- 1
    } else {
      enroll_per_hh[i] <- sample(1:hh_size[i], 1)
    }
  }

  # Infection state for all household members on all time steps
  complete_obs <- data.frame(
    t = numeric(),
    part_id = numeric(),
    enroll = numeric(),
    state = numeric(),
    hh_size = numeric(),
    hh_id = numeric()
  )

  last_x <- 0

  for (i in seq_along(hh_size)) {
    # Move HH members through SIR states
    for (d in 1:tmax) {
      if (d == 1) {
        new_obs <- make_new_obs(
          d,
          i,
          hh_size,
          part_ids,
          enroll_per_hh,
          start_prob
        )

        complete_obs <- complete_obs |>
          bind_rows(new_obs)

        prior <- new_obs$state
      } else {
        prior_inf <- sum(prior == 3)
        new_states <- rep(0, hh_size[i])
        for (part in 1:hh_size[i]) {
          if (prior[part] == 1) {
            eh_prob_x <- inv_logit(
              logit(eh_prob) +
                sum(x[last_x + part, ] * covs_eh)
            )
            ih_prob_x <- inv_logit(
              logit(ih_prob) +
                sum(x[last_x + part, ] * covs_ih)
            )
            no_inf_prob <- (1 - eh_prob_x) * (1 - ih_prob_x)^prior_inf
            new_states[part] <- sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(no_inf_prob, (1 - no_inf_prob), 0, 0)
            )
          } else if (prior[part] == 2) {
            new_states[part] = sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(0, 1 - sigma, sigma, 0)
            )
          } else if (prior[part] == 3) {
            new_states[part] = sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(0, 0, 1 - gamma, gamma)
            )
          } else {
            new_states[part] = sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(0, 0, 0, 1)
            )
          }
        }

        new_obs <- make_new_obs(
          d,
          i,
          hh_size,
          part_ids,
          enroll_per_hh,
          new_states
        )

        complete_obs <- complete_obs |>
          bind_rows(new_obs)

        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }

  complete_obs <- complete_obs |>
    arrange(hh_id, t, part_id)

  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  outcome_names <- paste0("y", seq_along(obs_prob))
  colnames(outcome) <- outcome_names
  for (i in seq_len(nrow(complete_obs))) {
    for (j in seq_along(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i, j] <- sample(c(0, 1), 1, prob = c(1 - p1, p1))
    }
  }

  complete_obs <- complete_obs |>
    bind_cols(as.data.frame(outcome))

  if (!complete_enroll) {
    obs <- complete_obs[complete_obs$enroll == 1, ]
  } else {
    obs <- complete_obs
  }

  list(obs = obs, complete_obs = complete_obs)
}
