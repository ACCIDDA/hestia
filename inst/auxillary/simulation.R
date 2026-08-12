#' @title Build one time step of household observations
#'
#' @description Helper used by the simulation functions to assemble the rows for
#'   a single household on a single time step. On the first time step states are
#'   drawn from the starting state probabilities; on later steps the already
#'   computed states are carried in directly.
#'
#' @param d Integer time step. When `d == 1` states are sampled from
#'   `state_info`, otherwise `state_info` is treated as the states themselves.
#' @param i Integer index of the household being built.
#' @param hh_size Integer vector of household sizes, indexed by household.
#' @param part_ids List of participant ID vectors, one element per household.
#' @param enroll_per_hh Integer vector giving the number of enrolled
#'   participants in each household.
#' @param state_info On the first time step, a vector of starting state
#'   probabilities. On later time steps, the vector of member states for the
#'   household.
#'
#' @return A data frame with one row per household member and the columns `t`,
#'   `part_id`, `enroll`, `state`, `hh_size`, and `hh_id`.
#'
#' @keywords internal
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

#' @title Simulate SIR household data
#'
#' @description Simulate longitudinal household outcome data for a pathogen that
#'   follows a susceptible-infected-recovered (SIR) process. Members move
#'   through three states (1 = susceptible, 2 = infected, 3 = recovered) and
#'   each true state generates one or more binary observations through an
#'   observation process.
#'
#' @details Each household is simulated independently for `tmax` time steps.
#'   Susceptible members can be infected from outside the household (extra
#'   household) or by infected household members (intra household); the
#'   probabilities can be shifted by covariates through `covs_eh` and `covs_ih`.
#'   Infected members recover at rate `gamma`. For each true state, every
#'   observation type in `obs_prob` produces a positive result with the listed
#'   probability.
#'
#' @param eh_prob Baseline per step probability that a susceptible member is
#'   infected from outside the household.
#' @param ih_prob Baseline per step probability that a susceptible member is
#'   infected by a single infected household member.
#' @param n_hh Number of households to simulate.
#' @param hh_size Vector of possible household sizes; each household size is
#'   drawn from this vector with replacement.
#' @param tmax Number of time steps to simulate per household.
#' @param gamma Per step probability that an infected member recovers.
#' @param covs_eh Numeric vector of coefficients for the extra household
#'   covariates, on the logit scale.
#' @param covs_ih Numeric vector of coefficients for the intra household
#'   covariates, on the logit scale.
#' @param obs_prob List of observation types. Each element is a vector giving
#'   the probability of a positive observation for that type in each true state.
#' @param start_prob Vector of starting state probabilities used to draw each
#'   member's state on the first time step.
#' @param complete_enroll If `TRUE`, return observations for every household
#'   member. If `FALSE`, return only enrolled members in `obs`.
#'
#' @return A list with `obs` (the observed rows, either all members or only
#'   enrolled members depending on `complete_enroll`) and `complete_obs` (every
#'   member at every time step). When any covariate coefficient is non zero, a
#'   third element `x` holds the simulated covariate matrix. Each observation
#'   data frame has columns `t`, `part_id`, `enroll`, `state`, `hh_size`,
#'   `hh_id`, and one `y` column per observation type.
#'
#' @examples
#' # Simulate a basic SIR study with two observation types
#' sim_sir(
#'   n_hh = 20,
#'   tmax = 50,
#'   obs_prob = list(c(0.05, 0.95, 0.05), c(0.01, 0.01, 0.8))
#' )
#'
#' @keywords internal
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
            eh_prob_x <- plogis(
              qlogis(eh_prob) +
                sum(x[last_x + part, ] * covs_eh)
            )
            ih_prob_x <- plogis(
              qlogis(ih_prob) +
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


#' @title Simulate SIIR household data
#'
#' @description Simulate longitudinal household outcome data for a pathogen with
#'   two parallel infectious compartments, for example symptomatic and
#'   asymptomatic infection. Members move through four states (1 = susceptible,
#'   2 = first infectious compartment, 3 = second infectious compartment,
#'   4 = recovered).
#'
#' @details On infection a susceptible member enters the first or second
#'   infectious compartment according to `split`. The two infectious
#'   compartments can have separate intra household infection probabilities and
#'   recovery rates. As in [sim_sir()], each true state generates binary
#'   observations through `obs_prob`.
#'
#' @param eh_prob Baseline per step probability that a susceptible member is
#'   infected from outside the household.
#' @param ih_prob Baseline per step intra household infection probability. May
#'   be a single value (shared across both infectious compartments) or a length
#'   two vector (one per infectious compartment).
#' @param n_hh Number of households to simulate.
#' @param hh_size Vector of possible household sizes; each household size is
#'   drawn from this vector with replacement.
#' @param tmax Number of time steps to simulate per household.
#' @param gamma Length two vector of per step recovery probabilities, one for
#'   each infectious compartment.
#' @param split Length two vector giving the proportion of new infections that
#'   enter each infectious compartment.
#' @param covs_eh Numeric vector of coefficients for the extra household
#'   covariates, on the logit scale.
#' @param covs_ih Numeric vector of coefficients for the intra household
#'   covariates, on the logit scale.
#' @param obs_prob List of observation types. Each element is a vector giving
#'   the probability of a positive observation for that type in each true state.
#' @param start_prob Vector of starting state probabilities used to draw each
#'   member's state on the first time step.
#' @param complete_enroll If `TRUE`, return observations for every household
#'   member. If `FALSE`, return only enrolled members in `obs`.
#'
#' @return A list with `obs` (the observed rows) and `complete_obs` (every
#'   member at every time step). Each observation data frame has columns `t`,
#'   `part_id`, `enroll`, `state`, `hh_size`, `hh_id`, and one `y` column per
#'   observation type.
#'
#' @examples
#' # Simulate a SIIR study splitting infections 70/30 across two compartments
#' sim_siir(
#'   n_hh = 20,
#'   tmax = 50,
#'   ih_prob = c(0.05, 0.025),
#'   gamma = c(1 / 5, 1 / 3),
#'   split = c(0.7, 0.3)
#' )
#'
#' @keywords internal
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
            eh_prob_x <- plogis(
              qlogis(eh_prob) +
                sum(x[last_x + part, ] * covs_eh)
            )
            ih_prob_x <- plogis(
              qlogis(ih_prob) +
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

#' @title Simulate SEIR household data
#'
#' @description Simulate longitudinal household outcome data for a pathogen that
#'   follows a susceptible-exposed-infected-recovered (SEIR) process. Members
#'   move through four states (1 = susceptible, 2 = exposed, 3 = infected,
#'   4 = recovered).
#'
#' @details Susceptible members can be infected from outside or inside the
#'   household and first enter the exposed state. Exposed members become
#'   infected at rate `sigma`, infected members recover at rate `gamma`, and
#'   only infected members drive intra household transmission. As in
#'   [sim_sir()], each true state generates binary observations through
#'   `obs_prob`.
#'
#' @param eh_prob Baseline per step probability that a susceptible member is
#'   infected from outside the household.
#' @param ih_prob Baseline per step probability that a susceptible member is
#'   infected by a single infected household member.
#' @param n_hh Number of households to simulate.
#' @param hh_size Vector of possible household sizes; each household size is
#'   drawn from this vector with replacement.
#' @param tmax Number of time steps to simulate per household.
#' @param sigma Per step probability that an exposed member becomes infected.
#' @param gamma Per step probability that an infected member recovers.
#' @param covs_eh Numeric vector of coefficients for the extra household
#'   covariates, on the logit scale.
#' @param covs_ih Numeric vector of coefficients for the intra household
#'   covariates, on the logit scale.
#' @param obs_prob List of observation types. Each element is a vector giving
#'   the probability of a positive observation for that type in each true state.
#' @param start_prob Vector of starting state probabilities used to draw each
#'   member's state on the first time step.
#' @param complete_enroll If `TRUE`, return observations for every household
#'   member. If `FALSE`, return only enrolled members in `obs`.
#'
#' @return A list with `obs` (the observed rows) and `complete_obs` (every
#'   member at every time step). Each observation data frame has columns `t`,
#'   `part_id`, `enroll`, `state`, `hh_size`, `hh_id`, and one `y` column per
#'   observation type.
#'
#' @examples
#' # Simulate a SEIR study with an exposed period before infectiousness
#' sim_seir(
#'   n_hh = 20,
#'   tmax = 50,
#'   sigma = 1 / 2,
#'   gamma = 1 / 5
#' )
#'
#' @keywords internal
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
            eh_prob_x <- plogis(
              qlogis(eh_prob) +
                sum(x[last_x + part, ] * covs_eh)
            )
            ih_prob_x <- plogis(
              qlogis(ih_prob) +
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


#' @title Simulate SIIR household data with competing infections
#'
#' @description Simulate longitudinal household outcome data for a pathogen with
#'   two competing infectious compartments, for example tweo strains of an
#'   infection where infection with one confers immunity against the other.
#'   Members move through four states (1 = susceptible, 2 = first infectious
#'   compartment, 3 = second infectious compartment, 4 = recovered).
#'
#' @details The two infectious compartments can have separate intra-household
#'   infection probabilities but share an extra-household infection probability.
#'   As in [sim_sir()], each true state generates binary observations through
#'   `obs_prob`.
#'
#' @param eh_prob Baseline per step probability that a susceptible member is
#'   infected from outside the household.
#' @param ih_prob Baseline per step intra household infection probability. May
#'   be a single value (shared across both infectious compartments) or a length
#'   two vector (one per infectious compartment).
#' @param n_hh Number of households to simulate.
#' @param hh_size Vector of possible household sizes; each household size is
#'   drawn from this vector with replacement.
#' @param tmax Number of time steps to simulate per household.
#' @param gamma Length two vector of per step recovery probabilities, one for
#'   each infectious compartment.
#' @param covs_eh Numeric vector of coefficients for the extra household
#'   covariates, on the logit scale.
#' @param covs_ih Numeric vector of coefficients for the intra household
#'   covariates, on the logit scale.
#' @param obs_prob List of observation types. Each element is a vector giving
#'   the probability of a positive observation for that type in each true state.
#' @param start_prob Vector of starting state probabilities used to draw each
#'   member's state on the first time step.
#' @param complete_enroll If `TRUE`, return observations for every household
#'   member. If `FALSE`, return only enrolled members in `obs`.
#'
#' @return A list with `obs` (the observed rows) and `complete_obs` (every
#'   member at every time step). Each observation data frame has columns `t`,
#'   `part_id`, `enroll`, `state`, `hh_size`, `hh_id`, and one `y` column per
#'   observation type.
#'
#' @examples
#' # Simulate a SIIR study splitting infections 70/30 across two compartments
#' sim_siir_compete(
#'   n_hh = 20,
#'   tmax = 50,
#'   ih_prob = c(0.05, 0.025),
#'   eh_prob = 0.01,
#'   gamma = c(1 / 5, 1 / 3),
#' )
#'
#' @keywords internal
sim_siir_compete <- function(
  eh_prob = 0.01,
  ih_prob = 0.05,
  n_hh = 100,
  hh_size = 1:5,
  tmax = 100,
  gamma = c(1 / 5, 1 / 10),
  covs_eh = c(0, 0),
  covs_ih = c(0, 0),
  obs_prob = list(
    c(0.05, 0.95, 0.10, 0.05),
    c(0.05, 0.10, 0.95, 0.05),
    c(0.01, 0.01, 0.01, 0.8)
  ),
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
            eh_prob_x <- plogis(
              qlogis(eh_prob) +
                sum(x[last_x + part, ] * covs_eh)
            )
            ih_prob_x <- plogis(
              qlogis(ih_prob) +
                sum(x[last_x + part, ] * covs_ih)
            )
            no_inf_prob <- (1 - eh_prob_x) * (1 - ih_prob_x)^prior_inf

            p_leave <- 1 - prod(no_inf_prob)
            p_leave_partition <- p_leave *
              (1 - no_inf_prob) /
              sum(1 - no_inf_prob)

            new_states[part] <- sample(
              x = c(1, 2, 3, 4),
              size = 1,
              prob = c(
                1 - p_leave,
                p_leave_partition[1],
                p_leave_partition[2],
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
