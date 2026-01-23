
logit <- function(x) qlogis(x)

inv_logit <- function(x) plogis(x)

#' @title Checking split parameter specification
#' @description Utility function for checking whether split parameter is properly defined
#'
#' @param to character vector giving the name(s) of the destination compartment
#' @param split numeric or character vector indicating how people moving out of the starting compartment are split between the destination compartments
#'
split_check <- function(to, split) {
  # Make sure splits are valid
  if(!(length(split == 1) & sum(is.na(split)) == 1)) { # single NA OK
    if(sum(is.na(split)) >= 1) { # otherwise NA not OK
      stop("Improper specification for split - cannot contain NA")
    } else if(is.numeric(split)) { # numeric
      if(length(to) == 1) { # single destination compartment
        if(length(split != 1)) {
          stop("Must have length(split) == 1 if length(to) == 1.")
        }
      } else { # multiple destination compartments
        if(length(split) == length(to)) {
          if(sum(split) != 1) {
            stop("Values for split must sum to 1 if providing for all destination compartments.")
          }
        } else { # length(split) != length(to)
          if(length(split) == length(to) - 1) {
            if(sum(split) > 1) {
              stop("split cannot sum to greater than 1")
            }
          } else {
            stop("length(split) must be equal to length(to) or length(to)-1")
          }
        }
      }
    } else { # not numeric
      if(length(to) == 1) {
        if(length(split) != 1) {
          stop("Must have length(split) == 1 if length(to) == 1.")
        }
      } else { # multiple destination compartments
        if(length(split) != length(to)-1) {
          stop("Number of parameter names for split must be one less than the number of destination compartments.")
        }
      }
    }
  }
}


#' @title Define a non-infection transition
#'
#' @description
#' Defines a state transition in the infection process model which does not
#' represent an transmission (infection) event
#'
#' @param from string giving name of origin compartment
#' @param to string or vector of strings giving the name of the destination compartment
#' @param split an optional character or numeric vector indicating what proportion of individuals
#' transition into each of the `to` compartments
#'
transit <- function(from, to, split = NA, ...) {

  .dots <- unlist(list(...))

  split_check(to, split)

  out <- list()
  for(i in 1:length(to)) {
    out[[i]] <- data.frame(from = from,
                           to = to[i],
                           source = NA,
                           rate_name = NA,
                           rate_value = NA,
                           split_name = NA,
                           split_value = NA)

    out[[i]]$rate_name <- names(.dots)
    out[[i]]$rate_value <- .dots

    if(i == 1) {
      if(is.numeric(split[i])) {
        out[[i]]$split_value <- split[i]
      } else {
        out[[i]]$split_name <- split[i]
      }
    } else {
      if(is.numeric(split)) {
        if(length(split) >= i) {
          out[[i]]$split_value <- split[i]
        } else {
          out[[i]]$split_value <- 1-sum(split[1:(i-1)])
        }
      } else {
        if(length(split) >= i) {
          out[[i]]$split_name <- split[i]
        } else {
          out[[i]]$split_name <- paste0("1-", paste(split[1:(i-1)], sep = "-", collapse = "-"))
        }
      }
    }

  }

  return(bind_rows(out))

}

#' @title Define a infection transition
#'
#' @description
#' Defines a state transition in the infection process model which is
#' the result of a transmission (infection) event
#'
#' @param from string giving name of origin compartment
#' @param to string giving the name of the destination compartment
#' @param source string (or vector of strings) designating which compartments are infectious. If NULL, the destination compartment is presumed to be the infectious compartment.
#'
transmit <- function(from, to, source = NA, split = NA) {

  split_check(to, split)

  out <- list()

  for(i in 1:length(to)) {

    out[[i]] <- data.frame(from = from,
                           to = to[i],
                           rate_name = NA,
                           rate_value = NA,
                           split_name = NA,
                           split_value = NA,
                           mult_inf_prob = NA)

    out[[i]]$source <- ifelse(sum(is.na(source))>0, to, list(source))

    if(i == 1) {
      if(is.numeric(split[i])) {
        out[[i]]$split_value <- split[i]
      } else {
        out[[i]]$split_name <- split[i]
      }
    } else {
      if(is.numeric(split)) {
        if(length(split) >= i) {
          out[[i]]$split_value <- split[i]
        } else {
          out[[i]]$split_value <- 1-sum(split[1:(i-1)])
        }
      } else {
        if(length(split) >= i) {
          out[[i]]$split_name <- split[i]
        } else {
          out[[i]]$split_name <- paste0("1-", paste(split[1:(i-1)], sep = "-", collapse = "-"))
        }
      }
    }
  }

  return(bind_rows(out))

}


#' Builds infection process model.
#'
#' @param ... a series of transit or transmit function calls
#' @param ih_cov indicator for whether intra-household infection risk is a function of covariates
#' @param eh_cov indicator for whether extra-household infection risk is a function of covariates
#'
make_infection_model <- function(..., mult_inf_probs = FALSE) {
  .dots <- list(...)

  out  <- dplyr::bind_rows(.dots)
  out$mult_inf_probs <- mult_inf_probs

  return(out)

}

#' Creates transmission probability matrix,
#'
#' @param inf_model infection process model object yielded by make_infection_model()
#'
get_transmission_details <- function(inf_model) {
  states <- unique(c(inf_model$from, inf_model$to))
  trans <- matrix(1e-10, nrow = length(states), ncol = length(states))
  rownames(trans) <- paste("to", states, sep = "_")
  colnames(trans) <- paste("from", states, sep = "_")
  mult <- matrix(1, nrow = length(states), ncol = length(states))
  rownames(mult) <- paste("to", states, sep = "_")
  colnames(mult) <- paste("from", states, sep = "_")

  trans_to_fit <- data.frame(from = character(),
                             to = character(),
                             trans_row = numeric(),
                             trans_col = numeric(),
                             source = list(),
                             rate_name = character(),
                             param = numeric())

  mult_to_fit <- data.frame(from = character(),
                            to = character(),
                            mult_row = numeric(),
                            mult_col = numeric(),
                            mult_name = character(),
                            param = numeric())

  for(i in 1:nrow(inf_model)) {

    # Transition probabilities
    if(!is.na(inf_model$rate_value[i])) {
      trans[states == inf_model$to[i],states == inf_model$from[i]] <- inf_model$rate_value[i]
    } else {
      temp <- data.frame(from = inf_model$from[i],
                         to = inf_model$to[i],
                         trans_row = which(states == inf_model$to[i]),
                         trans_col = which(states == inf_model$from[i]),
                         rate_name = inf_model$rate_name[i],
                         param = NA)
      temp$source <- ifelse(is.null(inf_model$source[i][[1]]), list(0) , list(which(states %in% inf_model$source[i][[1]])))
      trans_to_fit <- bind_rows(trans_to_fit, temp)
    }

    # Multipliers
    if(!is.na(inf_model$split_value[i])) {
      mult[states == inf_model$to[i],states == inf_model$from[i]] <- inf_model$split_value[i]
    } else if(!is.na(inf_model$split_name[i])) {
      temp <- data.frame(from = inf_model$from[i],
                         to = inf_model$to[i],
                         mult_row = which(states == inf_model$to[i]),
                         mult_col = which(states == inf_model$from[i]),
                         mult_name = inf_model$split_name[i],
                         param = NA)
      mult_to_fit <- bind_rows(mult_to_fit, temp)
    }
  }

  # Identify unique parameters to fit - transitions
  if(sum(!is.na(trans_to_fit$rate_name)) > 0) {
    fac_levels <- unique(trans_to_fit$rate_name[!is.na(trans_to_fit$rate_name)])
    trans_to_fit$param <- as.numeric(factor(trans_to_fit$rate_name, levels = fac_levels))
  }
  trans_to_fit <- trans_to_fit %>%
    mutate(param = replace_na(param, 0))

  # Identify unique parameters to fit - multipliers
  if(sum(!is.na(mult_to_fit$mult_name)) > 0) {
    fac_levels <- unique(mult_to_fit$mult_name[!is.na(mult_to_fit$mult_name) & !str_detect(mult_to_fit$mult_name, "1-")])
    mult_to_fit$param <- as.numeric(factor(mult_to_fit$mult_name, levels = fac_levels))
  }
  mult_to_fit <- mult_to_fit %>%
    mutate(param = replace_na(param, 0))

  if(nrow(mult_to_fit) > 0) {
    for(i in 1:nrow(mult_to_fit)) {
      if(grepl("1-", mult_to_fit$mult_name[i])) {
        to_match <- strsplit(mult_to_fit$mult_name[i], "-")[[1]]
        to_match <- to_match[to_match != "1"]
        matches <- mult_to_fit$param[mult_to_fit$mult_name %in% to_match]
        mult_to_fit$param[i] <- list(-1*matches)
      }
    }
  }

  return(list(states = states,
              trans_matrix = trans,
              mult_matrix = mult,
              trans_to_fit = trans_to_fit,
              mult_to_fit = mult_to_fit,
              inf_states = unique(unlist(trans_to_fit$source[is.na(trans_to_fit$rate_name)])),
              mult_inf_probs = inf_model$mult_inf_probs[1]))
}


make_observation_model <- function(...) {
  .dots <- list(...)

  ops <- list()
  for(i in 1:length(.dots)) {
    op <- matrix(nrow = 2, ncol = length(.dots[[i]]))
    op[1,] <- 1-.dots[[i]]
    op[2,] <- .dots[[i]]
    rownames(op) <- c("neg_obs", "pos_obs")
    colnames(op) <- names(.dots[[i]])
    ops[[i]] <- op
  }
  names(ops) <- names(.dots)
  return(ops)
}

make_stan_data <- function(inf_model, obs_model, data, init_probs, epsilon = 1e-10,
                           ih_cov = NULL, eh_cov = NULL) {


  inf_details <- get_transmission_details(inf_model)
  dat <- data %>%
    arrange(hh_id, t, part_id)

  dat$row_id <- 1:nrow(dat)
  hh_sum <- dat %>%
    group_by(hh_id, hh_size) %>%
    summarize(hh_start_ind = min(row_id),
              hh_end_ind = max(row_id),
              hh_tmin = min(t),
              hh_tmax = max(t),
              obs_per_hh = n()) %>%
    ungroup()

  source_state_matrix <- matrix(0, nrow = nrow(inf_details$trans_to_fit), ncol = length(inf_details$states))
  for(i in 1:nrow(inf_details$trans_to_fit)) {
    if(any(inf_details$trans_to_fit$source[[i]] == 0)) {
      next
    } else {
      source_state_matrix[i,inf_details$trans_to_fit$source[[i]]] <- 1
    }
  }

  obs_array <- array(dim = c(length(obs_model),2, length(inf_details$states)))
  for(i in 1:length(obs_model)) {
    obs_array[i,,] <- obs_process[[i]]
  }

  # Expand multipliers if needed
  if(is.list(inf_details$mult_to_fit$param)) {
    mult_info <- data.frame()
    for(i in 1:nrow(inf_details$mult_to_fit)) {
      if(length(inf_details$mult_to_fit$param[i][[1]]) == 1) {
        mult_info <- bind_rows(mult_info, inf_details$mult_to_fit[i,])
      } else {
        for(j in 1:length(inf_details$mult_to_fit$param[i][[1]])) {
          temp <- inf_details$mult_to_fit[i,] %>%
            unnest(param)
          mult_info <- bind_rows(mult_info, temp)
        }
      }
    }
  } else {
    mult_info <- inf_details$mult_to_fit
  }

  # TODO: deal with missing observations (change NA to -1)

  dat_stan <- list(n_states = length(inf_details$states),
                   trans = inf_details$trans_matrix,
                   n_inf_states = length(inf_details$inf_states),
                   inf_states = array(inf_details$inf_states),
                   n_trans_fit = nrow(inf_details$trans_to_fit),
                   param_index = array(inf_details$trans_to_fit$param),
                   trans_index = inf_details$trans_to_fit %>% select(trans_row, trans_col),
                   source_states = source_state_matrix,
                   multiplier = inf_details$mult_matrix,
                   n_mult_fit = nrow(mult_info),
                   n_mult_params = length(unique(abs(unlist(mult_info$param)))),
                   mult_param_index = unlist(mult_info$param),
                   mult_index = mult_info %>% select(mult_row, mult_col),
                   n_params = length(unique(inf_details$trans_to_fit$param[inf_details$trans_to_fit$param != 0])),
                   n_hh = max(dat$hh_id),
                   hh_size = hh_sum$hh_size,
                   n_obs = nrow(dat),
                   n_obs_type = length(obs_process),
                   n_unique_obs = 2, #TODO: allow multi-level outcomes
                   y = dat %>% select(names(obs_model)) + 1,
                   part_id = dat$part_id,
                   t_day = dat$t,
                   obs_per_hh = hh_sum$obs_per_hh,
                   hh_start_ind = hh_sum$hh_start_ind,
                   hh_end_ind = hh_sum$hh_end_ind,
                   hh_tmin = hh_sum$hh_tmin,
                   hh_tmax = hh_sum$hh_tmax,
                   obs_prob = obs_array,
                   init_probs = init_probs, #TODO: Toggle to fit
                   epsilon = epsilon,
                   n_inf_prob = ifelse(inf_details$mult_inf_probs, length(inf_details$inf_states), 1))

  # Get covariate information if applicable
  if(!(is.null(eh_cov) & is.null(ih_cov))) {
    if(!is.null(eh_cov) & !is.null(ih_cov)) {

      k_ih <- length(ih_cov)
      x_ih <- data[,ih_cov]

      k_eh <- length(eh_cov)
      x_eh <- data[,eh_cov]

      dat_stan <- append(dat_stan,
                         list(k_ih = k_ih,
                              x_ih = x_ih,
                              k_eh = k_eh,
                              x_eh = x_eh))

    } else {
      stop("Currently only support having covariates on both intra- and extra-household infeciton probabilities.")
    }
  }

  return(dat_stan)

}

# TODO: option to save state probabilities
# TODO: create processed model results
run_model <- function(inf_model, obs_model, data, init_probs, epsilon = 1e-10,
                      ih_cov = NULL, eh_cov = NULL,
                      file = "stan/hmm.stan", iter = 2000, chains = 4,
                      cores = getOption("mc.cores", 1L), init = NULL,
                      save_chains = TRUE, save_states = TRUE) {

  dat_stan <- make_stan_data(inf_model, obs_model, data, init_probs, epsilon, ih_cov, eh_cov)

  if(is.null(init)) {
    init = rep(list(list(logit_params = array(rep(logit(0.5), dat_stan$n_params)),
                         logit_mult_params = array(rep(logit(0.5), dat_stan$n_mult_params)),
                         beta_eh = logit(0.02),
                         beta_ih = array(rep(logit(0.02), dat_stan$n_inf_prob)))), 4)
  } else {
    init <- rep(list(init), 4)
  }

  if(save_states) {
    stan_fit <- stan(file = file,
                     data = dat_stan,
                     iter = iter,
                     chains = chains,
                     cores = cores,
                     init = init)
  } else {
    stan_fit <- stan(file = file,
                     data = dat_stan,
                     iter = iter,
                     chains = chains,
                     cores = cores,
                     init = init,
                     pars = "logalpha",
                     include = FALSE)
  }


  ch <- rstan::extract(stan_fit)

#   res <- bind_rows(summarize_chains(ch, "ih_prob"),
#                    summarize_chains(ch, "eh_prob"))

  return(list(# res = res,
              chains = if(save_chains){ch} else {NULL}))

}

# TODO: handle multi-dimensional parameters
summarize_chains <- function(ch, param, quantiles = c(0.025, 0.975)) {
  ch_param <- ch[[param]]

  # Get mean an median
  out <- data.frame(param = param, mean = mean(ch_param), median = median(ch_param))

  # Calculate additional quantiles
  qs <- list()
  for(i in 1:length(quantiles)) {
    qs[[i]] <- data.frame(x = quantile(ch_param, quantiles[i]))
    names(qs[[i]]) <- paste0("quartile_", quantiles[i])
  }
  qs <- bind_cols(qs)
  out <- bind_cols(out, qs)
  rownames(out) <- NULL

  return(out)
}



