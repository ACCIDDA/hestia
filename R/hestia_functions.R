#' @title Logit Transformation
#' @param x Numeric value between 0 and 1
#' @returns A numeric value.
#'
#' @keywords internal
logit <- function(x) {
  qlogis(x)
}

#' @title Inverse Logit Transformation
#' @param x Numeric value
#' @returns A numeric value between 0 and 1.
#'
#' @keywords internal
inv_logit <- function(x) {
  plogis(x)
}

#' @title Checking split parameter specification
#' @description Utility function for checking whether split parameter is
#'   properly defined
#'
#' @param to Character vector giving the name(s) of the destination compartment
#' @param split Numeric or character vector indicating how people moving out of
#'   the starting compartment are split between the destination compartments
#'
#' @keywords internal
split_check <- function(to, split) {
  # Make sure splits are valid
  if (length(split == 1) && sum(is.na(split)) == 1) {
    if (length(to) > 1) {
      stop("Need to specify a split if length(to) > 1")
    }
  } else {
    # single NA OK
    if (sum(is.na(split)) >= 1) {
      # otherwise NA not OK
      stop("Improper specification for split - cannot contain NA")
    } else if (!(is.character(split) || is.numeric(split))) {
      # Must be character or numeric
      stop("Must provide character or numeric values for split")
    } else if (is.numeric(split)) {
      # numeric
      if (sum(split < 0) > 0) {
        stop("Splits cannot be negative")
      }
      if (length(to) == 1) {
        stop("Cannot specify a split if length(to) == 1")
      } else {
        # multiple destination compartments
        if (length(split) == length(to)) {
          if (sum(split) != 1) {
            stop(
              "Values for split must sum to 1 if providing for all ",
              "destination compartments."
            )
          }
        } else {
          if (length(split) == length(to) - 1) {
            if (sum(split) > 1) {
              stop("Splits cannot sum to greater than 1")
            }
          } else {
            stop("length(split) must be equal to length(to) or length(to)-1")
          }
        }
      }
    } else {
      # not numeric
      if (length(to) == 1) {
        stop("Must have length(split) == 1 if length(to) == 1.")
      } else {
        # multiple destination compartments
        if (length(split) != length(to) - 1) {
          stop(
            "Number of parameter names for split must be one less than ",
            "the number of destination compartments."
          )
        }
      }
    }
  }
}

#' @title Checking to and from parameter specification
#' @description Utility function for checking whether to and from parameters
#'   are properly defined
#'
#' @param to Character vector giving the name(s) of the destination compartment
#' @param from string giving name of origin compartment
#'
#' @keywords internal
to_from_check <- function(to, from) {
  if (length(from) == 0 || length(to) == 0) {
    stop("'to' and 'from' arguments cannot be empty")
  }

  if (length(from) != 1) {
    stop("'from' must be of length one")
  }

  if (!is.character(from) || !is.character(to)) {
    stop("'to' and 'from' arguments must be characters")
  }

  if (length(intersect(to, from)) > 0) {
    stop("'to' and 'from' must be mutually exclusive")
  }
}


#' @title Define a non-infection transition
#'
#' @description Defines a state transition in the infection process model which
#'   does not represent an transmission (infection) event
#'
#' @param from string giving name of origin compartment
#' @param to string or vector of strings giving the name of the destination
#'   compartment
#' @param split an optional character or numeric vector indicating what
#'   proportion of individuals transition into each of the `to` compartments
#' @param ... Argument specifying the name and optionally the value of the
#'   transition rate parameter. For example if the parameter name is "gamma" and
#'   it is being fit this argument should be `gamma = NA`. If we want to provide
#'   a numeric value, x, for "gamma" this argument should be `gamma = x`.
#' @returns Data frame with a row for each unique transition and the following
#'   columns:
#'  - "from", the name of the origin compartment
#'  - "to", the name of the destination compartment
#'  - "source", NA (so output columns match with `transmit`)
#'  - "rate_name", name of transition rate parameter if fitting
#'  - "rate_value", numeric value of transition rate if not fitting
#'  - "split_name", name of split parameter if fitting
#'  - "split_value", numeric value of split if not fitting
#'
#' @examples
#' # Transition from compartment A to B which occurs at rate gamma where
#' # gamma is a user-specified value of 0.2
#' progress(from = "A", to = "B", gamma = 0.2)
#'
#' @examples
#' # Transition from compartment A to B which occurs at rate gamma where
#' # gamma is fit
#' progress(from = "A", to = "B", gamma = NA)
#'
#' @examples
#' # Transition from compartment A split into two destination compartments, B1
#' # and B2, were individuals leave A at rate delta (fit) and phi (fit) is the
#' # proportion of those who leave A who go into B1 (i.e. 1-phi go into B2).
#' # gamma is fit
#' progress(from = "A", to = c("B1", "B2"), split = "phi", delta = NA)
#'
#'
#' @export
progress <- function(from, to, split = NA, ...) {
  .dots <- unlist(list(...))

  # CHecks on parameter specification
  to_from_check(to, from)
  split_check(to, split)

  out <- list()
  for (i in seq_along(to)) {
    out[[i]] <- data.frame(
      from = from,
      to = to[i],
      rate_name = NA,
      rate_value = NA,
      split_name = NA,
      split_value = NA
    )

    out[[i]]$rate_name <- names(.dots)
    out[[i]]$rate_value <- .dots

    if (i == 1) {
      if (is.numeric(split[i])) {
        out[[i]]$split_value <- split[i]
      } else {
        out[[i]]$split_name <- split[i]
      }
    } else {
      if (is.numeric(split)) {
        if (length(split) >= i) {
          out[[i]]$split_value <- split[i]
        } else {
          out[[i]]$split_value <- 1 - sum(split[1:(i - 1)])
        }
      } else {
        if (length(split) >= i) {
          out[[i]]$split_name <- split[i]
        } else {
          out[[i]]$split_name <- paste0(
            "1-",
            paste(split[1:(i - 1)], sep = "-", collapse = "-")
          )
        }
      }
    }
  }

  out <- dplyr::bind_rows(out)

  if (!("rate_name" %in% colnames(out))) {
    stop("No paramter name provided for transition rate.")
  }

  out
}

#' @title Define a infection transition
#'
#' @description Defines a state transition in the infection process model which
#' is the result of a transmission (infection) event
#'
#' @param from string giving name of origin compartment
#' @param to string giving the name of the destination compartment
#' @param source string (or vector of strings) designating which compartments
#'   are infectious. If NA, the destination compartment is presumed to be the
#'   infectious compartment.
#' @param split an optional character or numeric vector indicating what
#'   proportion of individuals transition into each of the `to` compartments
#' @returns Data frame with a row for each unique transition and the following
#' columns:
#'  - "from", the name of the origin compartment
#'  - "to", the name of the destination compartment
#'  - "split_name", name of split parameter if fitting
#'  - "split_value", numeric value of split if not fitting
#'  - "source", names of compartment(s) that are the source of new infections
#'  i.e. the infectious compartments
#'
#' @examples
#' # Specify S > E transition in an SEIR model
#' transmit(from = "S", to = "E", source = "I")
#'
#' @examples
#' # Upon infection, individuals either enter Is (symptomatic infections) or
#' # Ia (asymptomatic infection) where 30% of infections are symptomatic.
#' transmit(from = "S", to = c("Is", "Ia"), split = 0.3)
#'
#' @export
transmit <- function(from, to, source = NA, split = NA) {
  # Checks on parameter specification
  to_from_check(to, from)
  split_check(to, split)

  out <- list()

  for (i in seq_along(to)) {
    out[[i]] <- data.frame(
      from = from,
      to = to[i],
      split_name = NA,
      split_value = NA
    )

    out[[i]]$source <- ifelse(sum(is.na(source)) > 0, list(to), list(source))

    if (i == 1) {
      if (is.numeric(split[i])) {
        out[[i]]$split_value <- split[i]
      } else {
        out[[i]]$split_name <- split[i]
      }
    } else {
      if (is.numeric(split)) {
        if (length(split) >= i) {
          out[[i]]$split_value <- split[i]
        } else {
          out[[i]]$split_value <- 1 - sum(split[1:(i - 1)])
        }
      } else {
        if (length(split) >= i) {
          out[[i]]$split_name <- split[i]
        } else {
          out[[i]]$split_name <- paste0(
            "1-",
            paste(split[1:(i - 1)], sep = "-", collapse = "-")
          )
        }
      }
    }
  }

  dplyr::bind_rows(out)
}

#' @title Build Infection Process Model
#'
#' @param ... a series of \link{progress} or \link{transmit} function calls
#' @param mult_inf_probs If FALSE then all infection probabilities are shared
#' across infectious compartments
#' @returns A data frame with a row for each unique transition in the infection
#' process model and the following columns:
#'  - "from", the name of the origin compartment
#'  - "to", the name of the destination compartment
#'  - "source", NA (so output columns match with `transmit`)
#'  - "rate_name", name of transition rate parameter if fitting
#'  - "rate_value", numeric value of transition rate if not fitting
#'  - "split_name", name of split parameter if fitting
#'  - "split_value", numeric value of split if not fitting
#'  - "source", names of compartment(s) that are the source of new infections
#'  - "mult_inf_probs", a logical indicator for whether infectious compartments
#'  have separate (`true`) or shared (`FALSE`) intra-household infection
#'  probabilities
#'
#' @examples
#' # Basic SIR model with recovery rate gamma to be fit
#' make_infection_model(transmit(from = "S", to = "I"),
#'                      progress(from = "I", to = "R", gamma = NA))
#'
#' @examples
#' # Split the infectious compartment into symptomatic (I_s) and asymptomatic
#' # (I_a) infections, with separate intra-household infection probabilities.
#' make_infection_model(
#'   transmit(from = "S",
#'            to = c("I_s", "I_a"),
#'            source = c("I_s", "I_a"),
#'            split = "phi"),
#'   progress(from = "I_s", to = "R", gamma_s = NA),
#'   progress(from = "I_a", to = "R", gamma_a = NA),
#'   mult_inf_probs = TRUE)
#'
#' @export
make_infection_model <- function(..., mult_inf_probs = FALSE) {
  .dots <- list(...)

  out <- dplyr::bind_rows(.dots)
  out$mult_inf_probs <- mult_inf_probs

  # Check that there is at least one transmit input
  if (!("source" %in% names(out))) {
    stop(
      "No infection process specified. Must have at least one transmit input."
    )
  }

  # Check for discontinuity
  states <- unique(c(out$from, out$to))
  graph <- matrix(0, nrow = length(states), ncol = length(states))
  rownames(graph) <- states
  colnames(graph) <- states
  for (i in seq_len(nrow(out))) {
    graph[out$from[i], out$to[i]] <- 1
    graph[out$to[i], out$from[i]] <- 1
  }

  if (!graph_connected(graph)) {
    stop("Compartmental model structure cannot be disconnected.")
  }

  # Check that source compartments are defined elsewhere in model
  sources <- unlist(out$source[!is.na(out$source)])
  if (length(which(!(sources %in% states))) > 0) {
    stop("Compartments specified in source argument do not exist in model.")
  }

  out
}

#' Check whether an undirected adjacency matrix is fully connected
#'
#' Breadth-first reachability from the first node. Replaces a single-use
#' dependency on Claddis::is_graph_connected (the only reason hestia imported
#' Claddis, whose heavy dependency tree failed to build on R-devel macOS).
#'
#' @param graph square symmetric 0/1 adjacency matrix.
#' @return TRUE if every node is reachable from the first, else FALSE.
#' @keywords internal
graph_connected <- function(graph) {
  n <- nrow(graph)
  if (n <= 1) {
    return(TRUE)
  }
  visited <- rep(FALSE, n)
  visited[1] <- TRUE
  queue <- 1L
  while (length(queue) > 0) {
    node <- queue[1]
    queue <- queue[-1]
    nbrs <- which(graph[node, ] != 0 & !visited)
    visited[nbrs] <- TRUE
    queue <- c(queue, nbrs)
  }
  all(visited)
}

#' @title Create Transmission Probability Matrix
#'
#' @param inf_model infection process model object yielded by
#'   make_infection_model()
#' @returns A list with the following entries:
#'  - "states", a character vector of state names
#'  - "trans_matrix", a matrix with transition rate values between origin
#'  (column) and destination (row) compartments where given. Where no values
#'  are provided for transition rates the entry is equal to 1e-10.
#'  - "mult_matrix", a matrix with split values for the transition between
#'  origin (column) and destination (row) compartments where given. Where no
#'  split values are provided for transition rates the entry is equal to 1.
#'  - "trans_to_fit", a data frame with information on all transitions within
#'  the model, with parameter indexing details if the transition rate is being
#'  fit
#'  - "mult_to_fit", a data frame with information on all splits specified in
#'  the model, with parameter indexing details if the split is being fit
#'  - "inf_states", numeric values corresponding to the indices in "states"
#'  which are the infectious states
#'  - "mult_inf_probs", indicator for whether infectious states have shared or
#'  separate inta-household infection probabilities. If FALSE then all infection
#'  probabilities are shared across infectious compartments.
#'
#' @examples
#' # Split the infectious compartment into symptomatic (I_s) and asymptomatic
#' # (I_a) infections, with separate intra-household infection probabilities.
#' inf_model <- make_infection_model(
#'   transmit(from = "S",
#'            to = c("I_s", "I_a"),
#'            source = c("I_s", "I_a"),
#'            split = "phi"),
#'   progress(from = "I_s", to = "R", gamma_s = NA),
#'   progress(from = "I_a", to = "R", gamma_a = NA),
#'   mult_inf_probs = TRUE)
#'
#' get_transmission_details(inf_model)
#'
#'
#' @global param
#'
#' @export
get_transmission_details <- function(inf_model) {
  states <- unique(c(inf_model$from, inf_model$to))
  trans <- matrix(1e-10, nrow = length(states), ncol = length(states))
  rownames(trans) <- paste("to", states, sep = "_")
  colnames(trans) <- paste("from", states, sep = "_")
  mult <- matrix(1, nrow = length(states), ncol = length(states))
  rownames(mult) <- paste("to", states, sep = "_")
  colnames(mult) <- paste("from", states, sep = "_")

  trans_to_fit <- data.frame(
    from = character(),
    to = character(),
    trans_row = numeric(),
    trans_col = numeric(),
    source = list(),
    rate_name = character(),
    param = numeric()
  )

  mult_to_fit <- data.frame(
    from = character(),
    to = character(),
    mult_row = numeric(),
    mult_col = numeric(),
    mult_name = character(),
    param = numeric()
  )

  if (!("rate_value" %in% colnames(inf_model))) {
    inf_model$rate_value <- NA
  }

  for (i in seq_len(nrow(inf_model))) {
    # Transition probabilities
    if (!is.na(inf_model$rate_value[i])) {
      trans[
        states == inf_model$to[i],
        states == inf_model$from[i]
      ] <- inf_model$rate_value[i]
    } else {
      temp <- data.frame(
        from = inf_model$from[i],
        to = inf_model$to[i],
        trans_row = which(states == inf_model$to[i]),
        trans_col = which(states == inf_model$from[i]),
        rate_name = inf_model$rate_name[i],
        param = NA
      )
      temp$source <- ifelse(
        is.null(inf_model$source[i][[1]]),
        list(0),
        list(which(states %in% inf_model$source[i][[1]]))
      )
      trans_to_fit <- dplyr::bind_rows(trans_to_fit, temp)
    }

    # Multipliers
    if (!is.na(inf_model$split_value[i])) {
      mult[
        states == inf_model$to[i],
        states == inf_model$from[i]
      ] <- inf_model$split_value[i]
    } else if (!is.na(inf_model$split_name[i])) {
      temp <- data.frame(
        from = inf_model$from[i],
        to = inf_model$to[i],
        mult_row = which(states == inf_model$to[i]),
        mult_col = which(states == inf_model$from[i]),
        mult_name = inf_model$split_name[i],
        param = NA
      )
      mult_to_fit <- dplyr::bind_rows(mult_to_fit, temp)
    }
  }

  # Identify unique parameters to fit - transitions
  if (sum(!is.na(trans_to_fit$rate_name)) > 0) {
    fac_levels <- unique(trans_to_fit$rate_name[!is.na(trans_to_fit$rate_name)])
    trans_to_fit$param <- as.numeric(factor(
      trans_to_fit$rate_name,
      levels = fac_levels
    ))
  }
  trans_to_fit <- trans_to_fit |>
    dplyr::mutate(param = tidyr::replace_na(param, 0))

  # Identify unique parameters to fit - multipliers
  if (sum(!is.na(mult_to_fit$mult_name)) > 0) {
    fac_levels <- unique(
      mult_to_fit$mult_name[
        !is.na(mult_to_fit$mult_name) &
          !grepl("1-", mult_to_fit$mult_name, fixed = TRUE)
      ]
    )
    mult_to_fit$param <- as.numeric(factor(
      mult_to_fit$mult_name,
      levels = fac_levels
    ))
  }
  mult_to_fit <- mult_to_fit |>
    dplyr::mutate(param = tidyr::replace_na(param, 0))

  if (nrow(mult_to_fit) > 0) {
    for (i in seq_len(nrow(mult_to_fit))) {
      if (grepl("1-", mult_to_fit$mult_name[i])) {
        to_match <- strsplit(mult_to_fit$mult_name[i], "-")[[1]]
        to_match <- to_match[to_match != "1"]
        matches <- mult_to_fit$param[mult_to_fit$mult_name %in% to_match]
        mult_to_fit$param[i] <- list(-1 * matches)
      }
    }
  }

  list(
    states = states,
    trans_matrix = trans,
    mult_matrix = mult,
    trans_to_fit = trans_to_fit,
    mult_to_fit = mult_to_fit,
    inf_states = unique(
      unlist(trans_to_fit$source[is.na(trans_to_fit$rate_name)])
    ),
    mult_inf_probs = inf_model$mult_inf_probs[1]
  )
}

#' @title Make Observation Process Model
#'
#' @param ... A series of named vectors. Each vector corresponds to an
#'   observation type. Each entry in the vector is named for the corresponding
#'   compartment in the infection process model. The value is the probability of
#'   observing a positive observation given the individual is in the
#'   compartment.
#' @returns A list with one entry for each observation type.
#'
#' @examples
#' # Observation process for an SIR model with two observations, PCR test
#' # results and presence of IgG antibodies.
#' make_observation_model(
#'   pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
#'   igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
#' )
#'
#' @export
make_observation_model <- function(...) {
  .dots <- list(...)

  ops <- list()
  for (i in seq_along(.dots)) {
    op <- matrix(nrow = 2, ncol = length(.dots[[i]]))
    op[1, ] <- 1 - .dots[[i]]
    op[2, ] <- .dots[[i]]
    rownames(op) <- c("neg_obs", "pos_obs")
    colnames(op) <- names(.dots[[i]])
    ops[[i]] <- op
  }
  names(ops) <- names(.dots)

  if (length(unique(names(.dots))) < length(names(.dots))) {
    stop("Observation names must be unique")
  }

  ops
}


#' @title Create Data for Stan Model
#'
#' @param inf_model Infection process model object generated by
#'   \link{make_infection_model}
#' @param obs_model Observation process model object generated by
#'   \link{make_observation_model}
#' @param data data frame with columns the following columns:
#'  - "t", time of observation
#'  - "hh_id", numeric household ID which must be 1, 2, ..., N where N is the
#'  total number of household
#'  - "part_id", participant ID within the household which must be 1, 2, ... n
#'  where n is the number of observed household members
#'  - "hh_size", number of household members (whether observed or not)
#'  - At least one column with binary outcome information where the column
#'  name matches the outcome name provided to \link{make_observation_model}
#' @param init_probs vector of initial probabilities for each infection process
#'   model state
#' @param epsilon very small number to use for zero probability transitions
#' @param ih_cov NULL for run without covariates, otherwise data frame with
#'   intra-household covariates for each participant
#' @param eh_cov NULL for run without covariates, otherwise data frame with
#'   extra-household covariates for each participant
#' @returns A list with data elements for input to stan model.
#'
#' @global hh_id
#' @global part_id
#' @global hh_size
#' @global row_id
#' @global trans_row
#' @global trans_col
#' @global mult_row
#' @global mult_col
#'
#' @keywords internal
make_stan_data <- function(
  inf_model,
  obs_model,
  data,
  init_probs,
  epsilon = 1e-10,
  ih_cov = NULL,
  eh_cov = NULL
) {
  inf_details <- get_transmission_details(inf_model)
  dat <- data |>
    dplyr::arrange(hh_id, t, part_id)

  dat$row_id <- seq_len(nrow(dat))
  hh_sum <- dat |>
    dplyr::group_by(hh_id, hh_size) |>
    dplyr::summarize(
      hh_start_ind = min(row_id),
      hh_end_ind = max(row_id),
      hh_tmin = min(t),
      hh_tmax = max(t),
      obs_per_hh = dplyr::n()
    ) |>
    dplyr::ungroup()

  source_state_matrix <- matrix(
    0,
    nrow = nrow(inf_details$trans_to_fit),
    ncol = length(inf_details$states)
  )
  for (i in seq_len(nrow(inf_details$trans_to_fit))) {
    if (any(inf_details$trans_to_fit$source[[i]] == 0)) {
      next
    } else {
      source_state_matrix[i, inf_details$trans_to_fit$source[[i]]] <- 1
    }
  }

  obs_array <- array(dim = c(length(obs_model), 2, length(inf_details$states)))
  for (i in seq_along(obs_model)) {
    obs_array[i, , ] <- obs_model[[i]]
  }

  # Expand multipliers if needed
  if (is.list(inf_details$mult_to_fit$param)) {
    mult_info <- data.frame()
    for (i in seq_len(nrow(inf_details$mult_to_fit))) {
      if (length(inf_details$mult_to_fit$param[i][[1]]) == 1) {
        mult_info <- dplyr::bind_rows(mult_info, inf_details$mult_to_fit[i, ])
      } else {
        for (j in seq_along(inf_details$mult_to_fit$param[i][[1]])) {
          temp <- inf_details$mult_to_fit[i, ] |>
            tidyr::unnest(param)
          mult_info <- dplyr::bind_rows(mult_info, temp)
        }
      }
    }
  } else {
    mult_info <- inf_details$mult_to_fit
  }

  # TODO: deal with missing observations (change NA to -1)

  dat_stan <- list(
    n_states = length(inf_details$states),
    trans = inf_details$trans_matrix,
    n_inf_states = length(inf_details$inf_states),
    inf_states = array(inf_details$inf_states),
    n_trans_fit = nrow(inf_details$trans_to_fit),
    param_index = array(inf_details$trans_to_fit$param),
    trans_index = inf_details$trans_to_fit |>
      dplyr::select(trans_row, trans_col),
    source_states = source_state_matrix,
    transition_multiplier = inf_details$mult_matrix,
    n_mult_fit = nrow(mult_info),
    n_mult_params = length(unique(abs(unlist(mult_info$param)))),
    mult_param_index = unlist(mult_info$param),
    mult_index = mult_info |> dplyr::select(mult_row, mult_col),
    n_params = length(unique(
      inf_details$trans_to_fit$param[inf_details$trans_to_fit$param != 0]
    )),
    n_hh = max(dat$hh_id),
    hh_size = hh_sum$hh_size,
    n_obs = nrow(dat),
    n_obs_type = length(obs_model),
    n_unique_obs = 2,
    y = dat |> dplyr::select(names(obs_model)) + 1,
    part_id = dat$part_id,
    t_day = dat$t,
    obs_per_hh = hh_sum$obs_per_hh,
    hh_start_ind = hh_sum$hh_start_ind,
    hh_end_ind = hh_sum$hh_end_ind,
    hh_tmin = hh_sum$hh_tmin,
    hh_tmax = hh_sum$hh_tmax,
    obs_prob = obs_array,
    init_probs = init_probs,
    #TODO: Toggle to fit
    epsilon = epsilon,
    n_inf_prob = ifelse(
      inf_details$mult_inf_probs,
      length(inf_details$inf_states),
      1
    )
  )

  # Get covariate information if applicable
  if (!(is.null(eh_cov) && is.null(ih_cov))) {
    if (!is.null(eh_cov) && !is.null(ih_cov)) {
      k_ih <- ncol(ih_cov)
      x_ih <- ih_cov

      k_eh <- ncol(eh_cov)
      x_eh <- eh_cov

      dat_stan <- append(
        dat_stan,
        list(
          k_ih = k_ih,
          x_ih = x_ih,
          k_eh = k_eh,
          x_eh = x_eh
        )
      )
    } else {
      stop(
        "Currently only support having covariates on both intra- and ",
        "extra-household infeciton probabilities."
      )
    }
  }

  dat_stan
}

#' @title Run Household Transmission Model
#'
#' @param inf_model Infection process model object generated by
#'   \link{make_infection_model}
#' @param obs_model Observation process model object generated by
#'   \link{make_observation_model}
#' @param data data frame with columns the following columns:
#'  - "t", time of observation
#'  - "hh_id", numeric household ID which must be 1, 2, ..., N where N is the
#'  total number of household
#'  - "part_id", participant ID within the household which must be 1, 2, ... n
#'  where n is the number of observed household members
#'  - "hh_size", number of household members (whether observed or not)
#'  - At least one column with binary outcome information where the column
#'  name matches the outcome name provided to \link{make_observation_model}
#' @param init_probs vector of initial probabilities for each infection process
#'   model state
#' @param epsilon very small number to use for zero probability transitions
#' @param ih_cov NULL for run without covariates, otherwise data frame with
#'   intra-household covariates for each participant
#' @param eh_cov NULL for run without covariates, otherwise data frame with
#'   extra-household covariates for each participant
#' @param file file path to stan model
#' @param iter number of MCMC iterations
#' @param chains number of MCMC chains
#' @param cores number of cores for parallelization
#' @param init initial conditions for MCMC chains
#' @param save_states indicator for whether to save state probabilities
#' @returns `draws_array` object with chains for each model parameter
#'
#' @examples
#' # Subset sir package data to 10 households (for speedy example)
#' sir_sub <- sir[sir$hh_id <= 10, ]
#'
#' # Make infetion process model for SIR model
#' inf_mod <- make_infection_model(
#'   transmit(from = "S", to = "I"),
#'   progress(from = "I", to = "R", gamma = NA))
#'
#' # Make observation process model
#' obs_mod <- make_observation_model(
#'   pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
#'   igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8))
#'
#' run_model(inf_model = inf_mod,
#'           obs_model = obs_mod,
#'           data = sir_sub,
#'           init_probs = c(1 - 2 * 1e-10, 1e-10, 1e-10))
#'
#' @export
run_model <- function(
  inf_model,
  obs_model,
  data,
  init_probs,
  epsilon = 1e-10,
  ih_cov = NULL,
  eh_cov = NULL,
  iter = 2000,
  chains = 4,
  cores = getOption("mc.cores", 1L),
  init = NULL,
  save_states = FALSE
) {
  dat_stan <- make_stan_data(
    inf_model,
    obs_model,
    data,
    init_probs,
    epsilon,
    ih_cov,
    eh_cov
  )

  is_cov <- !is.null(eh_cov) && !is.null(ih_cov)

  if (is.null(init)) {
    if (is_cov) {
      init = rep(
        list(
          list(
            logit_params = array(rep(logit(0.5), dat_stan$n_params)),
            logit_mult_params = array(rep(
              logit(0.5),
              dat_stan$n_mult_params
            )),
            beta_eh = rep(0, dat_stan$k_eh),
            beta_ih = rep(0, dat_stan$k_ih),
            beta0_eh = logit(0.02),
            beta0_ih = array(rep(logit(0.02), dat_stan$n_inf_prob))
          )
        ),
        4
      )
    } else {
      init = rep(
        list(
          list(
            logit_params = array(rep(logit(0.5), dat_stan$n_params)),
            logit_mult_params = array(rep(
              logit(0.5),
              dat_stan$n_mult_params
            )),
            beta_eh = logit(0.02),
            beta_ih = array(rep(logit(0.02), dat_stan$n_inf_prob))
          )
        ),
        4
      )
    }
  } else {
    init <- rep(list(init), 4)
  }

  if (save_states) {
    if (is_cov) {
      stan_fit <- rstan::sampling(
        stanmodels$hmm_cov,
        data = dat_stan,
        iter = iter,
        chains = chains,
        cores = cores,
        init = init
      )
    } else {
      stan_fit <- rstan::sampling(
        stanmodels$hmm,
        data = dat_stan,
        iter = iter,
        chains = chains,
        cores = cores,
        init = init
      )
    }
  } else {
    if (is_cov) {
      stan_fit <- rstan::sampling(
        stanmodels$hmm_cov,
        data = dat_stan,
        iter = iter,
        chains = chains,
        cores = cores,
        init = init,
        pars = "logalpha",
        include = FALSE
      )
    } else {
      stan_fit <- rstan::sampling(
        stanmodels$hmm,
        data = dat_stan,
        iter = iter,
        chains = chains,
        cores = cores,
        init = init,
        pars = "logalpha",
        include = FALSE
      )
    }
  }

  stan_out <- list(
    stan_fit = stan_fit,
    stan_data = dat_stan,
    eh_cov_names = colnames(eh_cov),
    ih_cov_names = colnames(ih_cov)
  )

  rename_chains(inf_model, stan_out)
}

#' @title Extract chains and renames with user-provided parameter names
#'
#' @param inf_model Infection process model object generated by
#'   \link{make_infection_model}
#' @param model_output A list generated within \link{run_model} which contains
#'   the stanfit, the stan input data, and covariate names
#' @returns A `draws_array` object with chains for each model parameter
#'
#' @keywords internal
rename_chains <- function(inf_model, model_output) {
  # Get parameter information
  inf_details <- get_transmission_details(inf_model)
  mult_names <- unique(inf_details$mult_to_fit$mult_name[
    inf_details$mult_to_fit$param >= 1
  ])
  trans_names <- unique(inf_details$trans_to_fit$rate_name[
    inf_details$trans_to_fit$param >= 1
  ])
  n_ih <- ifelse(
    inf_details$mult_inf_probs == FALSE,
    1,
    length(inf_details$inf_states)
  )

  # Extract chains
  draws_full <- posterior::as_draws_array(model_output$stan_fit)

  # Get variable names and subset to parameters of interest
  var_names <- posterior::variables(draws_full)

  if (length(grep("beta0", var_names)) > 0) {
    var_names_sub <- c(
      grep("beta0_eh", var_names, value = TRUE),
      grep("beta0_ih", var_names, value = TRUE)
    )
  } else {
    var_names_sub <- c(
      grep("beta_eh", var_names, value = TRUE),
      grep("beta_ih", var_names, value = TRUE)
    )
  }

  if (length(trans_names) > 0) {
    var_names_sub <- c(
      var_names_sub,
      grep("logit_params", var_names, value = TRUE)
    )
  }

  if (length(mult_names) > 0) {
    var_names_sub <- c(
      var_names_sub,
      grep("logit_mult_params", var_names, value = TRUE)
    )
  }

  # Pull covariate coefficients if covariate run
  if (length(grep("beta0", var_names)) > 0) {
    var_names_sub <- c(
      var_names_sub,
      grep("beta_eh", var_names, value = TRUE),
      grep("beta_ih", var_names, value = TRUE)
    )
  }

  # Make new variable names from user-specified parameter names
  if (n_ih == 1) {
    ih_names <- "ih_prob"
  } else {
    ih_names <- paste0("ih_prob_", inf_details$states[inf_details$inf_states])
  }

  var_names_new <- c("eh_prob", ih_names, trans_names, mult_names)

  if (length(grep("beta0", var_names)) > 0) {
    var_names_new <- c(
      var_names_new,
      paste0(model_output$eh_cov_names, "_eh"),
      paste0(model_output$ih_cov_names, "_ih")
    )
  }

  # Subset to variables of interest and rename
  draws <- posterior::subset_draws(draws_full, variable = var_names_sub)
  posterior::variables(draws) <- var_names_new

  draws
}
