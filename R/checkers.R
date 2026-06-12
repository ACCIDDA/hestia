# Input-validation helpers for the exported hestia functions.
#
# These run at the top of each exported function so that bad inputs fail fast
# with a clear message, rather than surfacing later as a cryptic Stan, dplyr,
# or indexing error. They lean on `checkmate` for the common type/shape checks
# and add a few hestia-specific value checks on top.

#' Validate the `from`/`to` arguments shared by `transmit()` and `progress()`
#'
#' @param from origin compartment, expected to be a single non-empty string.
#' @param to destination compartment(s), expected to be a non-empty character
#'   vector.
#' @importFrom checkmate assert_character
#' @keywords internal
check_to_from <- function(from, to) {
  checkmate::assert_character(
    from,
    min.len = 1,
    any.missing = FALSE,
    .var.name = "from"
  )
  checkmate::assert_character(
    to,
    min.len = 1,
    any.missing = FALSE,
    .var.name = "to"
  )
  if (any(from %in% to)) {
    stop("'from' and 'to' must be different compartments.")
  }
  invisible(NULL)
}

#' Validate the `split` argument shared by `transmit()` and `progress()`
#'
#' `split` may be left as a single `NA` (no split), a character vector of
#' parameter names, or a numeric vector of proportions. This only checks the
#' broad type and finiteness; the detailed length/sum rules stay in
#' `split_check()`.
#'
#' @param split the split specification to validate.
#' @importFrom checkmate assert_numeric assert_character
#' @keywords internal
check_split <- function(split) {
  if (length(split) == 1 && is.na(split)) {
    return(invisible(NULL))
  }
  if (!(is.character(split) || is.numeric(split))) {
    stop("'split' must be NA, a character vector, or a numeric vector.")
  }
  if (is.numeric(split)) {
    checkmate::assert_numeric(
      split,
      any.missing = FALSE,
      finite = TRUE,
      .var.name = "split"
    )
  } else {
    checkmate::assert_character(
      split,
      min.len = 1,
      min.chars = 1,
      any.missing = FALSE,
      .var.name = "split"
    )
  }
  invisible(NULL)
}

#' Validate the `gamma`-style transition-rate values passed to `progress()`
#'
#' Each rate value must be either `NA` (fit the parameter) or a positive
#' numeric (use the supplied value).
#'
#' @param dots the unlisted `...` rate arguments from `progress()`.
#' @keywords internal
check_rate_dots <- function(dots) {
  if (length(dots) == 0) {
    stop("A transition-rate parameter must be supplied (e.g. gamma = NA).")
  }
  vals <- unlist(dots)
  for (v in vals) {
    if (is.na(v)) {
      next
    }
    if (!is.numeric(v) || v <= 0) {
      stop("Transition-rate values must be NA or a positive number.")
    }
  }
  invisible(NULL)
}

#' Validate a single observation-process specification
#'
#' Each observation is a named numeric vector of detection probabilities, one
#' entry per compartment, with values in `[0, 1]`.
#'
#' @param obs a single observation specification.
#' @param nm its name, used in error messages.
#' @importFrom checkmate assert_numeric
#' @keywords internal
check_observation <- function(obs, nm = "observation") {
  checkmate::assert_numeric(
    obs,
    min.len = 1,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    names = "named",
    .var.name = nm
  )
  invisible(NULL)
}

#' Validate the assembled model lists passed to `run_model()`
#'
#' @param inf_model output of `make_infection_model()`.
#' @param obs_model output of `make_observation_model()`.
#' @importFrom checkmate assert_data_frame assert_list
#' @keywords internal
check_models <- function(inf_model, obs_model) {
  checkmate::assert_data_frame(
    inf_model,
    min.rows = 1,
    .var.name = "inf_model"
  )
  checkmate::assert_list(
    obs_model,
    min.len = 1,
    names = "named",
    .var.name = "obs_model"
  )
  invisible(NULL)
}

#' Validate the observation `data` frame passed to `run_model()`
#'
#' @param data the participant-level data frame.
#' @param obs_model output of `make_observation_model()`; its names must each
#'   appear as an outcome column in `data`.
#' @importFrom checkmate assert_data_frame
#' @keywords internal
check_run_data <- function(data, obs_model) {
  checkmate::assert_data_frame(data, min.rows = 1, .var.name = "data")
  required <- c("t", "hh_id", "part_id", "hh_size")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "'data' is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "."
    )
  }
  missing_obs <- setdiff(names(obs_model), names(data))
  if (length(missing_obs) > 0) {
    stop(
      "'data' is missing outcome column(s) named in the observation model: ",
      paste(missing_obs, collapse = ", "),
      "."
    )
  }
  invisible(NULL)
}

#' Validate the `init_probs` argument passed to `run_model()`
#'
#' Initial state probabilities must be a finite non-negative numeric vector
#' that sums to roughly one.
#'
#' @param init_probs the vector of initial state probabilities.
#' @importFrom checkmate assert_numeric
#' @keywords internal
check_init_probs <- function(init_probs) {
  checkmate::assert_numeric(
    init_probs,
    min.len = 1,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    .var.name = "init_probs"
  )
  if (abs(sum(init_probs) - 1) > 1e-6) {
    stop("'init_probs' must sum to 1.")
  }
  invisible(NULL)
}

#' Validate a positive-integer sampler argument
#'
#' Used by [stan_options()] to coerce and check scalar count arguments such as
#' `iter`, `chains`, and `cores`.
#'
#' @param x the value to validate.
#' @param name the argument name, used in the error message.
#' @returns `x` coerced to a length-one integer.
#' @importFrom checkmate assert_integerish
#' @keywords internal
check_positive_int <- function(x, name) {
  checkmate::assert_integerish(
    x,
    lower = 1,
    len = 1,
    any.missing = FALSE,
    .var.name = name,
    coerce = TRUE
  )
}

#' Validate the `seed` sampler argument passed to [stan_options()]
#'
#' `seed` must be a single value coercible to an integer.
#'
#' @param seed the value to validate.
#' @returns `seed` coerced to a length-one integer.
#' @keywords internal
check_seed <- function(seed) {
  if (length(seed) != 1L) {
    stop("'seed' must be a single value.", call. = FALSE)
  }
  val <- suppressWarnings(as.integer(seed))
  if (is.na(val)) {
    stop("'seed' must be coercible to an integer.", call. = FALSE)
  }
  val
}
