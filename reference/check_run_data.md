# Validate the observation `data` frame passed to `run_model()`

Validate the observation `data` frame passed to
[`run_model()`](https://accidda.github.io/hestia/reference/run_model.md)

## Usage

``` r
check_run_data(data, obs_model)
```

## Arguments

- data:

  the participant-level data frame.

- obs_model:

  output of
  [`make_observation_model()`](https://accidda.github.io/hestia/reference/make_observation_model.md);
  its names must each appear as an outcome column in `data`.
