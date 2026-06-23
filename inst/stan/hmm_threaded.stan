// reduce_sum-parallelised version of hmm.stan.
//
// The forward algorithm is household-separable: each household's likelihood is
// independent given the parameters, so households are sliced across threads via
// reduce_sum.
//
// Kept in inst/stan_threaded/ (not inst/stan/) on purpose: rstantools
// precompiles every inst/stan/*.stan as an rstan model at install time, which
// fails for this threaded model. It is instead compiled at runtime by cmdstanr.

functions {

  real get_diagonal_element(matrix m, int i) {
    real out = 1;
    for (j in 1:rows(m)) if (j != i) out -= m[j, i];
    return out;
  }

  int is_in(int pos, array[] int pos_var) {
    for (p in 1:size(pos_var)) if (pos_var[p] == pos) return 1;
    return 0;
  }

  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;
    for (i in 1:cols(m)) out[, i] = m[, i] / sum(m[, i]);
    return out;
  }

  matrix replace_zeroes(matrix m, real epsilon) {
    matrix[rows(m), cols(m)] out = m;
    for (i in 1:rows(m))
      for (j in 1:cols(m))
        if (m[i, j] == 0) out[i, j] = epsilon;
    return out;
  }

  // Summed forward-algorithm log-likelihood for households start:end.
  // slice_hh is required by reduce_sum's signature but unused; households are
  // indexed directly via start:end (hh_indices[h] == h).
  real partial_log_lik(
    array[] int slice_hh, int start, int end,
    // household layout
    array[] int hh_size, array[] int obs_per_hh,
    array[] int hh_start_ind, array[] int hh_end_ind,
    array[] int hh_tmin, array[] int hh_tmax,
    // observation data
    int n_obs_type, array[, ] int y, array[] int part_id, array[] int t_day,
    // transition structure
    int n_states, array[] int inf_states,
    int n_trans_fit, array[] int param_index,
    array[, ] int trans_index, array[, ] int source_states,
    int n_mult_fit, array[] int mult_param_index, array[, ] int mult_index,
    int n_inf_prob,
    matrix trans, matrix transition_multiplier,
    // fitted quantities (probability scale)
    array[] real params, array[] real mult_params,
    array[] real ih_prob, real eh_prob,
    // observation model + initialisation
    array[] matrix obs_prob, vector init_probs, real epsilon
  ) {

    real llik_sum = 0;
    int t_win = max(hh_tmax) - min(hh_tmin) + 1;

    for (h in start:end) {

      // Forward probabilities, local to this household (no cross-household state).
      matrix[hh_size[h] * n_states, t_win] alpha;
      matrix[hh_size[h] * n_states, t_win] logalpha;
      matrix[hh_size[h], t_win] llik = rep_matrix(0, hh_size[h], t_win);

      array[obs_per_hh[h], n_obs_type] int y_hh
        = y[hh_start_ind[h]:hh_end_ind[h], ];
      array[obs_per_hh[h]] int part_id_hh = part_id[hh_start_ind[h]:hh_end_ind[h]];
      array[obs_per_hh[h]] int t_day_hh = t_day[hh_start_ind[h]:hh_end_ind[h]];

      array[hh_size[h], n_states] int i_rows;
      int index = 1;
      int obs_switch;

      // ---- t = 1: initialise with starting probabilities ----
      for (i in 1:hh_size[h]) {
        array[n_states] int ref
          = linspaced_int_array(n_states, n_states*(i-1)+1, n_states*i);
        matrix[n_obs_type, n_states] obs;

        obs_switch = (t_day_hh[index] == 1 && part_id_hh[index] == i);
        if (obs_switch) {
          for (k in 1:n_obs_type) {
            if (y_hh[index, k] != -1)
              obs[k, ] = obs_prob[k][y_hh[index, k], ];
            else
              obs[k, ] = rep_row_vector(1, n_states);
          }
          index = min(index + 1, obs_per_hh[h]);
        } else {
          obs = rep_matrix(1, n_obs_type, n_states);
        }

        logalpha[ref, 1] = log(init_probs);
        for (k in 1:n_obs_type) logalpha[ref, 1] += to_vector(log(obs[k, ]));
        for (s in inf_states) i_rows[i, s] = n_states*(i-1) + s;

        llik[i, 1] = log_sum_exp(logalpha[ref, 1]);
        alpha[(n_states*(i-1)+1):(n_states*i), 1] = softmax(logalpha[ref, 1]);
      }

      // ---- t = 2, ..., household's last day: forward recursion ----
      for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {
        for (p in 1:hh_size[h]) {

          array[n_states] int ref
            = linspaced_int_array(n_states, n_states*(p-1)+1, n_states*p);
          vector[n_states] logalpha_temp = logalpha[ref, tt-1];
          matrix[n_obs_type, n_states] obs;
          // Rebuild the transition matrix from the base each step.
          matrix[n_states, n_states] trans_temp = trans;
          matrix[n_states, n_states] mult_temp = transition_multiplier;
          array[n_states] real no_inf_prob;
          matrix[hh_size[h], n_states] no_hh_inf_prob;

          obs_switch = (t_day_hh[index] == tt && part_id_hh[index] == p);
          if (obs_switch) {
            for (k in 1:n_obs_type) {
              if (y_hh[index, k] != -1)
                obs[k, ] = obs_prob[k][y_hh[index, k], ];
              else
                obs[k, ] = rep_row_vector(1, n_states);
            }
            index = min(index + 1, obs_per_hh[h]);
          } else {
            obs = rep_matrix(1, n_obs_type, n_states);
          }

          // Probability person p avoids infection from each household member.
          int ct = 1;
          for (s in 1:n_states) {
            if (is_in(s, inf_states)) {
              real p_ih = n_inf_prob == 1 ? ih_prob[1] : ih_prob[ct];
              no_hh_inf_prob[, s] =
                to_vector(alpha[i_rows[, s], tt-1]) * (1 - p_ih)
                + (1 - to_vector(alpha[i_rows[, s], tt-1]));
              no_hh_inf_prob[p, s] = 1;          // cannot infect self
              if (n_inf_prob != 1) ct += 1;
              no_inf_prob[s] = prod(no_hh_inf_prob[, s]);
            } else {
              no_inf_prob[s] = 1;
            }
          }

          // Fitted transitions: non-infection use a param, infection use the FoI.
          for (m in 1:n_trans_fit) {
            if (sum(source_states[m, ]) == 0) {
              trans_temp[trans_index[m, 1], trans_index[m, 2]] =
                params[param_index[m]];
            } else {
              real no_inf = 1;
              for (s in 1:n_states)
                if (source_states[m, s] == 1) no_inf *= no_inf_prob[s];
              trans_temp[trans_index[m, 1], trans_index[m, 2]] =
                1 - no_inf * (1 - eh_prob);
            }
          }

          // Fitted transition-split multipliers.
          for (m in 1:n_mult_fit) {
            if (mult_param_index[m] > 0)
              mult_temp[mult_index[m, 1], mult_index[m, 2]] =
                mult_params[mult_param_index[m]];
            else
              mult_temp[mult_index[m, 1], mult_index[m, 2]] +=
                -mult_params[-mult_param_index[m]];
          }
          trans_temp .*= mult_temp;

          // Make columns sum to 1, guard zeros, normalise.
          for (s in 1:n_states)
            trans_temp[s, s] = get_diagonal_element(trans_temp, s);
          trans_temp = normalize_cols(replace_zeroes(trans_temp, epsilon));

          logalpha[ref, tt] = log(trans_temp * exp(logalpha_temp));
          for (k in 1:n_obs_type) logalpha[ref, tt] += to_vector(log(obs[k, ]));

          alpha[(n_states*(p-1)+1):(n_states*p), tt] = softmax(logalpha[ref, tt]);
          llik[p, tt] = log_sum_exp(logalpha[ref, tt]);
        }
      }

      llik_sum += sum(llik[, hh_tmax[h] - hh_tmin[h] + 1]);
    }

    return llik_sum;
  }
}

data {

  // Transitions
  int n_states;
  matrix[n_states, n_states] trans;
  int n_inf_states;
  array[n_inf_states] int inf_states;
  int n_trans_fit;
  array[n_trans_fit] int param_index;
  array[n_trans_fit, 2] int trans_index;
  array[n_trans_fit, n_states] int source_states;
  int n_params;

  // Multipliers
  matrix[n_states, n_states] transition_multiplier;
  int n_mult_fit;
  int n_mult_params;
  array[n_mult_fit] int mult_param_index;
  array[n_mult_fit, 2] int mult_index;

  // Household information
  int n_hh;
  array[n_hh] int hh_size;

  // Data
  int n_obs;
  int n_obs_type;
  int n_unique_obs;
  array[n_obs, n_obs_type] int y;
  array[n_obs] int part_id;
  array[n_obs] int t_day;
  array[n_hh] int obs_per_hh;
  array[n_hh] int hh_start_ind;
  array[n_hh] int hh_end_ind;
  array[n_hh] int hh_tmin;
  array[n_hh] int hh_tmax;

  // Initial state and observation probabilities
  array[n_obs_type] matrix[n_unique_obs, n_states] obs_prob;
  vector[n_states] init_probs;

  real epsilon;
  int n_inf_prob;
}

parameters {
  array[n_params] real logit_params;
  array[n_mult_params] real logit_mult_params;
  real beta_eh;
  array[n_inf_prob] real beta_ih;
}

transformed parameters {
  // Only the cheap transforms live here; the forward algorithm runs inside
  // partial_log_lik so logalpha need not be stored.
  array[n_params] real params = inv_logit(logit_params);
  array[n_mult_params] real mult_params = inv_logit(logit_mult_params);
  array[n_inf_prob] real ih_prob = inv_logit(beta_ih);
  real eh_prob = inv_logit(beta_eh);
}

model {
  beta_eh ~ normal(-3, 3);
  beta_ih ~ normal(-3, 3);

  array[n_hh] int hh_indices;
  for (h in 1:n_hh) hh_indices[h] = h;

  target += reduce_sum(
    partial_log_lik, hh_indices, 1,
    hh_size, obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
    n_obs_type, y, part_id, t_day,
    n_states, inf_states, n_trans_fit, param_index, trans_index, source_states,
    n_mult_fit, mult_param_index, mult_index, n_inf_prob,
    trans, transition_multiplier,
    params, mult_params, ih_prob, eh_prob,
    obs_prob, init_probs, epsilon
  );
}

generated quantities {
  // Per-household log-likelihood (e.g. for LOO); replaces hmm.stan's llik_final.
  vector[n_hh] llik_final;
  for (h in 1:n_hh) {
    array[1] int hh_one = {h};
    llik_final[h] = partial_log_lik(
      hh_one, h, h,
      hh_size, obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
      n_obs_type, y, part_id, t_day,
      n_states, inf_states, n_trans_fit, param_index, trans_index, source_states,
      n_mult_fit, mult_param_index, mult_index, n_inf_prob,
      trans, transition_multiplier,
      params, mult_params, ih_prob, eh_prob,
      obs_prob, init_probs, epsilon
    );
  }
}
