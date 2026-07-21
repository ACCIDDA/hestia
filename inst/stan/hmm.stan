
functions {

  // Calculate diagonal element i of matrix m
  // Diagonal element = 1 - sum(non-diagonal elements)
  real get_diagonal_element(matrix m, int i){
    real out;
    out = 1;
    for(j in 1:rows(m)) {
      if(j != i) {
        out = out - m[j, i];
      }
    }
    return out;
  }

  // function that is comparable to R's %in% function
  // pos is the value to test for matches
  // array pos_var is the 1-dimensional array of possible matches
  // returns pos_match=1 if pos matches at least one element of pos_var and pos_match=0 otherwise
  // example code:
  // if(r_in(3,{1,2,3,4})) will evaluate as TRUE
  // from: https://discourse.mc-stan.org/t/stan-equivalent-of-rs-in-function/3849

  int is_in(int pos, array[] int pos_var) {
    int pos_match;
    array[size(pos_var)] int all_matches;

    for (p in 1:size(pos_var)) {
      all_matches[p] = (pos_var[p] == pos);
    }

    if (sum(all_matches) > 0) return 1;
    else return 0;
  }

  // Normailze the columns of a matrix to sum to 1
  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;

    for(i in 1:cols(m)) {
      out[,i] = m[,i]/sum(m[,i]);
    }
    return out;
  }

  matrix replace_zeroes(matrix m, real epsilon) {

    matrix[rows(m), cols(m)] out;
    out = m;
    for(i in 1:rows(m)) {
      for(j in 1:cols(m)) {
        if(m[i,j] == 0) {
          out[i,j] = epsilon;
        }
      }
    }
    return(out);
}

  // Per-household forward algorithm returning the log forward-probability matrix
  // (logalpha) for household h.
  matrix hh_logalpha(int h,
                     array[] int hh_size,
                     int n_obs_type,
                     int n_states,
                     array[,] int y,
                     array[] int part_id,
                     array[] int t_day,
                     array[] int obs_per_hh,
                     array[] int hh_start_ind,
                     array[] int hh_end_ind,
                     array[] int hh_tmin,
                     array[] int hh_tmax,
                     array[] int inf_states,
                     int n_inf_prob,
                     int n_trans_fit,
                     array[] int param_index,
                     array[,] int trans_index,
                     array[,] int source_states,
                     int n_mult_fit,
                     array[] int mult_param_index,
                     array[,] int mult_index,
                     matrix trans,
                     matrix transition_multiplier,
                     array[] real params,
                     array[] real mult_params,
                     array[] real ih_prob,
                     real eh_prob,
                     array[] matrix obs_prob,
                     vector init_probs,
                     real epsilon) {

    matrix[hh_size[h]*n_states, max(hh_tmax)-min(hh_tmin) + 1] alpha; // forward prob, normalized
    matrix[hh_size[h]*n_states, max(hh_tmax)-min(hh_tmin) + 1] logalpha; // log forward probability
    array[obs_per_hh[h], n_obs_type] int y_hh;
    array[obs_per_hh[h]] int part_id_hh;
    array[obs_per_hh[h]] int t_day_hh;
    int index; // index for next observation
    array[hh_size[h], n_states] int i_rows; // rows in alpha corresponding to infectious states
    int obs_switch; // indicator for whether there is an observation corresponding to this time step

    // subset to data only for the given HH
    y_hh = y[(hh_start_ind[h]):(hh_end_ind[h]),];
    t_day_hh = t_day[(hh_start_ind[h]):(hh_end_ind[h])];
    part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])];

    index = 1;

    { // START FORWARD ALGORITHM

    // fill first column of alpha using starting probabilities
    for(i in 1:hh_size[h]) {
      array[n_states] int ref;
      matrix[n_obs_type, n_states] obs; // observation component for enrolled memebrs, set to 1 if no observation for this time step

      ref = linspaced_int_array(n_states, n_states*(i-1)+1, n_states*(i-1)+n_states);

      obs_switch = 0;

      if(t_day_hh[index] == 1) {
        if(part_id_hh[index] == i) {
          obs_switch = 1;
        }
      }

      if(obs_switch == 1) {
        for(k in 1:n_obs_type) {
          if(y_hh[index, k] != -1) {
            obs[k, ] = obs_prob[k, y_hh[index, k],];
          } else {
            obs[k, ] = rep_row_vector(1, n_states);
          }
        }
      } else {
        obs = rep_matrix(1, n_obs_type, n_states);
      }

      if(obs_switch == 1) {
       index = min(index + 1, hh_end_ind[h]-hh_start_ind[h]+1);
      }

      // Fill in starting probability for SIR states
      logalpha[ref, 1] = log(init_probs);
      for(k in 1:n_obs_type) {
        logalpha[ref, 1] = logalpha[ref, 1] + to_vector(log(obs[k,]));
      }
      for(s in inf_states) {
        i_rows[i, s] = n_states*(i-1)+s;
      }

      // normalize and convert to the probability scale
      alpha[(n_states*(i-1)+1):(n_states*(i-1)+n_states), 1] = softmax(logalpha[ref,1]);

    } // end participant loop - t=1, update logalpha with observation probability
    for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {

      for(p in 1:hh_size[h]) {
        array[n_states] real no_inf_prob; // probability of avoiding all infections
        matrix[hh_size[h], n_states] no_hh_inf_prob; // probability of avoiding infection from each HH member
        array[n_states] int ref;
        vector[n_states] logalpha_temp; // log forward probability
        matrix[n_obs_type, n_states] obs;
        matrix[n_states, n_states] trans_temp = trans; // rebuild from the base each step
        matrix[n_states, n_states] mult_temp;

        ref = linspaced_int_array(n_states, n_states*(p-1)+1, n_states*(p-1)+n_states);

        logalpha_temp = logalpha[ref,tt-1];

        obs_switch = 0;

        if(t_day_hh[index] == tt) {
          if(part_id_hh[index] == p) {
            obs_switch = 1;
          }
        }

        if(obs_switch == 1) {
          for(k in 1:n_obs_type) {
            if(y_hh[index, k] != -1) {
              obs[k,] = obs_prob[k, y_hh[index, k], ];
            } else {
              obs[k,] = rep_row_vector(1, n_states);
            }
          }
        } else {
          obs = rep_matrix(1, n_obs_type, n_states);
        }

        if(obs_switch == 1) {
          index = min(index + 1, hh_end_ind[h]-hh_start_ind[h]+1);
        }

        int ct = 1;
        for(s in 1:n_states) {
          if(is_in(s, inf_states)) {
            if(n_inf_prob == 1) {
              no_hh_inf_prob[,s] = to_vector(alpha[i_rows[, s], tt-1])*(1-ih_prob[1]) + (1 - to_vector(alpha[i_rows[,s], tt-1])); // Pr of avoiding infection from each household member
            } else {
              no_hh_inf_prob[,s] = to_vector(alpha[i_rows[, s], tt-1])*(1-ih_prob[ct]) + (1 - to_vector(alpha[i_rows[,s], tt-1])); // Pr of avoiding infection from each household member
              ct += 1;
            }
            no_hh_inf_prob[p, s] = 1; // Particpant can't infect themselves
          } else {
            no_hh_inf_prob[,s] = rep_vector(1, hh_size[h]);
          }
          no_inf_prob[s] = prod(no_hh_inf_prob[,s]); // Probability of avoiding infection from all household members
        }

        // fill in tranistions that are being fit
        for(m in 1:n_trans_fit) {
          if(sum(source_states[m,]) == 0) {
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = params[param_index[m]];
          } else {
            real no_inf = 1;
            for(s in 1:n_states) {
              if(source_states[m,s] == 1) {
                no_inf = no_inf*no_inf_prob[s];
              }
            }
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = 1-(no_inf*(1-eh_prob));
          }
        }

        // fill in multipliers that are being fit
        mult_temp = transition_multiplier; // need to reset mult_temp since it is self-referential
        for(m in 1:n_mult_fit) {
          if(mult_param_index[m] > 0) {
            mult_temp[mult_index[m, 1],mult_index[m, 2]] = mult_params[mult_param_index[m]];
          } else {
            mult_temp[mult_index[m, 1], mult_index[m, 2]] += -1*mult_params[-mult_param_index[m]];
          }
        }

        // transition splits
        trans_temp .*= mult_temp;

        // fill in diagonals (columns must sum to one)
        for(i in 1:cols(trans_temp)) {
          trans_temp[i,i] = get_diagonal_element(trans_temp, i);
        }

        // replace zeroes with epsilon and normalize
        trans_temp = replace_zeroes(trans_temp, epsilon);
        trans_temp = normalize_cols(trans_temp);

        // Compute the probability of each epidemiological state
        logalpha[ref, tt] = log(trans_temp*exp(logalpha_temp));
        for(k in 1:n_obs_type) {
          logalpha[ref, tt] = logalpha[ref, tt] + to_vector(log(obs[k,]));
        }

        // normalize and convert to probability scale
        alpha[(n_states*(p-1)+1):(n_states*(p-1)+n_states), tt] = softmax(logalpha[ref,tt]);

      } // end participant loop - update logalpha with observation probability

    } // end time loop
    } // END FORWARD ALGORITHM

    return logalpha;
  }

  // Household log-likelihood: sum over participants of the final-column
  // log_sum_exp of logalpha (household h's forward-algorithm contribution).
  real hh_llik(matrix logalpha, int hh_size_h, int n_states, int Th) {
    real s = 0;
    for(p in 1:hh_size_h) {
      array[n_states] int ref = linspaced_int_array(n_states, n_states*(p-1)+1, n_states*(p-1)+n_states);
      s += log_sum_exp(logalpha[ref, Th]);
    }
    return s;
  }

  // Summed forward-algorithm log-likelihood for households start:end. Body is
  // hmm.stan's transformed-parameters household loop; slice_hh is required by
  // the reduce_sum signature but unused (households indexed via start:end).
  real partial_log_lik(array[] int slice_hh, int start, int end,
                       array[] int hh_size,
                       int n_obs_type, 
                       int n_states,
                       array[,] int y, 
                       array[] int part_id, 
                       array[] int t_day,
                       array[] int obs_per_hh,
                       array[] int hh_start_ind, 
                       array[] int hh_end_ind,
                       array[] int hh_tmin, 
                       array[] int hh_tmax,
                       array[] int inf_states, 
                       int n_inf_prob,
                       int n_trans_fit, 
                       array[] int param_index,
                       array[,] int trans_index, 
                       array[,] int source_states,
                       int n_mult_fit, 
                       array[] int mult_param_index,
                       array[,] int mult_index,
                       matrix trans, 
                       matrix transition_multiplier,
                       array[] real params, 
                       array[] real mult_params,
                       array[] real ih_prob, 
                       real eh_prob,
                       array[] matrix obs_prob, 
                       vector init_probs, 
                       real epsilon) {

  real llik_sum = 0;

  for(h in start:end) {
    matrix[hh_size[h]*n_states, max(hh_tmax)-min(hh_tmin) + 1] logalpha =
      hh_logalpha(h, hh_size, n_obs_type, n_states, y, part_id, t_day,
                  obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
                  inf_states, n_inf_prob, n_trans_fit, param_index, trans_index,
                  source_states, n_mult_fit, mult_param_index, mult_index,
                  trans, transition_multiplier, params, mult_params,
                  ih_prob, eh_prob, obs_prob, init_probs, epsilon);
    llik_sum += hh_llik(logalpha, hh_size[h], n_states, hh_tmax[h] - hh_tmin[h] + 1);
  }

  return llik_sum;
  }

}


data {

  // Transitions
  int  n_states; // number of latent states
  matrix[n_states, n_states] trans; // transition probabilities, cols = starting state, rows = ending states
  int n_inf_states; // number of infectious states
  array[n_inf_states] int inf_states; // infectious states
  int n_trans_fit; // number of transitions to fit
  array[n_trans_fit] int param_index; // parameter corresponding with each non-infection transition to fit, 0 if an infection transition
  array[n_trans_fit, 2] int trans_index; // row/col indices of transition matrix corresponding to each parameter to be fit
  array[n_trans_fit, n_states] int source_states; // states that are the source of infecetion of the transition (0 if non-infection transition)
  int n_params; // number of additional (non-infection) parameters to fit

  // Multipliers
  matrix[n_states, n_states] transition_multiplier; // transition multiplier - allows for transition splits
  int n_mult_fit; // number of multipliers to fit
  int n_mult_params; // number of unique multipliers parameters to fit
  array[n_mult_fit] int mult_param_index; // parameter corresponding with each non-infection transition to fit, negative if a 1- situation
  array[n_mult_fit, 2] int mult_index; // row/col indices of transition matrix corresponding to each parameter to be fit

  // Household information
  int n_hh; // number of households
  array[n_hh] int hh_size; // household size

  // Data
  int n_obs; // number of observations
  int n_obs_type; // number of observation types
  int n_unique_obs; // number of unique outcomes for enrolled individuals
  array[n_obs, n_obs_type] int y; // outcome vector, ordered by 1) household, then 2) time, then 3) individual
  array[n_obs] int part_id; // particpants associated with the observation
  array[n_obs] int t_day; // observation times for enrolled (day)
  array[n_hh] int obs_per_hh; // total number of observations per HH for enrolled
  array[n_hh] int hh_start_ind; // starting index for the HH
  array[n_hh] int hh_end_ind; // ending index for the HH
  array[n_hh] int hh_tmin; // minimum day for which there is a HH observation
  array[n_hh] int hh_tmax; // maximum day for which there is a HH observation

  // Initial state and observation probabilities
  array[n_obs_type] matrix[n_unique_obs, n_states] obs_prob; // observation process for SIR states
  vector[n_states] init_probs; // starting state probabilities

  real epsilon; // small number to avoid log(0) issues
  int n_inf_prob; // number of infection probabilties to fit, either 1 or equal to the number of infectious compartments

  int<lower=0, upper=1> save_llik; // if 1, write per-household log-likelihood (llik_final) in generated quantities
  int<lower=0, upper=1> save_states; // if 1, write per-participant log forward probabilities (logalpha) in generated quantities
}

parameters {
  array[n_params] real logit_params;
  array[n_mult_params] real logit_mult_params;
  real beta_eh; // monthly intercepts for extra-household probabilities
  array[n_inf_prob] real beta_ih; // intra-household probability

}

transformed parameters {
  array[n_inf_prob] real ih_prob;
  real eh_prob;
  array[n_params] real params;
  array[n_mult_params] real mult_params;

  params = inv_logit(logit_params);
  mult_params = inv_logit(logit_mult_params);
  ih_prob = inv_logit(beta_ih);
  eh_prob = inv_logit(beta_eh);
}

model {

  beta_eh ~ normal(-3,3);
  beta_ih ~ normal(-3,3);

  // Parallelised forward algorithm (households are independent given params)
  array[n_hh] int hh_indices;
  for(h in 1:n_hh) hh_indices[h] = h;

  target += reduce_sum(partial_log_lik, hh_indices, 1,
                       hh_size,
                       n_obs_type, n_states,
                       y, part_id, t_day,
                       obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
                       inf_states, n_inf_prob,
                       n_trans_fit, param_index, trans_index, source_states,
                       n_mult_fit, mult_param_index, mult_index,
                       trans, transition_multiplier,
                       params, mult_params, ih_prob, eh_prob,
                       obs_prob, init_probs, epsilon);

}

generated quantities {
  // Per-household log-likelihood for loo/waic; empty (length 0) unless save_llik = 1.
  vector[save_llik ? n_hh : 0] llik_final;
  
  // Per-participant log forward probabilities; empty (0x0) unless save_states = 1.
  matrix[save_states ? sum(hh_size)*n_states : 0, save_states ? max(hh_tmax)-min(hh_tmin)+1 : 0] logalpha;
  
  if(save_llik || save_states) {
    int row_offset = 0;
    if(save_states) logalpha = rep_matrix(0, rows(logalpha), cols(logalpha));
    for(h in 1:n_hh) {
      int Th = hh_tmax[h] - hh_tmin[h] + 1;
      matrix[hh_size[h]*n_states, max(hh_tmax)-min(hh_tmin)+1] la =
        hh_logalpha(h, hh_size, n_obs_type, n_states, y, part_id, t_day,
                    obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
                    inf_states, n_inf_prob, n_trans_fit, param_index, trans_index,
                    source_states, n_mult_fit, mult_param_index, mult_index,
                    trans, transition_multiplier, params, mult_params,
                    ih_prob, eh_prob, obs_prob, init_probs, epsilon);
      if(save_llik) llik_final[h] = hh_llik(la, hh_size[h], n_states, Th);
      if(save_states) logalpha[(row_offset+1):(row_offset+hh_size[h]*n_states), 1:Th] = la[, 1:Th];
      row_offset += hh_size[h]*n_states;
    }
  }
}
