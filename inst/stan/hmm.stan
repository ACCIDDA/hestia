
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
  int is_in(int pos,int[] pos_var) {
    int pos_match;
    int all_matches[size(pos_var)];
    
    for (p in 1:(size(pos_var))) {
      all_matches[p] = (pos_var[p]==pos);
    }
    
    if(sum(all_matches)>0) {
      pos_match = 1;
      return pos_match;
    } else {
      pos_match = 0;
      return pos_match;
    }
  }
  
  // Normailze the columns of a matrix to sum to 1
  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;
  
    for(i in 1:cols(m)) {
      out[,i] = m[,i]/sum(m[,1]);
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
  
}


data {
  
  // Transitions
  int  n_states; // number of latent states
  matrix[n_states, n_states] trans; // transition probabilities, cols = starting state, rows = ending states
  int n_inf_states; // number of infectious states
  int inf_states[n_inf_states]; // infectious states
  int n_trans_fit; // number of transitions to fit
  int param_index[n_trans_fit]; // parameter corresponding with each non-infection transition to fit, 0 if an infection transition
  int trans_index[n_trans_fit, 2]; // row/col indices of transition matrix corresponding to each parameter to be fit
  int source_states[n_trans_fit, n_states]; // states that are the source of infecetion of the transition (0 if non-infection transition)
  int n_params; // number of additional (non-infection) parameters to fit
  
  // Multipliers
  matrix[n_states, n_states] multiplier; // transition multiplier - allows for transition splits
  int n_mult_fit; // number of multipliers to fit
  int n_mult_params; // number of unique multipliers parameters to fit
  int mult_param_index[n_mult_fit]; // parameter corresponding with each non-infection transition to fit, negative if a 1- situation
  int mult_index[n_mult_fit, 2]; // row/col indices of transition matrix corresponding to each parameter to be fit
  
  // Household information
  int n_hh; // number of households
  int hh_size[n_hh]; // household size
    
  // Data 
  int n_obs; // number of observations
  int n_obs_type; // number of observation types
  int n_unique_obs; // number of unique outcomes for enrolled individuals
  int y[n_obs, n_obs_type]; // outcome vector, ordered by 1) household, then 2) time, then 3) individual
  int part_id[n_obs]; // particpants associated with the observation
  int t_day[n_obs]; // observation times for enrolled (day)
  int obs_per_hh[n_hh]; // total number of observations per HH for enrolled
  int hh_start_ind[n_hh]; // starting index for the HH
  int hh_end_ind[n_hh]; // ending index for the HH
  int hh_tmin[n_hh]; // minimum day for which there is a HH observation
  int hh_tmax[n_hh]; // maximum day for which there is a HH observation

  // Initial state and observation probabilities
  matrix[n_unique_obs, n_states] obs_prob[n_obs_type]; // observation process for SIR states
  vector[n_states] init_probs; // starting state probabilities
  
  real epsilon; // small number to avoid log(0) issues
  int n_inf_prob; // number of infection probabilties to fit, either 1 or equal to the number of infectious compartments
}

parameters {
  real logit_params[n_params];
  real logit_mult_params[n_mult_params];
  real beta_eh; // monthly intercepts for extra-household probabilities
  real beta_ih[n_inf_prob]; // intra-household probability

}

transformed parameters {
  vector[n_hh] llik_final; // sum of logalpha for final timestep
  real ih_prob[n_inf_prob];
  real eh_prob;
  matrix[sum(hh_size)*n_states, max(hh_tmax)-min(hh_tmin) + 1] logalpha; // log forward probability
  matrix[n_states, n_states] trans_temp;
  real params[n_params];
  real mult_params[n_mult_params];
  
  params = inv_logit(logit_params);
  mult_params = inv_logit(logit_mult_params);
  ih_prob = inv_logit(beta_ih);
  eh_prob = inv_logit(beta_eh);
  trans_temp = trans;
  
  
  for(h in 1:n_hh) { // loop through household
    
    matrix[hh_size[h]*n_states, max(hh_tmax)-min(hh_tmin) + 1] alpha; // forward prob, normalized
    matrix[hh_size[h], max(hh_tmax)-min(hh_tmin) + 1] llik; // lik contribution for enrolled per participant and time
    int y_hh[obs_per_hh[h], n_obs_type];
    int part_id_hh[obs_per_hh[h]];
    int t_day_hh[obs_per_hh[h]];
    int index; // index for next observation
    int i_rows[hh_size[h], n_states];// rows in alpha corresponding to infectious states
    int last_lik;
    int obs_switch; // indicator for whether there is an observation corresponding to this time step
    
    llik = rep_matrix(0, hh_size[h], max(hh_tmax)-min(hh_tmin) + 1);
    
    if(h == 1) {
      last_lik = 0;
    } else {
      last_lik = sum(hh_size[1:(h-1)]);
    }
    
    // subset to data only for the given HH
    y_hh = y[(hh_start_ind[h]):(hh_end_ind[h]),];
    t_day_hh = t_day[(hh_start_ind[h]):(hh_end_ind[h])];
    part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])];
    
    index = 1;
    
    { // START FORWARD ALGORITHM
    
    // fill first column of alpha using starting probabilities
    for(i in 1:hh_size[h]) {
      int ref[n_states];
      matrix[n_obs_type, n_states] obs; // observation component for enrolled memebrs, set to 1 if no observation for this time step
    
      ref = linspaced_int_array(n_states, n_states*last_lik+n_states*(i-1)+1, n_states*last_lik+n_states*(i-1)+n_states);
      
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
      
      llik[i, 1] = log_sum_exp(logalpha[ref,1]);
      
      // normalize and convert to the probability scale
      alpha[(n_states*(i-1)+1):(n_states*(i-1)+n_states), 1] = softmax(logalpha[ref,1]);
    
    } // end participant loop - t=1, update logalpha with observation probability
    for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {
      
      for(p in 1:hh_size[h]) {
        real no_inf_prob[n_states]; // probability of avoiding all infections
        matrix[hh_size[h], n_states] no_hh_inf_prob; // probability of avoiding infection from each HH member
        int ref[n_states];
        vector[n_states] logalpha_temp; // log forward probability
        matrix[n_obs_type, n_states] obs;
        matrix[n_states, n_states] mult_temp;
        
        ref = linspaced_int_array(n_states, n_states*last_lik+n_states*(p-1)+1, n_states*last_lik+n_states*(p-1)+n_states);
        
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
        mult_temp = multiplier; // need to reset mult_temp since it is self-referential
        for(m in 1:n_mult_fit) {
          if(mult_param_index[m] > 0) {
            mult_temp[mult_index[m, 1],mult_index[m, 2]] = mult_params[mult_param_index[m]];
          } else {
            mult_temp[mult_index[m, 1], mult_index[m, 2]] += -1*mult_params[-mult_param_index[m]];
          }
        }
        
        // transition splits
        trans_temp = trans_temp .* mult_temp;

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
        
        llik[p, tt] = log_sum_exp(logalpha[ref,tt]);
        
      } // end participant loop - update logalpha with observation probability
      
      if(tt == (hh_tmax[h] - hh_tmin[h] + 1)) {
        llik_final[h] = sum(llik[,tt]); 
      }
        
    } // end time loop
    } // END FORWARD ALGORITHM
    
  } // end household loop
  
}

model {
  
  beta_eh ~ normal(-3,3);
  beta_ih ~ normal(-3,3);
  
  // Only increment by final alpha
  target += sum(llik_final);

}

