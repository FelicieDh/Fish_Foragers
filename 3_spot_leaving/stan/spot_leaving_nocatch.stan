//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// 

functions {

  real partial_sum_lpmf(
    array[] int slice_n,
    int start, int end,
    array[] int trip_start, // row index where each trip starts
    array[] int trip_end, // row index where each trip ends
    array[] int leave,
    vector dist_sonar,
    vector fish3,
    array[] int compID,
    array[] int partID,
    array[] int tripID,
    vector beta,
    matrix v_comp,
    matrix v_id,
    matrix v_trip
  ) {
    real lp = 0;

    for (trip in start:end) {
            // loop through each timestamp of angling trip
    for (n in trip_start[trip]:trip_end[trip]) {
      real eta;

                eta = beta[1]+v_trip[tripID[n], 1]+v_comp[compID[n], 1] + v_id[partID[n], 1] +
                          (beta[2]+v_trip[tripID[n], 2]+v_comp[compID[n],2] + v_id[partID[n], 2])*fish3[n]+
                          (beta[3]+v_trip[tripID[n], 3]+v_comp[compID[n], 3] + v_id[partID[n], 3])*dist_sonar[n];
           

          lp += bernoulli_logit_lpmf( leave[n] | eta );

        }

    }

    return lp;
    }
}
data {
  int<lower = 1> N; // length of time series
  array[N] int<lower = 0, upper = 1> leave;

  
  // predictors
  vector[N] fish3; 
  //array[N] int dist_last_spot;
  vector[N] dist_sonar;
  // 

   // number of foraging trips
  int Ntrips;
  int Ncomp;
  int Nids;
  
  array[Ntrips] int slice_trips;
  array[N] int tripID;
  array[N] int compID;
  array[N] int partID;
  
  // row indices to subset trips
  array[Ntrips] int trip_start; // contains first index of trip
  array[Ntrips] int trip_end; // contains last index of trip
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // regression parameter
  vector[3] beta;


  // trip level offsets
  matrix[3, Ntrips] z_trip;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[3] sigma_trip;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[3] Rho_trip; //Cholesky factor for covariance of learning parameters among individuals
  
  // lake level offsets
  matrix[3, Ncomp] z_comp;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[3] sigma_comp;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[3] Rho_comp; //Cholesky factor for covariance of learning parameters among individuals
  
  // id level offsets
  matrix[3, Nids] z_id;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[3] sigma_id;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[3] Rho_id; //Cholesky factor for covariance of learning parameters among individuals
}

transformed parameters {

  // non centered parameterization
  matrix[Ntrips, 3] v_trip;
  v_trip = (diag_pre_multiply(sigma_trip, Rho_trip) * z_trip)';

  // non centered parameterization
  matrix[Ncomp, 3] v_comp;
  v_comp = (diag_pre_multiply(sigma_comp, Rho_comp) * z_comp)';

  // non centered parameterization
  matrix[Nids, 3] v_id;
  v_id = (diag_pre_multiply(sigma_id, Rho_id) * z_id)';

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.

model {

  // offsets
  to_vector(z_trip) ~ normal(0, 1);
  sigma_trip ~ exponential(1);
  Rho_trip ~ lkj_corr_cholesky(4);

  // offsets
  to_vector(z_comp) ~ normal(0, 1);
  sigma_comp ~ exponential(1);
  Rho_comp ~ lkj_corr_cholesky(4);

  //  offsets
  to_vector(z_id) ~ normal(0, 1);
  sigma_id ~ exponential(1);
  Rho_id ~ lkj_corr_cholesky(4);

  // regression priors
  beta ~ normal(0, 1);

  
  target += reduce_sum(partial_sum_lpmf,
                      slice_trips,
                      1,
                      trip_start, // row index where each trip starts
                      trip_end, // row index where each trip ends
                      leave,
                      dist_sonar,
                      fish3,
                      compID,
                      partID,
                      tripID,
                      beta,
                      v_comp,
                      v_id,
                      v_trip

  );
            
                        
}
generated quantities {
  
  vector[N] log_lik;
  
   for (trip in 1:Ntrips) {
            // loop through each timestamp of angling trip
    for (n in trip_start[trip]:trip_end[trip]) {
    real eta;

                eta = beta[1]+v_trip[tripID[n], 1]+v_comp[compID[n], 1] + v_id[partID[n], 1] +
                          (beta[2]+v_trip[tripID[n], 2]+v_comp[compID[n],2] + v_id[partID[n], 2])*fish3[n]+
                          (beta[3]+v_trip[tripID[n], 3]+v_comp[compID[n], 3] + v_id[partID[n], 3])*dist_sonar[n];
                          
          log_lik[n] = bernoulli_logit_lpmf( leave[n] | eta );

        }

    }
  
}



