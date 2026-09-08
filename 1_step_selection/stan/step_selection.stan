//
// Full model
// 

functions {
  real partial_sum_lpmf(
    array[] int slice_n,
    int start, int end,
    int N_spots,
    array[] real time_since_start,
    array[] vector locality,
    array[] vector sonar_distance,
    array[] int columns_success,
    int max_columns_success,
    array[] matrix spatial_success,
    array[] matrix time_success,
    array[] int columns_loss,
    int max_columns_loss,
    array[] matrix spatial_loss,
    array[] matrix time_loss,
    array[] int columns_social,
    int max_columns_social,
    array[] matrix distance_social,
    array[] int compID,
    array[] int partID,
    array[] int tripID,
    vector beta,
    vector bandwidth_log,
    matrix v_compID,
    matrix v_partID,
    matrix v_tripID
  ) {
    real lp = 0;
    for (n in start:end) {

      vector[N_spots] social_feature;

      // Social kernel (success)
      social_feature = (
        exp(- (distance_social[n][, 1:columns_social[n]] .^2) /
          (2 * square(exp(bandwidth_log[1] + v_compID[compID[n],7]))))
            .* (1.0 ./ (1 + exp(-20 * (distance_social[n][, 1:columns_social[n]] - 5)))) //5m buffer
      ) * rep_vector(1.0, columns_social[n]);
      

      // Success kernel
      

    matrix[N_spots,columns_success[n]] success_spatial_feature;
    matrix[N_spots,columns_success[n]] success_time_feature;

    vector[N_spots] success_feature;
   
     //compute long social memory
     if (spatial_success[n][1, 1]<0){
        success_feature =rep_vector(0, N_spots);
     } else {
       success_spatial_feature = (
        exp(
          - (spatial_success[n][, 1:columns_success[n]] .^2) / (2 * square(exp(bandwidth_log[2] + v_compID[compID[n],8] )))
        )
      ); 
      
        success_time_feature = (
        exp(
          - (time_success[n][, 1:columns_success[n]] .^2) / (2 * square(exp(bandwidth_log[3] + v_compID[compID[n],9])))
        )
      ); 
      
      success_feature = ( success_spatial_feature .*  success_time_feature)* rep_vector(1.0, columns_success[n]); // element wise mulltiplication is a dot
      
      };
      

        // Loss kernel


    matrix[N_spots,columns_loss[n]] loss_spatial_feature;
    matrix[N_spots,columns_loss[n]] loss_time_feature;

    vector[N_spots] loss_feature;
   
     //compute long social memory
     if (spatial_loss[n][1, 1]<0){
        loss_feature =rep_vector(0, N_spots);
     } else {
       loss_spatial_feature = (
        exp(
          - (spatial_loss[n][, 1:columns_loss[n]] .^2) / (2 * square(exp(bandwidth_log[4] + v_compID[compID[n],10] )))
        )
      ); 
      
        loss_time_feature = (
        exp(
          - (time_loss[n][, 1:columns_loss[n]] .^2) / (2 * square(exp(bandwidth_log[5] + v_compID[compID[n],11])))
        )
      ); 
      
      loss_feature = ( loss_spatial_feature .*  loss_time_feature)* rep_vector(1.0, columns_loss[n]); // element wise mulltiplication is a dot
      
      };
      

      vector[N_spots] eta;
      eta =
         (beta[1] + v_compID[compID[n],1] + v_partID[partID[n],1]  + v_tripID[tripID[n],1]) * social_feature
        + (beta[2] + v_compID[compID[n],2] + v_partID[partID[n],2] + v_tripID[tripID[n],2]) * success_feature
        + (beta[3] + v_compID[compID[n],3] + v_partID[partID[n],3] + v_tripID[tripID[n],3]) * loss_feature
        + (beta[4] + v_compID[compID[n],4] + v_partID[partID[n],4] + v_tripID[tripID[n],4]) * sonar_distance[n]
        + (beta[5] + v_compID[compID[n],5] + v_partID[partID[n],5] + v_tripID[tripID[n],5]) * locality[n]
        + (beta[6] + v_compID[compID[n],6] + v_partID[partID[n],6] + v_tripID[tripID[n],6]) * locality[n]*time_since_start[n];

      lp += categorical_logit_lpmf(1 | eta);
    }
    return lp;
  }
}


data {
  int<lower=1> N; // number of patch choices observed
  int<lower=1> N_spots; // number of options per choice (including the chosen one and the synthetic atlernatives)

  array[N] real<lower=0, upper=1> time_since_start; //time variable
  array[N] vector[N_spots] locality; //locality matrix
  array[N] vector[N_spots] sonar_distance; //edge distance matrix
  
  
  //the success matrix
  array[N] int columns_success; //how many columns to subset from distance matrix
  int<lower=1> max_columns_success;
  array[N] matrix[N_spots, max_columns_success] spatial_success;
  array[N] matrix[N_spots, max_columns_success] time_success;


  //the loss matrix
  array[N] int columns_loss; //how many columns to subset from distance matrix
  int<lower=1> max_columns_loss;
  array[N] matrix[N_spots, max_columns_loss] spatial_loss;
  array[N] matrix[N_spots, max_columns_loss] time_loss;

  //the social matrix data 
  array[N] int columns_social; //number of columns in the social dataset
  int<lower=1> max_columns_social;
  array[N] matrix[N_spots, max_columns_social] distance_social;

  
  //the compID random offset
  int N_compID; //this is the number of competition ID
  array[N] int compID; //this is the ID of compID (I need to change my compID to integer in my data prep step)
  
  //the anglerID random offset
  int N_partID; //this is the number of participant ID
  array[N] int partID; //this is the ID of participant ID (I need to change my compID to integer in my data prep step)

    //the tripID random offset
  int N_tripID; //this is the number of participant ID
  array[N] int tripID; //this is the ID of participant ID (I need to change my compID to integer in my data prep step)
}

parameters {
  //real intercept;

  // Fixed effects
  vector[6] beta; //one of these effects if the bandwidth
  
  vector[5] bandwidth_log; //a bandwidth for each calculated feature
  
  matrix[11, N_compID] z_compID; //this creates the matrix for the random offsets that I need
  vector<lower = 0>[11] sigma_compID;  //Here I have 1 sigma for each of my parameters for the compID facotr
  cholesky_factor_corr[11] Rho_compID; //Cholesky factor for covariance of learning parameters among compID, with 7 parameters (it's a 7x6 matrix)
  
  matrix[6, N_partID] z_partID; //this creates the matrix for the random offsets of participant that I need
  vector<lower = 0>[6] sigma_partID;  //Here I have 1 sigma for each of my parameters for the compID facotr
  cholesky_factor_corr[6] Rho_partID; //Cholesky factor for covariance of learning parameters among compID, with 7 parameters (it's a 7x6 matrix)

  matrix[6, N_tripID] z_tripID; //this creates the matrix for the random offsets of participant that I need
  vector<lower = 0>[6] sigma_tripID;  //Here I have 1 sigma for each of my parameters for the compID facotr
  cholesky_factor_corr[6] Rho_tripID; //Cholesky factor for covariance of learning parameters among compID, with 7 parameters (it's a 7x6 matrix)
}

transformed parameters { //I use this to transfmor my z_compID and others, out of my sigma and rho
  // non centered parameterization
  matrix[N_compID, 11] v_compID; //7 parameters
  v_compID = (diag_pre_multiply(sigma_compID, Rho_compID) * z_compID)'; //this gives me my values on the correct scale
  
    // non centered parameterization
  matrix[N_partID, 6] v_partID; //7 parameters
  v_partID = (diag_pre_multiply(sigma_partID, Rho_partID) * z_partID)'; //this gives me my values on the correct scale
  
      // non centered parameterization
  matrix[N_tripID, 6] v_tripID; //7 parameters
  v_tripID = (diag_pre_multiply(sigma_tripID, Rho_tripID) * z_tripID)'; //this gives me my values on the correct scale
  
  //I can also pre transform my bandwith here if I want instead of in the model block
  
}

model {

  bandwidth_log ~ normal(3, 1); //(3,1)
  beta ~ normal(0, 2); 
  
  to_vector(z_compID) ~ normal(0, 1); 
  sigma_compID ~ exponential(1);  //variances
  Rho_compID~ lkj_corr_cholesky(4); //How much weight for extreme correlations. 4 seems reasonable in Schakowski et al. 2025

  to_vector(z_partID) ~ normal(0, 1); 
  sigma_partID ~ exponential(1);  //variances
  Rho_partID~ lkj_corr_cholesky(4); //How much weight for extreme correlations. 4 seems reasonable in Schakowski et al. 2025

  to_vector(z_tripID) ~ normal(0, 1); 
  sigma_tripID ~ exponential(0.5);  //variances
  Rho_tripID~ lkj_corr_cholesky(4); //How much weight for extreme correlations. 4 seems reasonable in Schakowski et al. 2025

 
 target += reduce_sum(
  partial_sum_lpmf,                // function
  rep_array(1, N),                // indexes 1..N (or an actual index array if needed)
  1,  // grainsize: tweak this for performance (e.g., 50 or 100)
  N_spots,
  time_since_start,
  locality,
  sonar_distance,
  columns_success,
  max_columns_success,
  spatial_success,
  time_success,
  columns_loss,
  max_columns_loss,
  spatial_loss,
  time_loss,
  columns_social,
  max_columns_social,
  distance_social,
  compID,
  partID,
  tripID,
  beta,
  bandwidth_log,
  v_compID,
  v_partID,
  v_tripID
);

  }

//generated quantities {
//   vector[N] log_lik;
//   //real bandwidth = 10.0;
// 
//   for (n in 1:N) {
//      vector[N_spots] social_feature;
// 
// 
//       // Social kernel
//       social_feature = (
//         exp(- (Distance_social[n][, 1:columns_social[n]] .^2) /
//           (2 * square(exp(bandwidth_log[1] + v_compID[compID[n],7]))))
//         .* (1.0 ./ (1 + exp(-20 * (Distance_social[n][, 1:columns_social[n]] - 5))))
//       ) * rep_vector(1.0, columns_social[n]);
// 
//       // Success kernel
//       
//     matrix[N_spots,columns_success[n]] success_spatial_feature;
//     matrix[N_spots,columns_success[n]] success_time_feature;
// 
//     vector[N_spots] success_feature;
//    
//      //compute long social memory
//      if (spatial_success[n][1, 1]<0){
//         success_feature =rep_vector(0, N_spots);
//      } else {
//        success_spatial_feature = (
//         exp(
//           - (spatial_success[n][, 1:columns_success[n]] .^2) / (2 * square(exp(bandwidth_log[2] + v_compID[compID[n],8] )))
//         )
//       ); 
//       
//         success_time_feature = (
//         exp(
//           - (time_success[n][, 1:columns_success[n]] .^2) / (2 * square(exp(bandwidth_log[3] + v_compID[compID[n],9])))
//         )
//       ); 
//       
//       success_feature = ( success_spatial_feature .*  success_time_feature)* rep_vector(1.0, columns_success[n]); // element wise mulltiplication is a dot
//       
//       };
//       
//       
// 
//         // Loss kernel
//       
//     matrix[N_spots,columns_loss[n]] loss_spatial_feature;
//     matrix[N_spots,columns_loss[n]] loss_time_feature;
// 
//     vector[N_spots] loss_feature;
//    
//      //compute long social memory
//      if (spatial_loss[n][1, 1]<0){
//         loss_feature =rep_vector(0, N_spots);
//      } else {
//        loss_spatial_feature = (
//         exp(
//           - (spatial_loss[n][, 1:columns_loss[n]] .^2) / (2 * square(exp(bandwidth_log[4] + v_compID[compID[n],10] )))
//         )
//       ); 
//       
//         loss_time_feature = (
//         exp(
//           - (time_loss[n][, 1:columns_loss[n]] .^2) / (2 * square(exp(bandwidth_log[5] + v_compID[compID[n],11])))
//         )
//       ); 
//       
//       loss_feature = ( loss_spatial_feature .*  loss_time_feature)* rep_vector(1.0, columns_loss[n]); // element wise mulltiplication is a dot
//       
//       };
// 
//       vector[N_spots] eta;
//       eta =
//          (beta[1] + v_compID[compID[n],1] + v_partID[partID[n],1]  + v_tripID[tripID[n],1]) * success_feature
//         + (beta[2] + v_compID[compID[n],2] + v_partID[partID[n],2] + v_tripID[tripID[n],2]) * loss_feature
//         + (beta[3] + v_compID[compID[n],3] + v_partID[partID[n],3] + v_tripID[tripID[n],3]) * social_feature
//         + (beta[4] + v_compID[compID[n],4] + v_partID[partID[n],4] + v_tripID[tripID[n],4]) * sonar_distance[n, ]'
//         //+ (beta[5] + v_compID[compID[n],5] + v_partID[partID[n],5] + v_tripID[tripID[n],5]) * density3[n, ]'
//         + (beta[5] + v_compID[compID[n],5] + v_partID[partID[n],5] + v_tripID[tripID[n],5]) * locality[n, ]'
//         + (beta[6] + v_compID[compID[n],6] + v_partID[partID[n],6] + v_tripID[tripID[n],6]) * locality[n, ]'*time_since_start[n];
// 
//       log_lik[n] = categorical_logit_lpmf(1 | eta);
//   }}
  
