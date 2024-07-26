// TENSOR PRODUCT FOR REAL DATA

// QUARTERLY SCALE
// AGE-DEPENDENT DIAGNOSES
// NO UNDER-REPORTING

// AGE DEPENDENT DIAGNOSES 
// In usual age-dependent setting (agediag and agediag1, agediag2) 
// diag probs are modelled with rw
// Here we use tensor product splines
// for each of the states... fewer parameters and continuous age modelling

functions{

  vector row_sums(matrix X) {
    vector[rows(X)] s ;  
    s = X * rep_vector(1, cols(X)) ;
    return s ;
  }
  
  row_vector col_sums(matrix X) {
    row_vector[cols(X)] s ;
    s = rep_row_vector(1, rows(X)) * X ;
    return s ;
  }
  
  // Q_fct creates quarterly transition (only, no diag) matrices 
  // Assume diagnosis occurs first and if not then progression
  matrix Q_fct(matrix[] d, matrix q, int t, int a, int a0){
    matrix[4,4] Q;
    Q[1,1] = (1 - d[1,t,a]) * (1 - q[1,a0]) ;
    Q[1,2] = (1 - d[1,t,a]) * q[1,a0] ;
    Q[1,3] = 0 ;
    Q[1,4] = 0 ;
    Q[2,1] = 0 ;
    Q[2,2] = (1 - d[2,t,a]) * (1 - q[2,a0]) ;
    Q[2,3] = (1 - d[2,t,a]) * q[2,a0] ;
    Q[2,4] = 0 ;
    Q[3,1] = 0 ;
    Q[3,2] = 0 ;
    Q[3,3] = (1 - d[3,t,a]) * (1 - q[3,a0]) ;
    Q[3,4] = (1 - d[3,t,a]) * q[3,a0] ;
    Q[4,1] = 0 ;
    Q[4,2] = 0 ;
    Q[4,3] = 0 ;
    Q[4,4] = (1 - d[4,t,a]) * (1 - q[4,a0]) ;
    return Q;
  }
  
  int age_fct(int a0, int t, int t0, int nquar, int nage){
    // Assume individuals age at quarters 5,9,13...
    // So when modulus(t,4) == 1
    real rem;
    int a; // current age
    int addyr;
  
    a = a0; // at t0, current age =a0
    for(t1 in min(t0+1,nquar):t){
      rem = t1 % 4; 
       if(rem==1) addyr = 1;
       else addyr = 0;
      a = a + addyr;
    }
    
    return min(a,nage);
  }
  
  matrix[,] diag_arr_fct(matrix[,] diag_arr, matrix e_m, int nquar, int nage){
    matrix[nquar,5] diag_arr1[nquar,nage];
    row_vector[5] out[nquar,nage];
    row_vector[5] x;
    int a; 
  
    // Setting all entries of diag_arr1 to 0
    // (default: no new arrival)
    // Needed: if not get nan
    for(t0 in 1:nquar){
      for(a0 in 1:nage){
         diag_arr1[t0,a0] = e_m;
      }
    }
    
    // Specific t0,a0,t for which new arrivals
    // are non zero
    for(t0 in 1:nquar){
      for(a0 in 1:nage){
         for(t in t0:nquar){
            a = age_fct(a0, t, t0, nquar, nage);
            x = diag_arr1[t,a,t0,]; // To deal with last cumulative age
            diag_arr1[t,a,t0,] = x + diag_arr[t,a0,t0,]; 
         }
      }
    }
    return diag_arr1;
  }

}


data{
  int<lower=1> nquar; // number of years 
  int<lower=1> nage; // number of ages 
  int<lower=1> ninfpars; // number of infection parameters - here equal to number of diagnosis parameters 
  int<lower=0> HIV[nquar,nage] ; // data with new HIV diag by year and age
  int<lower=0> AIDS[nquar,nage] ; // data with new AIDS diag by year and age
  int<lower=0> CD4[nquar,nage,4] ; // CD4 data by number of years and ages and CD4 states (4)
  matrix<lower=0>[4,nage] init_prev ; // Initial prevalence by age at inf and CD4 states (4) 
  matrix<lower=0,upper=1>[4,nage] q;// matrix of progressions for different ages.. fixed
  matrix[nquar*nage,ninfpars] X; // design matrix of spline, note this is the same for infections and 4 diagnoses states - if we assume same number of parameters
  matrix[ninfpars-1,2*ninfpars-2] S1; // penalty matrix
}

transformed data{
  int ndiagpars = ninfpars;
  vector[3] empty_vec = rep_vector(0,3);
  row_vector[5] empty_vec1 = rep_row_vector(0,5);
  vector[4] empty_vec2 = rep_vector(0,4);
  matrix[nquar,5] empty_mat = rep_matrix(0,nquar,5);
  vector[ninfpars-1] zero_vec = rep_vector(0,ninfpars-1);
}

parameters{
  vector[ninfpars] beta; //beta are the infection parameters
  vector<lower=0>[2] sigma; // smoothing parameters
  vector<lower=0>[2] sigma_d1; // smoothing parameters, diag st 1
  vector<lower=0>[2] sigma_d2; // smoothing parameters, diag st 2
  vector<lower=0>[2] sigma_d3; // smoothing parameters, diag st 3
  vector<lower=0>[2] sigma_d4; // smoothing parameters, diag st 4
  vector[ndiagpars] delta1; 
  vector[ndiagpars] delta2; 
  vector[ndiagpars] delta3; 
  vector[ndiagpars] delta4; 
}

transformed parameters{
  // NOTE: In STAN all variable declarations at beginnig of block
  // before any manipulation occurs
  vector<lower=0>[nquar*nage] H; //vector of expected log responses
  vector<lower=0, upper=1>[nquar*nage] d1; // logit-diagnosis probabilities, state 1
  vector<lower=0, upper=1>[nquar*nage] d2; // logit-diagnosis probabilities, state 2
  vector<lower=0, upper=1>[nquar*nage] d3; // logit-diagnosis probabilities, state 3
  vector<lower=0, upper=1>[nquar*nage] d4; // logit-diagnosis probabilities, state 4
  matrix[nquar,nage] d[4];// diagnoses for different diag state, times and current ages
  matrix<lower=0, upper=1>[4,4] P; // transition matrix
  vector<lower=0>[4] lat_arr[nquar]; // expected number of diag at time t age a state k
  vector<lower=0>[4] lat_p;
  row_vector<lower=0>[5] diag_p;
  matrix[nquar, 5] diag_arr[nquar,nage]; // expected number of diag at time t age a state k
  matrix<lower=0>[nquar, 5] diag_arr1[nquar,nage]; // expected number of diag at time t age a state k
  vector<lower=0>[4] temp; // expected number of diag at time t age a state k
  row_vector<lower=0>[5] exp_diag[nquar,nage]; // expected number of diag at time t age a state k
  matrix[nquar,nage] exp_HIV; // expected number of HIV diag at time t age a
  matrix[nquar,nage] exp_AIDS; // expected number of HIV diag at time t age a
  simplex[4] p_tak[nquar,nage]; // proportion of HIV diagnoses in state k (time t and age a)
  matrix[ninfpars-1,ninfpars-1] K1_inv; // precision matrix of spline prior for infections
  matrix[ndiagpars-1,ndiagpars-1] D1_inv; // precision matrix of spline prior for diagnosis from state 1 
  matrix[ndiagpars-1,ndiagpars-1] D2_inv; // precision matrix of spline prior for diagnosis from state 2
  matrix[ndiagpars-1,ndiagpars-1] D3_inv; // precision matrix of spline prior for diagnosis from state 3
  matrix[ndiagpars-1,ndiagpars-1] D4_inv; // precision matrix of spline prior for diagnosis from state 4
  vector<lower=0>[2] lambda; // smoothing parameters for infection spline
  vector<lower=0>[2] lambda_d1; // smoothing parameters for diagnosis spline, state 1
  vector<lower=0>[2] lambda_d2; // smoothing parameters for diagnosis spline, state 2
  vector<lower=0>[2] lambda_d3; // smoothing parameters for diagnosis spline, state 3
  vector<lower=0>[2] lambda_d4; // smoothing parameters for diagnosis spline, state 4

  lambda = inv(square(sigma));
  lambda_d1 = inv(square(sigma_d1));
  lambda_d2 = inv(square(sigma_d2));
  lambda_d3 = inv(square(sigma_d3));
  lambda_d4 = inv(square(sigma_d4));
  
  // PRIORS for the tensor product
  // incidence
  K1_inv = S1[1:(ninfpars-1),1:(ninfpars-1)] * lambda[1]  + S1[1:(ninfpars-1),ninfpars:(2*ninfpars-2)] * lambda[2];
  
  //diag probs
  D1_inv = S1[1:(ninfpars-1),1:(ninfpars-1)] * lambda_d1[1]  + S1[1:(ninfpars-1),ninfpars:(2*ninfpars-2)] * lambda_d1[2];
  D2_inv = S1[1:(ninfpars-1),1:(ninfpars-1)] * lambda_d2[1]  + S1[1:(ninfpars-1),ninfpars:(2*ninfpars-2)] * lambda_d2[2];
  D3_inv = S1[1:(ninfpars-1),1:(ninfpars-1)] * lambda_d3[1]  + S1[1:(ninfpars-1),ninfpars:(2*ninfpars-2)] * lambda_d3[2];
  D4_inv = S1[1:(ninfpars-1),1:(ninfpars-1)] * lambda_d4[1]  + S1[1:(ninfpars-1),ninfpars:(2*ninfpars-2)] * lambda_d4[2];

  // INCIDENCE
  H = exp(X*beta); // These are expected infections
  
  // DIAGNOSES
  // Using a spline
  // d1 to d4 are diagnosis probabilities on the 0-1 scale
  d1 = inv_logit(X*delta1);
  d2 = inv_logit(X*delta2);
  d3 = inv_logit(X*delta3);
  d4 = inv_logit(X*delta4);
  
 // To make notation easier
 // define array d[nstate, nt, nage]
 // for dianosis probabilities
  for(t in 1:nquar){
    for(a in 1:nage){
      d[1,t,a] = d1[(a-1)*nquar+t];
      d[2,t,a] = d2[(a-1)*nquar+t];
      d[3,t,a] = d3[(a-1)*nquar+t];
      d[4,t,a] = d4[(a-1)*nquar+t];
    }
  }

  // NOW NEED TO CREATE THE TRANSITION MATRICES
  
  // Useful to define a zero-matrix
  // For later calculations
  // Wrt to other models I have removed the PA
  // Save memory
  
  for(t0 in 1:nquar){
    for(a0 in 1:nage){
      // At time = time of inf
      // New infs enter model in state 1
      // a = a0 clearly; t=t0
      
      // At time 1 we need to consider initial prevalence
      // Prevalence is latent in 1994Q4 (t = 0)
      // So consider diagnosis and progressions
      lat_p =empty_vec2;
      diag_p =empty_vec1;
      
      if(t0 == 1){
        P = Q_fct(d, q, 1, a0, a0);
        lat_p = P' * init_prev[,a0];
        diag_p[1] = init_prev[1,a0] * d[1,1,a0] ;
        diag_p[2] = init_prev[2,a0] * d[2,1,a0] ;
        diag_p[3] = init_prev[3,a0] * d[3,1,a0] ;
        diag_p[4] = init_prev[4,a0] * d[4,1,a0] ;
        diag_p[5] = init_prev[4,a0] * (1-d[4,1,a0]) * q[4,a0] ;
      }
      // No diagnosis or progression in first quarter...
      lat_arr[t0,1] = H[(a0-1)*nquar+t0] + lat_p[1];
      lat_arr[t0,2:4] = empty_vec + lat_p[2:4] ;
      diag_arr[t0,a0,t0,] = empty_vec1 + diag_p ;
      
      for(t in (t0+1):nquar){ // might worry if t = nquar, but in STAN loops in (nquar+1):nquar are not evaluated (unlike R, like JAGS) --- saving an if loop
        int a = age_fct(a0, t, t0, nquar, nage);
        P = Q_fct(d, q, t, a, a0);
        
        temp = lat_arr[t-1,];
        lat_arr[t,] = P' * temp;
        // mapping a0,t and t0 to corresponding current age via age_fct
        diag_arr[t,a0,t0,1] = temp[1] * d[1,t,a] ;
        diag_arr[t,a0,t0,2] = temp[2] * d[2,t,a] ;
        diag_arr[t,a0,t0,3] = temp[3] * d[3,t,a] ;
        diag_arr[t,a0,t0,4] = temp[4] * d[4,t,a] ;
        diag_arr[t,a0,t0,5] = temp[4] * (1-d[4,t,a]) * q[4,a0] ;
      }
    }
  }
  
  diag_arr1 = diag_arr_fct(diag_arr, empty_mat, nquar, nage);
  
  for(t in 1:nquar){
    for(a in 1:nage){
      exp_diag[t,a,] = col_sums(diag_arr1[t,a,,]);
      // expected HIV diagnoses
      exp_HIV[t,a] = exp_diag[t,a,1] + exp_diag[t,a,2] + exp_diag[t,a,3] + exp_diag[t,a,4]; 
      // Expected AIDS diagnoses
      exp_AIDS[t,a] = exp_diag[t,a,5]; 
    }
  }
  
  /// For a = 1 and some times (t=5, t=9 ...)
  // there are 0 diagnes arrivals
  // This is due to model assumptions
  // So removing age 1 from likelihood
  // Fixing a=1 to get a valid simplex
  for(t in 1:nquar){
    p_tak[t,1,1] = 1;
    p_tak[t,1,2] = 0;
    p_tak[t,1,3] = 0;
    p_tak[t,1,4] = 0;
  }
  
  for(t in 1:nquar){
    for(a in 2:nage){ 
      // Now we want the proportion of CD4 that
      // is diangosed with HIV in state k at time t and age a (p_kta)
      p_tak[t,a,1] = exp_diag[t,a,1] / exp_HIV[t,a];
      p_tak[t,a,2] = exp_diag[t,a,2] / exp_HIV[t,a];
      p_tak[t,a,3] = exp_diag[t,a,3] / exp_HIV[t,a];
      p_tak[t,a,4] = exp_diag[t,a,4] / exp_HIV[t,a];
    }
  }

}

model{
  // Priors
  // Splines
  // sigma ~ inv_gamma(.05,.005); // s.d of observations. Manual recommends using flat priors, following Gelman's paper
  // If commented default un-normslized uniform prior
  
  sigma ~ student_t(4,0,20); // Weakly informative priors
  sigma_d1 ~ student_t(4,0,20); // Weakly informative priors
  sigma_d2 ~ student_t(4,0,20); // Weakly informative priors
  sigma_d3 ~ student_t(4,0,20); // Weakly informative priors
  sigma_d4 ~ student_t(4,0,20); // Weakly informative priors
  
  beta[1] ~ normal(0,5) ;// intercept (age and time specific infections) on the log scale - exp(10) = 22026.47 - fairly uninformative prior
  beta[2:ninfpars] ~ multi_normal_prec(zero_vec,K1_inv) ;
  
  delta1[1] ~ normal(0,3) ;// intercept from diagnosis probabilities (state 1) on the inverse logit scale - N(0,3) fairly uninformative prior
  delta1[2:ndiagpars] ~ multi_normal_prec(zero_vec,K1_inv) ;

  delta2[1] ~ normal(0,3) ;// intercept from diagnosis probabilities (state 2) on the inverse logit scale - N(0,3) fairly uninformative prior
  delta2[2:ndiagpars] ~ multi_normal_prec(zero_vec,K1_inv) ;

  delta3[1] ~ normal(0,3) ;// intercept from diagnosis probabilities (state 3) on the inverse logit scale - N(0,3) fairly uninformative prior
  delta3[2:ndiagpars] ~ multi_normal_prec(zero_vec,K1_inv) ;

  delta4[1] ~ normal(0,3) ;// intercept from diagnosis probabilities (state 4) on the inverse logit scale - N(0,3) fairly uninformative prior
  delta4[2:ndiagpars] ~ multi_normal_prec(zero_vec,K1_inv) ;

  // Likelihood of model
  for(t in 1:nquar){ // a is in 2 ages, because for age = 1 -> at some t get Po(0)
    // Individuals age at times 5,9,13...
    // So if individuals inf at t=5 at age=1
    // Cannot be dx at t=5 at age age=1 (dx t > 6)
    // remove age = 1 from likelihood
     for(a in 2:nage){ 
      HIV[t,a] ~ poisson(exp_HIV[t,a]); 
      AIDS[t,a] ~ poisson(exp_AIDS[t,a]);
      CD4[t,a,] ~ multinomial(p_tak[t,a,]); 
      }
    }
    
}
  
generated quantities{
// Log-lik calculations to get WAIC / LOO/ DIC
  vector[nquar*nage] log_lik; // has to be vector to be consistent with R pckg
    for(t in 1:nquar){
      log_lik[(1-1)*nquar+t] = poisson_lpmf(HIV[t,1] | exp_HIV[t,1]); 
      for(a in 2:nage){
          log_lik[(a-1)*nquar+t] = poisson_lpmf(HIV[t,a] | exp_HIV[t,a]) + poisson_lpmf(AIDS[t,a] | exp_AIDS[t,a]) +  multinomial_lpmf(CD4[t,a,] | p_tak[t,a,]) ; 
      }
    }
}
