// NOTE FOR EQUIVALENCE WITH JAGS MODEL
// JAGS parametrizes normal in terms of precision = 1 / var
// And in Paul's paper priors were on variances
// STAN parametrizes normal in terms of sd ( = sqrt(var) ) 

functions{
  // Function to create transition matrices - i.e. one large transition matrix
  matrix tmat_fct(int t, vector q, matrix d, vector ones) {
    matrix[13, 13] Lambda;
    
    Lambda = diag_matrix(ones); // Setting diagonal entries to one
    
    Lambda[1,1] = 0;
    Lambda[1,2] = (1 - q[1]) * (1 - d[1,t]);
    Lambda[1,3] = q[1] * (1 - d[1,t]);
    Lambda[1,9] = d[1,t];
    
    Lambda[2,2] = 0;
    Lambda[2,4] = (1 - q[1]) * (1 - d[1,t]);
    Lambda[2,5] = q[1] * (1 - d[1,t]);
    Lambda[2,9] = d[1,t];
    
    Lambda[3,3] = 0;
    Lambda[3,5] = (1 - q[2]) * (1 - d[1,t]);
    Lambda[3,6] = q[2] * (1 - d[1,t]);
    Lambda[3,9] = d[1,t];
    
    Lambda[4,4] = (1 - q[1]) * (1 - d[2,t]);
    Lambda[4,5] = q[1] * (1 - d[2,t]);
    Lambda[4,10] = d[2,t];
    
    Lambda[5,5] = (1 - q[2]) * (1 - d[3,t]);
    Lambda[5,6] = q[2] * (1 - d[3,t]);
    Lambda[5,11] = d[3,t];
    
    Lambda[6,6] = (1 - q[3]) * (1 - d[4,t]);
    Lambda[6,7] = q[3] * (1 - d[4,t]);
    Lambda[6,12] = d[4,t];
    
    Lambda[7,7] = (1 - q[4]) * (1 - d[5,t]);
    Lambda[7,8] = q[4] * (1 - d[5,t]);
    Lambda[7,13] = d[5,t];
    
    return Lambda;
  }
}

data{
  int<lower=1> nquar; // number of quarters
  int<lower=0> HIV[nquar]; //Number of HIV diagnoses
  int<lower=0> AIDS[nquar] ; // AIDS diag data
  int<lower=0> CD4[nquar,5] ; // diagnosis data by states
  vector<lower=0, upper=1>[4] q; // Fixed progression probabilities
  vector<lower=0>[7] init_prev; // Initial prevalence (as for age-indept model)
  int<lower=1> ninfpars;
  matrix[nquar,ninfpars] X;
}

transformed data{
  vector[13] ones = rep_vector(1, 13);
  vector[12] zeroes = rep_vector(0, 12);
}

parameters{
  vector[ninfpars] beta;
  real<lower=0> sigma;
  vector[nquar] delta_raw1;  // random walk diagnoses
  vector[nquar] delta_raw2; 
  vector[nquar] delta_raw3; 
  vector[nquar] delta_raw4; 
  vector[nquar] delta_raw5; 
  vector<lower=0>[5] vardelta; // variances of diagnoses
  real<lower=0, upper=1> under_rep;
}

transformed parameters{
  vector[nquar] gamma; // log-infections (random walk)
  real<lower=0> lambda;
  vector<lower=0>[nquar] h; // number of infs at time t (exp(gamma))
  matrix<lower=0, upper=1>[5, nquar] d; // matrix of diag probs
  matrix<lower=0>[13,nquar] cum_diags ; // expected cum number of individuals in state s at time t 
  
  row_vector<lower=0>[nquar] incidence_state8; // expected incidence in state 8 (AIDS)
  row_vector<lower=0>[nquar] incidence_state9; // expected incidence in state 9 (RITA)
  row_vector<lower=0>[nquar] incidence_state10; // expected incidence in state 10 (CD4>500)
  row_vector<lower=0>[nquar] incidence_state11; // expected incidence in state 11 (CD4 350-500)
  row_vector<lower=0>[nquar] incidence_state12; // expected incidence in state 12 (CD4 200-350)
  row_vector<lower=0>[nquar] incidence_state13; // expected incidence in state 13 (CD4<350)
  
  simplex[5] exp_p[nquar]; // proportion of HIV diagnoses in state k (at time t)
  row_vector<lower=0>[nquar] exp_HIV_dx; // expected number of HIV dx
  row_vector<lower=0>[nquar] exp_AIDS_dx; // expected number of HIV dx
  matrix<lower=0,upper=1>[13,13] PA[nquar]; // collection of model transition matrices
  vector[13] tmpvec; 
  vector[13] tmpvec1; 
  vector[13] init_vec;

  vector[nquar] delta1;  // random walk diagnoses
  vector[nquar] delta2; 
  vector[nquar] delta3; 
  vector[nquar] delta4; 
  vector[nquar] delta5; 
  matrix<lower=0>[7,nquar] prev_mat; // expected prevalence

  // INCIDENCE
  
  // mean of spline log incidence
  gamma = X*beta;
  h = exp(gamma);
  
  lambda = inv(square(sigma));
  
  // Initializing delta
  delta1[1] = -3.2 + 0.2*delta_raw1[1]; // Diag in 1995 - st1
  delta2[1] = -3.2 + 0.2*delta_raw2[1]; // Diag in 1995 - st2
  delta3[1] = -3.2 + 0.2*delta_raw3[1]; // Diag in 1995 - st3
  delta4[1] = -3.0 + 0.2*delta_raw4[1]; // Diag in 1995 - st4
  delta5[1] = -2.5 + 0.3*delta_raw5[1]; // Diag in 1995 - st5
  
  // Non-centered reparametrization
  // Can probably optimize that 
  // using tail and head
  for(t in 2:nquar){
    delta1[t] = delta1[t-1] + sqrt(vardelta[1])*delta_raw1[t];
    delta2[t] = delta2[t-1] + sqrt(vardelta[2])*delta_raw2[t];
    delta3[t] = delta3[t-1] + sqrt(vardelta[3])*delta_raw3[t];
    delta4[t] = delta4[t-1] + sqrt(vardelta[4])*delta_raw4[t];
    delta5[t] = delta5[t-1] + sqrt(vardelta[5])*delta_raw5[t];
  }

  // DIAGNOSES 
  // Post 1984 (as model starts in 1995)
  for(t in 1:nquar){
    d[1,t] = inv_logit(delta1[t]);
    d[2,t] = inv_logit(delta2[t]);
    d[3,t] = inv_logit(delta3[t]);
    d[4,t] = inv_logit(delta4[t]);
    d[5,t] = inv_logit(delta5[t]);
  }
  
  // EPIDEMIC EVOLUTION
  
  // Creating transition matrices
  for(t in 1:nquar){
    PA[t,,] = tmat_fct(t, q, d, ones);
  }
  
  // Initialization at time 1
  init_vec = append_row(h[1],zeroes); // New infections at time 1
  tmpvec1 = append_row(init_prev,zeroes[1:6]); // 0 initial prevalence in dx states... augmenting vec
  cum_diags[,1] = PA[1,,]' * tmpvec1 + init_vec; // Initial prevalence advancing
  
  // Times 2, ..., nquar
  for(t in 2:nquar){
    tmpvec = append_row(h[t],zeroes);
    tmpvec1 = cum_diags[,t-1];
    cum_diags[,t] = PA[t,,]' * tmpvec1 + tmpvec;
  }
  
  // Incidence
  // Difference in cumulative diagnoses
  
  // No incidence at t=1
  incidence_state8[1] = cum_diags[8,1];
  incidence_state9[1] = cum_diags[9,1];
  incidence_state10[1] = cum_diags[10,1];
  incidence_state11[1] = cum_diags[11,1];
  incidence_state12[1] = cum_diags[12,1];
  incidence_state13[1] = cum_diags[13,1];
  
  // t = 2, ..., nquar
  incidence_state8[2:nquar] = tail(cum_diags[8,], nquar - 1) - head(cum_diags[8,], nquar - 1);
  incidence_state9[2:nquar] = tail(cum_diags[9,], nquar - 1) - head(cum_diags[9,], nquar - 1);
  incidence_state10[2:nquar] = tail(cum_diags[10,], nquar - 1) - head(cum_diags[10,], nquar - 1);
  incidence_state11[2:nquar] = tail(cum_diags[11,], nquar - 1) - head(cum_diags[11,], nquar - 1);
  incidence_state12[2:nquar] = tail(cum_diags[12,], nquar - 1) - head(cum_diags[12,], nquar - 1);
  incidence_state13[2:nquar] = tail(cum_diags[13,], nquar - 1) - head(cum_diags[13,], nquar - 1);

  // Prevalence matrix
  prev_mat = cum_diags[1:7,];
  
  // Expected number of HIV diagnoses
  // Sum of incidences in states 9 to 13
  exp_HIV_dx = incidence_state9 + incidence_state10 + incidence_state11 + incidence_state12 + incidence_state13;
   
  // Expected number of AIDS diagnoses
  // Include under-reporting from year 2000
  // Avoiding if-else loop
  // exp_AIDS_dx[1:nquar] = incidence_state8[1:nquar]*under_rep; 
  exp_AIDS_dx[1:nquar] = incidence_state8[1:nquar]; 

  // Including CD4 counts
  // Reliable CD4 data only available from 1991 (quarter=53) onwards
  for(t in 1:nquar){
    exp_p[t,1] = incidence_state9[t] / exp_HIV_dx[t];
    exp_p[t,2] = incidence_state10[t] / exp_HIV_dx[t];
    exp_p[t,3] = incidence_state11[t] / exp_HIV_dx[t];
    exp_p[t,4] = incidence_state12[t] / exp_HIV_dx[t];
    exp_p[t,5] = incidence_state13[t] / exp_HIV_dx[t];
  }
  
}

model{
  
  // PRIORS
  
  // Infection
  // Infection
  beta[1] ~ normal(5, 1.5); // Vaguely informative priors... corresponds to quarter inf in (exp(2),exp(8)) = (8,4390)
  beta[2:ninfpars] ~ normal(0,sigma);
  sigma ~ student_t(4, 0, 20) T[0,]; // More intuitive to put a prior on sd rather than precision (lambda)
  
  // LIKELIHOOD 
  
  // Non centered parametrization for diagnoses
  delta_raw1 ~ normal(0, 1);
  delta_raw2 ~ normal(0, 1);
  delta_raw3 ~ normal(0, 1);
  delta_raw4 ~ normal(0, 1);
  delta_raw5 ~ normal(0, 1);
  vardelta ~ gamma(1,32); // prior for logit random walk
  
  // Under reporting
  under_rep ~ beta(236, 118);
  
  // Poisson likelihood HIV data
  HIV ~ poisson(exp_HIV_dx);
  
  // Poisson likelihood AIDS data
  AIDS ~ poisson(exp_AIDS_dx);
  
  // Multinomial term CD4
  for(t in 1:nquar){
    CD4[t,] ~ multinomial(exp_p[t,]);
  }
}


generated quantities{
// Log-lik calculations to get WAIC / LOO/ DIC
  vector[nquar] log_lik; 
  for(t in 1:nquar){
    log_lik[t] = poisson_lpmf(HIV[t] | exp_HIV_dx[t]) + poisson_lpmf(AIDS[t] | exp_AIDS_dx[t]) +  multinomial_lpmf(CD4[t,] | exp_p[t,]) ; 
  }
}
  
