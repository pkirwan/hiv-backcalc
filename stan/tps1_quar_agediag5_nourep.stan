// THIN-PLATE REGRESSION SPLINE FOR REAL DATA   
  


// QUARTERLY SCALE  
// AGE-DEPENDENT DIAGNOSES  
// UNDER-REPORTING  
  
  


// AGE DEPENDENT DIAGNOSES   
// Age and state specific intercept that affects all the diag rw  
  
  


// HERE: Have age and state specific intercept  
// ie. the contribution of age can change within state  
  
  


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
matrix Q_fct(matrix d, matrix q, int t, int a){  
matrix[4,4] Q;  
Q[1,1] = (1 - d[1,t]) * (1 - q[1,a]) ;  
Q[1,2] = (1 - d[1,t]) * q[1,a] ;  
Q[1,3] = 0 ;  
Q[1,4] = 0 ;  
Q[2,1] = 0 ;  
Q[2,2] = (1 - d[2,t]) * (1 - q[2,a]) ;  
Q[2,3] = (1 - d[2,t]) * q[2,a] ;  
Q[2,4] = 0 ;  
Q[3,1] = 0 ;  
Q[3,2] = 0 ;  
Q[3,3] = (1 - d[3,t]) * (1 - q[3,a]) ;  
Q[3,4] = (1 - d[3,t]) * q[3,a] ;  
Q[4,1] = 0 ;  
Q[4,2] = 0 ;  
Q[4,3] = 0 ;  
Q[4,4] = (1 - d[4,t]) * (1 - q[4,a]) ;  
return Q;  
}  
  
  


int age_fct(int a0, int t, int t0, int nquar, int nage){  
// Assume individuals age at quarters 5,9,13…  
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
  
  


// Links age to age group  
// age 1-10 -&gt; ag1  
// age 11-10 -&gt; ag2  
// age 21-30 -&gt; ag3  
// age 31-52 -&gt; ag4  
int age_to_agroup(int a){  
int ag;  
if(a &lt;= 10){  
ag = 1;  
}  
else if(a &gt; 10 &amp;&amp; a &lt;= 20){  
ag = 2;  
}  
else if(a &gt; 20 &amp;&amp; a &lt;= 30){  
ag = 3;  
}  
else {  
ag = 4;  
}  
return ag;  
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
int<lower=1> ninfpars; // number of infection parameters   
int<lower=0> nmissing; // number of quarters for which we want the posterior predictive distribution.  
int<lower=0, upper="1"> sigma_inform; // should we use an informative prior on the smoothing parameter of the incidence spline  
int<lower=0> HIV[nquar,nage] ; // data with new HIV diag by year and age  
int<lower=0> AIDS[nquar,nage] ; // data with new AIDS diag by year and age  
int<lower=0> CD4[nquar,nage,4] ; // CD4 data by number of years and ages and CD4 states(3)  
matrix<lower=0>[4,nage] init_prev ; // Initial prevalence by age at inf and CD4 states(3)   
matrix<lower=0,upper=1>[4,nage] q;// matrix of progressions for different ages.. fixed  
matrix[(nquar+nmissing)*nage,ninfpars] X; // design matrix … 50 is number of inf params  
}  
  
  


transformed data{  
vector[3] empty_vec = rep_vector(0,3);  
row_vector[5] empty_vec1 = rep_row_vector(0,5);  
vector[4] empty_vec2 = rep_vector(0,4);  
matrix[nquar+nmissing,5] empty_mat = rep_matrix(0,nquar+nmissing,5);  
vector[ninfpars-1] zero_vec = rep_vector(0,ninfpars-1);  
int nagroup = 4; // number of age groups  
}  
  
  


parameters{  
vector[ninfpars] beta; //beta are the infection parameters  
vector<lower=0>[2] sigma; // smoothing parameters  
vector<lower=0>[4] vardelta;  
vector[nquar+nmissing] delta_raw1;   
vector[nquar+nmissing] delta_raw2;   
vector[nquar+nmissing] delta_raw3;   
vector[nquar+nmissing] delta_raw4;   
vector[nagroup*4] alpha; // state (k) and age (ac) specific intercept for dx probs —   
}  
  
  


transformed parameters{  
// NOTE: In STAN all variable declarations at beginnig of block  
// before any manipulation occurs  
vector<lower=0>[(nquar+nmissing)*nage] H; //vector of expected log responses  
matrix<lower=0, upper="1">[4,4] PA[nquar+nmissing, nage, nagroup];  
vector<lower=0>[4] lat_arr[nquar+nmissing]; // expected number of diag at time t age a state k  
vector<lower=0>[4] lat_p;  
row_vector<lower=0>[5] diag_p;  
matrix[nquar+nmissing, 5] diag_arr[nquar+nmissing,nage]; // expected number of diag at time t age a state k  
matrix<lower=0>[nquar+nmissing, 5] diag_arr1[nquar+nmissing,nage]; // expected number of diag at time t age a state k  
vector<lower=0>[4] temp; // expected number of diag at time t age a state k  
row_vector<lower=0>[5] exp_diag[nquar+nmissing,nage]; // expected number of diag at time t age a state k  
matrix[nquar+nmissing,nage] exp_HIV; // expected number of HIV diag at time t age a  
matrix[nquar+nmissing,nage] exp_AIDS; // expected number of HIV diag at time t age a  
simplex[4] p_tak[nquar+nmissing,nage]; // proportion of HIV diagnoses in state k (time t and age a)  
matrix[4,nquar+nmissing] d[4];// matrix of diagnoses for different ages and age-groups  
matrix[ninfpars-1,ninfpars-1] K1_inv;   
matrix[nquar+nmissing,nagroup] delta1;   
matrix[nquar+nmissing,nagroup] delta2;   
matrix[nquar+nmissing,nagroup] delta3;   
matrix[nquar+nmissing,nagroup] delta4;   
vector<lower=0>[2] lambda; // smoothing parameters  
  
  


lambda = inv(square(sigma));  
  
  


// INCIDENCE  
H = exp(X*beta); // These are log infection…  
  
  


// DIAGNOSES  
  
  


// Initializing delta  
// including the age specific intercept  
for(ag in 1:nagroup){  
delta1[1,ag] = -3.2 + 0.2*delta_raw1[1] + alpha[(ag-1)*nagroup+1]; // worth reconsidering? Originally for diag in 85  
delta2[1,ag] = -3.2 + 0.2*delta_raw2[1] + alpha[(ag-1)*nagroup+2]; // worth reconsidering? Originally for diag in 85  
delta3[1,ag] = -3.0 + 0.2*delta_raw3[1] + alpha[(ag-1)*nagroup+3]; // worth reconsidering? Originally for diag in 85  
delta4[1,ag] = -2.5 + 0.3*delta_raw4[1] + alpha[(ag-1)*nagroup+4]; // worth reconsidering? Originally for diag in 85  
}  
  
  


// Non-centered reparametrization  
// Can probably optimize that   
// using tail and head  
for(ag in 1:nagroup){  
for(t in 2:(nquar+nmissing)){  
delta1[t,ag] = delta1[t-1,ag] + vardelta[1]*delta_raw1[t];  
delta2[t,ag] = delta2[t-1,ag] + vardelta[2]*delta_raw2[t];  
delta3[t,ag] = delta3[t-1,ag] + vardelta[3]*delta_raw3[t];  
delta4[t,ag] = delta4[t-1,ag] + vardelta[4]*delta_raw4[t];  
}  
}  
  
  


// diag probs are inv logit of delta  
// so that are normalised btwn 0 and 1  
// inv logit not vectorized -&gt; use loop  
for(ag in 1:nagroup){  
// Adding the age-specific intercept  
for(i in 1:(nquar+nmissing)){  
d[ag,1,i] = inv_logit(delta1[i,ag]);  
d[ag,2,i] = inv_logit(delta2[i,ag]);  
d[ag,3,i] = inv_logit(delta3[i,ag]);  
d[ag,4,i] = inv_logit(delta4[i,ag]);  
}  
}  
  
  


// NOW NEED TO CREATE THE TRANSITION MATRICES  
  
  


// Useful to define a zero-matrix  
// For later calculations  
  
  


// Looping over t and a to get trans mat  
// for all ages and times  
for (t in 1:(nquar+nmissing)) {  
for (a in 1:nage) {  
for(ag in 1:nagroup){  
//QUARTERLY progression probabilities to observed states  
PA[t,a,ag,,] = Q_fct(d[ag,,], q, t, a);  
}  
}  
}  
  
  


for(t0 in 1:(nquar+nmissing)){  
for(a0 in 1:nage){  
// At time = time of inf  
// New infs enter model in state 1  

    
    
      // At time 1 we need to consider initial prevalence
      // Prevalence is latent in 1994Q4 (t = 0)
      // So consider diagnosis and progressions
      lat_p =empty_vec2;
      diag_p =empty_vec1;
            
      if(t0 == 1){
        int ag0 = age_to_agroup(a0);
        lat_p = (PA[1,a0,ag0,]' * init_prev[,a0]);
        diag_p[1] = init_prev[1,a0] * d[ag0,1,1] ;
        diag_p[2] = init_prev[2,a0] * d[ag0,2,1] ;
        diag_p[3] = init_prev[3,a0] * d[ag0,3,1] ;
        diag_p[4] = init_prev[4,a0] * d[ag0,4,1] ;
        diag_p[5] = init_prev[4,a0] * (1-d[ag0,4,1]) * q[4,a0] ;
      }
      lat_arr[t0,1] = H[(a0-1)*(nquar+nmissing)+t0] + lat_p[1];
      lat_arr[t0,2:4] = empty_vec + lat_p[2:4] ;
      diag_arr[t0,a0,t0,] = empty_vec1 + diag_p ;
    
      for(t in (t0+1):(nquar+nmissing)){ // might worry if t = nquar+nmissing, but in STAN loops in (nquar+1):nquar are not evaluated (unlike R, like JAGS) --- saving an if loop
        int a = age_fct(a0, t, t0, nquar+nmissing, nage);
        int ag = age_to_agroup(a);
        
        temp = lat_arr[t-1,] ;
        lat_arr[t,] = (PA[t,a0,ag,]' * temp);
        // mapping a0,t and t0 to corresponding current age via age_fct
        diag_arr[t,a0,t0,1] = temp[1] * d[ag,1,t] ;
        diag_arr[t,a0,t0,2] = temp[2] * d[ag,2,t] ;
        diag_arr[t,a0,t0,3] = temp[3] * d[ag,3,t] ;
        diag_arr[t,a0,t0,4] = temp[4] * d[ag,4,t] ;
        diag_arr[t,a0,t0,5] = temp[4] * (1-d[ag,4,t]) * q[4,a0] ;
      }
    }

}  
  
  


diag_arr1 = diag_arr_fct(diag_arr, empty_mat, nquar+nmissing, nage);  
  
  


for(t in 1:(nquar+nmissing)){  
for(a in 1:nage){  
exp_diag[t,a,] = col_sums(diag_arr1[t,a,,]);  
// expected HIV diagnoses  
exp_HIV[t,a] = exp_diag[t,a,1] + exp_diag[t,a,2] + exp_diag[t,a,3] + exp_diag[t,a,4];   
// Expected AIDS diagnoses  
exp_AIDS[t,a] = exp_diag[t,a,5];   
}  
}  
  
  


/// For a = 1 and some times (t=5, t=9 …)  
// there are 0 diagnes arrivals  
// This is due to model assumptions  
// So removing age 1 from likelihood  
// Fixing a=1 to get a valid simplex  
for(t in 1:(nquar+nmissing)){  
p_tak[t,1,1] = 1;  
p_tak[t,1,2] = 0;  
p_tak[t,1,3] = 0;  
p_tak[t,1,4] = 0;  
}  
  
  


for(t in 1:(nquar+nmissing)){  
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
// sigma ~ inv_gamma(.05,.005); // s.d of observations. Manual recommends using flat priors, following Gelman’s paper  
// If commented default un-normslized uniform prior  
  
  
if(sigma_inform)  
sigma ~ inv_gamma(1.125,1.125); // More informative prior  
else  
sigma ~ student_t(4,0,200); // Weakly informative priors  
  
  


beta[1] ~ normal(0,30) ;// 30 s.d corresponding to precision given in JAGS  
beta[2:(ninfpars-2)] ~ normal(0,sigma[1]) ;  
beta[(ninfpars-1):ninfpars] ~ normal(0,sigma[2]) ;  
  
  


vardelta ~ gamma(1, 8); // prior for logit random walk  
  
  


alpha ~ normal(0,1); // Priors for age specific intercepts  
  
  


// Non centered parametrization for diagnoses  
delta_raw1 ~ normal(0, 1);  
delta_raw2 ~ normal(0, 1);  
delta_raw3 ~ normal(0, 1);  
delta_raw4 ~ normal(0, 1);  
  
  


// Likelihood of model  
for(t in 1:nquar){ // a is in 2 ages, because for age = 1 -&gt; at some t get Po(0)  
// Individuals age at times 5,9,13…  
// So if individuals inf at t=5 at age=1  
// Cannot be dx at t=5 at age age=1 (dx t &gt; 6)  
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
int<lower=0> HIV_sim[nmissing,nage]; // simulated missing HIV diagnosis data. Should this be declared as a matrix somehow?  
int<lower=0> AIDS_sim[nmissing,nage]; // simulated missing AIDS diagnosis data.  
for(t in 1:nquar){  
log_lik[(1-1)*nquar+t] = poisson_lpmf(HIV[t,1] | exp_HIV[t,1]);   
for(a in 2:nage){  
log_lik[(a-1)*nquar+t] = poisson_lpmf(HIV[t,a] | exp_HIV[t,a]) + poisson_lpmf(AIDS[t,a] | exp_AIDS[t,a]) + multinomial_lpmf(CD4[t,a,] | p_tak[t,a,]) ;   
}  
}  
if(nmissing) {   
for(t in 1:nmissing){   
for(a in 1:nage){   
HIV_sim[t,a] = poisson_rng(exp_HIV[nquar+t,nage]);   
AIDS_sim[t,a] = poisson_rng(exp_AIDS[nquar+t,nage]);   
}  
}  
}  
}  

