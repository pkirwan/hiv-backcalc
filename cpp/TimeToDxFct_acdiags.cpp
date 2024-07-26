// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using std::vector;
using namespace arma;


// Useful function to construct a continuous sequence
// taken from 
// https://github.com/coatless/r-to-armadillo/blob/master/src/seq.cpp
arma::vec seq_int(int a, int b){
  int d = std::abs(b-a)+1;
  return arma::linspace(a, b, d);
}


// [[Rcpp::export]]
int age_class_fct(int a){
  int ac;
  if(a <= 10){ ac=1;}
  else if(a > 10 && a <= 20){ ac=2;}
  else if(a > 20 && a <= 30){ ac=3;}
  else ac = 4;
  return ac;
} 


// probs_move_fct_t_to_d calculates the probability of moving
// out of the current state, for an infection at (t0,a0), at time t (>=t0)
// The idea is to asses when the (t0,a0) infection will be diagnosed
// Hence simulating from t0 to the end of epidemic (nt)
// prog depends on a0, whereas diags from a and t (current time)

// [[Rcpp::export]]
vec probs_move_fct_t_to_d(int t0, int a0, int st, int t, int nt, int nage, const cube& d, const mat& q){
  if(t < t0) {throw std::range_error("t cannot be smaller than t0");}
  if(st==1 & t!=t0) {throw std::range_error("In the first state t must be equal to t0");}
  int time_to_end = nt - t;
  vec probs(time_to_end);
  vec probs0(time_to_end);
  probs.fill(0);
  probs0.fill(0);
  for(int i=0; i<time_to_end; ++i){
    t = t + 1; // + 1 for consistency with R
    int aa = a0 + floor((t-1)/4);
    int a = std::min(aa,nage);
    int ac = age_class_fct(a);
    // Defining probabilities of moving (either dx or prog)
    // In state k, and at time t (t0 + i)
    probs0(i) = d(ac-1,st-1,t-1) + (1-d(ac-1,st-1,t-1))*q(st-1,a0-1); // all the -1 are for consistency with R
  }
  const vec cp = cumprod(1-probs0);
  vec one_vec=vec(1);
  one_vec(0)=1;
  vec cp1 = one_vec;
  if(time_to_end > 1) cp1 = join_cols(one_vec, cp.subvec(0,time_to_end-2)); // -2 because of R/C++ diff
  probs = probs0 % cp1;
  return(probs);
}

// probs_move_fct_snap calculates the snapshot probabilities
// The idea is to asses when the diagnosis levels at (t0,a0) 
// assuming that they will carry on infinetly
// Infinity is considered by choosing a particularly large censoring time 
// Hence simulating from t0 to the end of epidemic (nt)
// prog depends on a0, whereas diags from a and t0 (time at infection)

// [[Rcpp::export]]
vec probs_move_fct_snap(int t0, int a0, int st, int t, int nt, int nage, int t_max, const cube& d, const mat& q){
  if(t < t0) {throw std::range_error("t cannot be smaller than t0");}
  if(st==1 & t!=t0) {throw std::range_error("In the first state t must be equal to t0");}
  vec probs(t_max);
  vec probs0(t_max);
  probs.fill(0);
  probs0.fill(0);
  for(int i=0; i<t_max; ++i){
    t = t + 1; // + 1 for consistency with R
    int aa = a0 + floor((t-1)/4);
    int a = std::min(aa,nage);
    int ac = age_class_fct(a);
    // Defining probabilities of moving (either dx or prog)
    // In state k, and at time t (t0 + i)
    probs0(i) = d(ac-1,st-1,t0-1) + (1-d(ac-1,st-1,t0-1))*q(st-1,a0-1); // all the -1 are for consistency with R
  }
  const vec cp = cumprod(1-probs0);
  vec one_vec=vec(1);
  one_vec(0)=1;
  vec cp1 = join_cols(one_vec, cp.subvec(0,t_max-2)); // -2 because of R/C++ diff
  probs = probs0 % cp1;
  return(probs);
}

// probs_diag_fct_t_to_d function to calculate the probability
// of whether the move event is a diagnosis or progression
// assuming t_to_d dynamics
// i.e. probs depend on current time and age, and progression on age at inf

// [[Rcpp::export]]
double probs_diag_fct_t_to_d(int t0, int a0, int st, int t, int nt, int nage, const cube& d, const mat& q){
  int aa = a0 + floor((t-1)/4);
  int a = std::min(aa,nage);
  int ac = age_class_fct(a);
  double dd = d(ac-1,st-1,t-1);
  double d_star = dd + (1-dd)*q(st-1,a0-1); // all the -1 are for consistency with R
  double out = dd / d_star;
  return(out);
}

// probs_diag_fct_t_to_d function to calculate the probability
// of whether the move event is a diagnosis or progression
// assuming t_to_d dynamics
// i.e. probs depend on current age and time at dx, and progression on age at inf

// [[Rcpp::export]]
double probs_diag_fct_snap(int t0, int a0, int st, int t, int nt, int nage, const cube& d, const mat& q){
  int aa = a0 + floor((t-1)/4);
  int a = std::min(aa,nage);
  int ac = age_class_fct(a);
  double dd = d(ac-1,st-1,t0-1);
  double d_star = dd + (1-dd)*q(st-1,a0-1); // all the -1 are for consistency with R
  double out = dd / d_star;
  return(out);
}

// one_sim_t_to_d simulates the epidemic evolution
// of one infection at (t0,a0)
// interested in time and state of diagnosis

// [[Rcpp::export]]
vec one_sim_t_to_d(int t0, int a0, int nt, int nage, const cube& d, const mat& q){
  // Initializing quantities
  int dx_time = NA_INTEGER;
  int st_dx = 0; // NA would be better but creates issues later ... st_dx=0 means no diag
  int k = 1;
  int t = t0; // Time to diagnsis - initialized to t0
  int dx_ind = 0; // indicator of whether diagnosis has occurred
  vec out(2);
  vec prob_d(2);
  vec sum_to_one_vec=vec(1);
  vec t_m(1);
  vec dx(1);

  // probability vector
  while(dx_ind == 0){
//    Rcpp::Rcout << "t is " << std::endl << t << std::endl;
    vec prob_m = probs_move_fct_t_to_d(t0,a0,k,t,nt,nage,d,q);
    int t_max = prob_m.size();
    // Adding the "censoring probability"
    vec to_sam = seq_int(1,t_max+1);
    sum_to_one_vec(0)=1-sum(prob_m);
    // To avoid loss of precision
    if(sum_to_one_vec(0) < 1e-9) sum_to_one_vec(0) = 0;
    vec prob = join_cols(prob_m,sum_to_one_vec);
    // sampling the waiting time
    t_m = Rcpp::RcppArmadillo::sample(to_sam, 1, true, prob); // note sample does not work without Rcpp::RcppArmadillo
//    Rcpp::Rcout << "t_m is " << std::endl << t_m(0) << std::endl;
    // Below censoring event
    if(t_m(0)==(t_max+1) ){
      dx_time = nt + 1; // NA ideal but causes problems in analysis later
      break;
    }
    // If not censored update the time
    t = t + t_m(0);
//    Rcpp::Rcout << "t (updated) is " << std::endl << t << std::endl;
    // Now want to know whether event is a diagnosis
    double pr_d = probs_diag_fct_t_to_d(t0,a0,k,t,nt,nage,d,q);
    prob_d(0) = 1-pr_d;
    prob_d(1) = pr_d;
    // Rcout << "The value of prob_d is " << prob_d << std::endl;
    dx = Rcpp::RcppArmadillo::sample(seq_int(0,1), 1, true, prob_d); // 0 progression, 1 diagnosis
    dx_ind = dx(0); // need to convert sampling to integer otherwise does not work
//    Rcpp::Rcout << "dx_ind is " << std::endl << dx_ind << std::endl;
    // Rcout << "The value of dx_ind is " << dx_ind << std::endl;
    // Store diagnosis time and state if have diagnosis
    if(dx_ind == 1){
      dx_time = t - t0;
      st_dx = k;
    }
    // this is the special case that there is a progression event in last quarter
    // CENSORING need to break loop and do not count this as dx
    if(dx_ind == 0 && t==nt){
      dx_time = nt + 1; // NA ideal but causes problems in analysis later
      break;
    }
    k = k + 1;
    if(k == 5){ // If reaches state 5 - AIDS
      dx_time = t - t0;
      st_dx = 5;
      break;
    }
  }

  out(0) = dx_time;
  out(1) = st_dx;

  return(out);
}

// [[Rcpp::export]]
vec one_sim_snap(int t0, int a0, int nt, int nage, int t_max, const cube& d, const mat& q){
  // Initializing quantities
  int dx_time = NA_INTEGER;
  int st_dx = 0; // NA would be better but creates issues later ... st_dx=0 means no diag
  int k = 1;
  int t = t0; // Time to diagnsis - initialized to t0
  int dx_ind = 0; // indicator of whether diagnosis has occurred
  vec out(2);
  
  // For snapshot case we know for how long the epidemic will be simulated
  // know a priori - i.e. t_m (declare here and not in loop, efficiency)
  vec prob_m(t_max);
  vec to_sam = seq_int(1,t_max+1);
  vec t_m;
  vec sum_to_one_vec(1);
  vec prob_d(2);
  vec prob(t_max+1);
  vec dx(1);
  // probability vector
  while(dx_ind == 0){
//    Rcpp::Rcout << "t is " << std::endl << t << std::endl;
    prob_m = probs_move_fct_snap(t0,a0,k,t,nt,nage,t_max,d,q);
    // Adding the "censoring probability" -- this should never happen
    // If it does t_max need to be increased...
    sum_to_one_vec(0)=1-accu(prob_m);
    // To avoid loss of precision
    if(sum_to_one_vec(0) < 1e-9) sum_to_one_vec(0) = 0;
    prob = join_cols(prob_m,sum_to_one_vec);
//    Rcout << "The value of sum_to_one_vec is " << sum_to_one_vec << std::endl;
    // sampling the waiting time
    t_m = Rcpp::RcppArmadillo::sample(to_sam, 1, true, prob); // note sample does not work without Rcpp::RcppArmadillo
//    Rcpp::Rcout << "t_m is " << std::endl << t_m(0) << std::endl;
    // Below censoring event
    if(t_m(0)==(t_max+1) ){
      dx_time = 999; // NA ideal but causes problems in analysis later
      break;
    }
    // If not censored update the time
    t = t + t_m(0);
//    Rcpp::Rcout << "t (updated) is " << std::endl << t << std::endl;
    // Now want to know whether event is a diagnosis
    double pr_d = probs_diag_fct_snap(t0,a0,k,t,nt,nage,d,q);
    prob_d(0) = 1-pr_d;
    prob_d(1) = pr_d;
//    Rcout << "The value of prob_d is " << prob_d << std::endl;
    dx = Rcpp::RcppArmadillo::sample(seq_int(0,1), 1, true, prob_d); // 0 progression, 1 diagnosis
    dx_ind = dx(0); // need to convert sampling to integer otherwise does not work
//    Rcpp::Rcout << "dx_ind is " << std::endl << dx_ind << std::endl;
    // Store diagnosis time and state if have diagnosis
    if(dx_ind == 1){
      dx_time = t - t0;
      st_dx = k;
    }
    // this is the special case that there is a progression event in last quarter
    // CENSORING need to break loop and do not count this as dx
    if(dx_ind == 0 && t_m(0)==(t_max+1)){
      dx_time = 999; // NA ideal but causes problems in analysis later
      break;
    }
    k = k + 1;
//    Rcpp::Rcout << "k (updated) is " << std::endl << dx_ind << std::endl;
    if(k == 5){ // If reaches state 5 - AIDS
      dx_time = t - t0;
      st_dx = 5;
      break;
    }
  }

  out(0) = dx_time;
  out(1) = st_dx;

  return(out);
}



// [[Rcpp::export]]
List t_to_d(int nt, int nage, int nst, int n_ind_sim, const mat& d_mat, const mat& q){
  int niter = d_mat.n_rows;
  int nac=4;
  cube t_dx(nt-1,nage,n_ind_sim*niter); // nt-1, because if someone infected at nt... by assumption can not be dx
  cube st_dx(nt-1,nage,n_ind_sim*niter);
  for(int it=0; it < niter; ++it){
    vec d1 = d_mat.row(it).t();
    // converting the diagnosis probabilties as a cube (3D array) d[ac,st,t]
    arma::cube d(d1.begin(), nac, nst, nt);
    for(int t0=0; t0<(nt-1); ++t0){
      for(int a0=0; a0<nage; ++a0){
        for(int ind=0; ind<n_ind_sim; ++ind){
          vec sim = one_sim_t_to_d(t0+1, a0+1, nt, nage, d, q);
          t_dx(t0,a0,n_ind_sim*it+ind) = sim(0);
          st_dx(t0,a0,n_ind_sim*it+ind) = sim(1);
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("time_dx") = t_dx,
                            Rcpp::Named("st_dx") = st_dx);
}

// [[Rcpp::export]]
List snap(int nt, int nage, int nst, int n_ind_sim, int t_max,  const mat& d_mat, const mat& q){
  int niter = d_mat.n_rows;
  int nac=4;
  // Useful constants... could and should be passed as arguments
  cube t_dx(nt,nage,n_ind_sim*niter); // nt-1, because if someone infected at nt... by assumption can not be dx
  cube st_dx(nt,nage,n_ind_sim*niter);
  for(int it=0; it < niter; ++it){
    vec d1 = d_mat.row(it).t();
    // converting the diagnosis probabilties as a cube (3D array) d[ac,st,t]
    arma::cube d(d1.begin(), nac, nst, nt);
    for(int t0=0; t0<(nt); ++t0){
      for(int a0=0; a0<nage; ++a0){
        for(int ind=0; ind<n_ind_sim; ++ind){
          vec sim = one_sim_snap(t0+1, a0+1, nt, nage, t_max, d, q);
          t_dx(t0,a0,n_ind_sim*it+ind) = sim(0);
          st_dx(t0,a0,n_ind_sim*it+ind) = sim(1);
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("time_dx") = t_dx,
                            Rcpp::Named("st_dx") = st_dx);
}
