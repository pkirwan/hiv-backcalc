// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// CCS
using std::vector;


// [[Rcpp::export]]
mat normalize_mat(arma::mat m){
  mat out=m;
  vec row_sums = sum(m, 1);
  for(unsigned int i=0; i<m.n_rows; ++i){
    out.row(i) = (m.row(i) / row_sums[i]);
  }
  return out;
}

// [[Rcpp::export]]
mat normalize_mat_1(arma::mat m){
  mat out=m;
  double s = accu(m);
  for(unsigned int i=0; i<m.n_rows; ++i){
    out.row(i) = (m.row(i) / s);
  }
  return out;
}

// [[Rcpp::export]]
cube trans_mat(const arma::vec& d1, const arma::vec& d2, const arma::vec& d3, const arma::vec& d4, const arma::vec& q){
  int t_max = d1.n_elem; // epidemic times
  int ncd4_st = q.n_elem;
  int nst = 2*ncd4_st+1;

  cube trans_mats(nst, nst, t_max);
  cube prevalence(t_max, nst, t_max);

  mat empty_mat(nst, nst, fill::eye);
  mat empty_mat1(t_max, nst, fill::zeros);

  for(int t=0; t<t_max;++t){
    trans_mats.slice(t) = empty_mat; // Initializing...

    // Creating vector of dx probabilities at time t
    vec d(ncd4_st);
    d(0) = d1(t); d(1) = d2(t); d(2) = d3(t); d(3) = d4(t);

    for(int j=0; j<ncd4_st;++j){
      trans_mats(j,j,t) = (1-q(j))*(1-d(j));
      trans_mats(j,j+1,t) = q(j)*(1-d(j));
      trans_mats(j, j + ncd4_st + 1,t) = d(j);
    }

  }
  return trans_mats;
}

// [[Rcpp::export]]
mat prev_t(int t, const arma::vec& CD4_prev, const arma::vec& h, const cube& trans_mats){
  int t_max = h.n_elem;
  int ncd4_st = CD4_prev.n_elem;
  int nst = 2*ncd4_st+1;

  vec empty_vec(ncd4_st+1, fill::zeros);
  vec prev0 = join_cols(CD4_prev, empty_vec);
  mat prev(t_max,nst,fill::zeros);

  vec infs(nst, fill::zeros);
  infs(0) = h(t);


  prev.row(t) = infs.t(); // Can change to row/col vector to avoid .t()
  if(t==0){ // In case of initial prevalence
    prev.row(t) = (infs + prev0).t();
  }

  for(int i=t+1; i<t_max;++i){
    prev.row(i) = prev.row(i-1) * trans_mats.slice(i);
  }
  return prev;
}

// [[Rcpp::export]]
List forward_calculate(const arma::vec& CD4_prev, const arma::vec& h, const arma::vec& d1, const arma::vec& d2, const arma::vec& d3, const arma::vec& d4, const arma::vec& q){
  int t_max = h.n_elem;
  int ncd4_st = CD4_prev.n_elem;
  int nst = 2*ncd4_st+1;

  cube out(t_max, nst, t_max, fill::zeros);
  cube tmats = trans_mat(d1, d2, d3, d4, q);
  rowvec empty_vec(nst-ncd4_st, fill::zeros);
  //rowvec diags0t = diags0.t();

  mat HIV_diags(t_max, t_max);
  mat AIDS_diags(t_max, t_max);
  cube exp_p(nst-ncd4_st-1,t_max,t_max, fill::zeros);

  for(int t=0; t<t_max;++t){
     out.slice(t) = prev_t(t, CD4_prev, h, tmats);
    // Obtaining diagnoses
    mat X = out.slice(t);
    mat X1 = X(span::all, span(ncd4_st, nst-1));
    mat diags = join_cols(empty_vec,diff(X1));
    mat HIV_d = diags.cols(1, ncd4_st);
    HIV_diags.row(t) = sum(HIV_d, 1).t();
    AIDS_diags.row(t) = sum(diags.col(0), 1).t();
    exp_p.slice(t) = normalize_mat(HIV_d).t();
  }

  // Now I need to return all quantities of interest
  return List::create(
    _["prevalence"]  = out,
    _["HIV.diagnoses"]  = HIV_diags,
    _["AIDS.diagnoses"] = AIDS_diags,
    _["expected.p"] = exp_p
  ) ;
}

// [[Rcpp::export]]
cube temp_diags_fun(List diagnoses){
  mat diag = diagnoses(1);
  cube prev = diagnoses(3);
  int nst = prev.n_rows;
  int nquar = diag.n_rows;
  cube out(nquar,nst,nquar, fill::zeros);

  for(int i=0; i<nquar; ++i){
    for(int j=0; j<nquar; ++j){
      for(int k=0; k<nst; ++k){
       out(j,k,i) = diag(i,j)*prev(k,j,i);
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
mat diag_to_infec_agg(int diag_time,  const cube& temp_HIV){
  int t_max = temp_HIV.n_rows;
  mat diag_mat(temp_HIV.n_cols,diag_time-1);

  for(int i=0; i<(diag_time-1); ++i){
    mat tmp = temp_HIV.slice(i);
    diag_mat.col(i) = tmp.row(diag_time-1).t(); // -1 R/C++ different indices
  }

  // Adding 0 from diag_time onwards
  // So that all objects have the same size
  mat zero_mat(temp_HIV.n_cols,t_max-diag_time+1,fill::zeros);

  // cbind matrices
  mat out = join_rows(diag_mat, zero_mat);

  return out;
}

// [[Rcpp::export]]
cube diag_to_infec_wrap(const cube& temp_HIV){
  int t_max = temp_HIV.n_rows;
  cube out(temp_HIV.n_cols,t_max,t_max,fill::zeros);

  for(int i=1; i<t_max; ++i){
    out.slice(i) = normalize_mat_1(diag_to_infec_agg(i+1,  temp_HIV));
  }
  return out;
}


// [[Rcpp::export]]
vec diag_to_infec_agg1(int diag_time,  const cube& prev, int state){
  int t_max = prev.n_slices;
  vec out(t_max, fill::zeros);

  for(int i=0; i<(diag_time-1); ++i){
    out[i] = prev(diag_time-2,state-1,i); //-1 for different R/C++ conventions  //t-2 emulates what Paul did with adding zeros (prevalence 1 quarter before dx)
  }

  return out;
}

// [[Rcpp::export]]
mat diag_to_infec_wrap1(const cube& prev, int state){
  int t_max = prev.n_slices;
  mat out(t_max, t_max, fill::zeros);

  for(int i=(state-1); i<(t_max-1); ++i){  // -1 for different R/C++ conventions
    rowvec tmp = diag_to_infec_agg1(i+2, prev, state).t(); //-1 for different R/C++ conventions  //t-2 emulates what Paul did with adding zeros (prevalence 1 quarter before dx)
    out.row(i-state+1) = tmp/accu(tmp); // Normalization step
  }

  return out;
}


// The idea is the following
// m1[yr.inf,yr.dx] gives imputed yr.inf of AIDS dx
// so consider all yr.inf < t (1:y) and remove all dx that have happened
// before t (-1:-y)... these are the actually undiagnosed prevalent

// [[Rcpp::export]]
vec prev_fun(mat m1){
  int n_c = m1.n_cols-1;
  int n_r = m1.n_rows;
  vec out(n_r, fill::zeros);
  for(int i=0; i < n_r-1; ++i){
    out(i) = accu(m1.submat(0, i+1, i, n_c));
  }
  return out;
}

// The idea is the same as for prev_fun
// Just that these have cube
// because HIV diagnoses are further stratified by state

// [[Rcpp::export]]
vec prev_fun1(cube c1){
  int n_c = c1.n_cols-1;
  int n_r = c1.n_rows-1;
  int n_s = c1.n_slices-1;
  vec out(n_c+1, fill::zeros);
  for(int i=0; i < n_c; ++i){
    out(i) = accu( c1.subcube( 0, 0, i+1, n_r, i, n_s ) );
  }
  return out;
}
