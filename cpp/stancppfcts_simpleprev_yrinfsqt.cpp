// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using std::vector;
using namespace arma;

// [[Rcpp::export]]
mat Q1_fct(const mat& d, const mat& q, int t, int a){
  mat Q(4,4);
  Q(0,0) = (1 - d(0,t)) * (1 - q(0,a)) ;
  Q(0,1) = (1 - d(0,t)) * q(0,a) ;
  Q(0,2) = 0 ;
  Q(0,3) = 0 ;
  Q(1,0) = 0 ;
  Q(1,1) = (1 - d(1,t)) * (1 - q(1,a)) ;
  Q(1,2) = (1 - d(1,t)) * q(1,a) ;
  Q(1,3) = 0 ;
  Q(2,0) = 0 ;
  Q(2,1) = 0 ;
  Q(2,2) = (1 - d(2,t)) * (1 - q(2,a)) ;
  Q(2,3) = (1 - d(2,t)) * q(2,a) ;
  Q(3,0) = 0 ;
  Q(3,1) = 0 ;
  Q(3,2) = 0 ;
  Q(3,3) = (1 - d(3,t)) * (1 - q(3,a)) ;
  return Q;
}

// [[Rcpp::export]]
mat P1_fct(const mat& Q, const mat& d, const mat& q, int t, int a, int l){
  mat Pl(4,4);
  mat PA(4,5);
  mat I4 = eye(4,4); //variable declaration for identity mat
  mat zero_mat(4, 4, fill::zeros);
  
  if(l==0){
    Pl = I4 + Q + (Q * Q) + (Q * Q * Q);
  }
  else if(l==1){
    Pl = I4 + Q + (Q * Q);
  }
  else if(l==2){
    Pl = I4 + Q ;
  }
  else if(l==3){
    Pl = I4 ;
  }
  else{
    Pl = zero_mat ;
  }
  
  for(int i=0; i<4; ++i){
    for(int k=0; k<4; ++k){
      // HIV diagnosis (k=1:4)
      PA(i,k) = Pl(i,k) * d(k,t) ;
    }
    // Progression to AIDS (k=5)
    PA(i,4) = Pl(i,3) * (1 - d(3,t)) * q(3,a) ;
  }
  return PA;
}

// Function to make any [t][a0][t0] array 
// into a [t][a][t0] array

vector<vector<vector<arma::vec> > > a0toa_fct(const vector<vector<vector<arma::vec> > >& array, int nyr, int nage){
  vector<vector<vector<arma::vec> > > out(nyr, vector<vector<arma::vec> >(nage, vector<arma::vec>(nyr)));
  
  vec emptyvec(array[0][0][0].n_elem,fill::zeros);
  
  for(int t=0; t<nyr; ++t){
    for(int a =0; a<nage; ++a){
      for(int t0=0; t0<nyr; ++t0){
        out[t][a][t0] = emptyvec;
      }
    }
  }
  
  for(int t0=0; t0<nyr; ++t0){
    for(int a0 =0; a0<nage; ++a0){
      for(int t=t0; t<nyr; ++t){
        int a = std::min(a0+t-t0,nage-1);
        // Rcpp::Rcout << "a:" << std::endl << a << std::endl;
        out[t][a][t0] += array[t][a0][t0];
      }
    }
  }
  
  return out;
}


// [[Rcpp::export]]
cube gof_iter_fct(const vec& H, const mat& d, const mat& q, const mat& prev, int nyr, int nage){
  vector<vector<arma::mat> > QA(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > QA1(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > QA2(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > QA3(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > PA(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > PA1(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > PA2(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > PA3(nyr, vector<arma::mat>(nage));
  vector<vector<vector<arma::vec> > > lat_arr(nyr, vector<vector<arma::vec> >(nage, vector<arma::vec>(nyr)));
  vector<vector<vector<arma::vec> > > diag_arr(nyr, vector<vector<arma::vec> >(nage, vector<arma::vec>(nyr)));
  vector<vector<vector<arma::vec> > > diag_arr1(nyr, vector<vector<arma::vec> >(nage, vector<arma::vec>(nyr)));
  
  cube exp_diag(nyr,nage,5);
  const vec empty_vec(5,fill::zeros);
  const vec empty_vec1(4,fill::zeros);
  
  for(int t=0; t<nyr; ++t){
    for(int a=0; a<nage; ++a){
      mat Q = Q1_fct(d,q,t,a);
      QA[t][a] = Q * Q * Q * Q; // Infection in previous year... 4 transitions
      QA1[t][a] = Q * Q * Q; // Infection in quarter 1 ... 3 transitions
      QA2[t][a] = Q * Q; // Infection in quarter 2 ... 2 transitions
      QA3[t][a] = Q; // Infection in quarter 3 ... 1 transitions
      PA[t][a] = P1_fct(Q, d, q, t, a, 0); // Infection in previous year... 4 transitions
      PA1[t][a] = P1_fct(Q, d, q, t, a, 1); // Infection in quarter 1 ... 3 transitions
      PA2[t][a] = P1_fct(Q, d, q, t, a, 2); // Infection in quarter 2 ... 2 transitions
      PA3[t][a] = P1_fct(Q, d, q, t, a, 3); // Infection in quarter 3 ... 1 transitions
      // If infection in last quarter, by assumption no diagnoses
    }
  }
  

  for(int t0=0; t0<nyr; ++t0){
    for(int a0=0; a0<nage; ++a0){
      // Adding initial prevalence at t==0
      vec diag_p(5, fill::zeros);
      vec lat_p(4, fill::zeros);
      vec infs4(4, fill::zeros);
      vec tmp(4);

      if(t0==0){
        tmp(0) = prev(0,a0);
        tmp(1) = prev(1,a0);
        tmp(2) = prev(2,a0);
        tmp(3) = prev(3,a0);

        // Progressed initial prevalence
        lat_p = lat_p + (QA[0][a0].t() * tmp);

        // Diagnosed initial prevalence
        diag_p = diag_p + (PA[0][a0].t() * tmp);

      }
      else{ // Probably dont need this
        lat_p = empty_vec1;
        diag_p = empty_vec;
      }

      // This initialization step is needed
      lat_arr[t0][a0][t0] = empty_vec1;
      diag_arr[t0][a0][t0] = empty_vec;
      
      // Rcpp::Rcout << "lat_p:" << std::endl << lat_p << std::endl;
      // Rcpp::Rcout << "dim:" << std::endl << lat_arr[t0][a0][t0].n_elem << std::endl;
      infs4[0] = H((a0*nyr)+t0)/4;

      lat_arr[t0][a0][t0] = lat_p + (QA1[t0][a0].t() * infs4 + QA2[t0][a0].t() * infs4 + QA3[t0][a0].t() * infs4 + infs4);
      
      diag_arr[t0][a0][t0] = diag_p + (PA1[t0][a0].t() * infs4 + PA2[t0][a0].t() * infs4 + PA3[t0][a0].t() * infs4);

      for(int t=(t0+1); t<nyr; ++t){

        // Rcpp::Rcout << "t:" << std::endl << t << std::endl;

        lat_arr[t][a0][t0] = QA[t][a0].t()* lat_arr[t-1][a0][t0];

        diag_arr[t][a0][t0] = PA[t][a0].t()* lat_arr[t-1][a0][t0];

      }
    }
  }
  
  diag_arr1 = a0toa_fct(diag_arr, nyr, nage);

  for(int t=0; t<nyr; ++t){
    for(int a=0; a<nage; ++a){
      exp_diag(t,a,0) =  0;
      exp_diag(t,a,1) =  0;
      exp_diag(t,a,2) =  0;
      exp_diag(t,a,3) =  0;
      exp_diag(t,a,4) =  0;

      for(int t0=0; t0<=t; ++t0){
        exp_diag(t,a,0) +=  diag_arr1[t][a][t0][0];
        exp_diag(t,a,1) +=  diag_arr1[t][a][t0][1];
        exp_diag(t,a,2) +=  diag_arr1[t][a][t0][2];
        exp_diag(t,a,3) +=  diag_arr1[t][a][t0][3];
        exp_diag(t,a,4) +=  diag_arr1[t][a][t0][4];
      }
    }
  }
  
  return exp_diag;
}


// [[Rcpp::export]]
cube prev_iter_fct(const vec& H, const mat& d, const mat& q, const mat& prev, int nyr, int nage){
  vector<vector<arma::mat> > QA(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > QA1(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > QA2(nyr, vector<arma::mat>(nage));
  vector<vector<arma::mat> > QA3(nyr, vector<arma::mat>(nage));
  vector<vector<vector<arma::vec> > > lat_arr(nyr, vector<vector<arma::vec> >(nage, vector<arma::vec>(nyr)));
  vector<vector<vector<arma::vec> > > lat_arr1(nyr, vector<vector<arma::vec> >(nage, vector<arma::vec>(nyr)));
  
  cube exp_diag(nyr,nage,4);
  const vec empty_vec1(4,fill::zeros);
  
  for(int t=0; t<nyr; ++t){
    for(int a=0; a<nage; ++a){
      mat Q = Q1_fct(d,q,t,a);
      QA[t][a] = Q * Q * Q * Q; // Infection in previous year... 4 transitions
      QA1[t][a] = Q * Q * Q; // Infection in quarter 1 ... 3 transitions
      QA2[t][a] = Q * Q; // Infection in quarter 2 ... 2 transitions
      QA3[t][a] = Q; // Infection in quarter 3 ... 1 transitions
      // If infection in last quarter, by assumption no diagnoses
    }
  }
  
  for(int t0=0; t0<nyr; ++t0){
    for(int a0=0; a0<nage; ++a0){
      // Adding initial prevalence at t==0
      vec diag_p(5, fill::zeros);
      vec lat_p(4, fill::zeros);
      vec infs4(4, fill::zeros);
      vec tmp(4);
      
      // This initialization step is needed
      lat_arr[t0][a0][t0] = empty_vec1;

      if(t0==0){
        tmp(0) = prev(0,a0);
        tmp(1) = prev(1,a0);
        tmp(2) = prev(2,a0);
        tmp(3) = prev(3,a0);
        
        // Progressed initial prevalence
        lat_p = lat_p + (QA[0][a0].t() * tmp);
        
      }
      else{ // Probably dont need this
        lat_p = empty_vec1;
      }
      
      
      // Rcpp::Rcout << "lat_p:" << std::endl << lat_p << std::endl;
      // Rcpp::Rcout << "dim:" << std::endl << lat_arr[t0][a0][t0].n_elem << std::endl;
      infs4[0] = H((a0*nyr)+t0)/4;

      lat_arr[t0][a0][t0] = lat_p + (QA1[t0][a0].t() * infs4 + QA2[t0][a0].t() * infs4 + QA3[t0][a0].t() * infs4 + infs4);

      for(int t=(t0+1); t<nyr; ++t){

        // Rcpp::Rcout << "t:" << std::endl << t << std::endl;

        lat_arr[t][a0][t0] = QA[t][a0].t()* lat_arr[t-1][a0][t0];
        
      }
    }
  }
  
  lat_arr1 = a0toa_fct(lat_arr, nyr, nage);

  for(int t=0; t<nyr; ++t){
    for(int a=0; a<nage; ++a){
      exp_diag(t,a,0) =  0;
      exp_diag(t,a,1) =  0;
      exp_diag(t,a,2) =  0;
      exp_diag(t,a,3) =  0;

      for(int t0=0; t0<=t; ++t0){
        exp_diag(t,a,0) +=  lat_arr1[t][a][t0][0];
        exp_diag(t,a,1) +=  lat_arr1[t][a][t0][1];
        exp_diag(t,a,2) +=  lat_arr1[t][a][t0][2];
        exp_diag(t,a,3) +=  lat_arr1[t][a][t0][3];
      }
    }
  }
  
  return exp_diag;
}
