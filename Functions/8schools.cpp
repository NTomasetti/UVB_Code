// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace arma;
using namespace boost::math;
using namespace std;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// Autodiff for log p(theta, mu, tau, y) for the reparam gradient estimator at the intial time
struct schoolsInitial {
  const vec epsilon;
  const vec ybar;
  const vec sigmaj;
  const int J;
  schoolsInitial(const vec& epsIn, const vec& ybarIn, const vec& sigmajIn, const int& jIn) :
    epsilon(epsIn), ybar(ybarIn), sigmaj(sigmajIn), J(jIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::lgamma;
    
    int dim = J + 2;
    double nu = 4;
    
    Matrix<T, Dynamic, 1> theta(dim);
    
    for(int i = 0; i < dim; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(10*(i+1) + j) * epsilon(j);
      }
    }
    T tau = exp(theta(0));
    T mu = theta(1);
    
    T logdetJ = theta(0);
    for(int i = 0; i < dim; ++i){
      logdetJ += log(fabs(lambda(10*(i+1) + i)));
    }
    
    T loglik = 0;
    for(int i = 0; i < J; ++i){
      loglik += -log(tau) - 0.5 * (nu + 1) * log(1 + pow((theta(2+i) - mu)/tau, 2) / nu);
      loglik += -log(sigmaj(i)) - pow(ybar(i) - theta(2+i), 2) / (2 * sigmaj(i) * sigmaj(i));
    }
 
    return loglik + logdetJ;
  }
};

// Autodiff for log p(theta, mu, tau, y) for the reparam gradient estimator at the update
struct schoolsUpdate {
  const vec epsilon;
  const vec ybar;
  const vec sigmaj;
  const int J;
  const vec mean;
  const mat Linv;
  schoolsUpdate(const vec& epsIn, const vec& ybarIn, const vec& sigmajIn, const int& jIn, const vec& meanIn, const mat& LinvIn) :
    epsilon(epsIn), ybar(ybarIn), sigmaj(sigmajIn), J(jIn), mean(meanIn), Linv(LinvIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt; using std::lgamma;
    
    int oldJ = mean.n_elem - 2;
    int dim = 2 + J;
    double nu = 4;
    
    Matrix<T, Dynamic, 1> theta(dim);
    
    for(int i = 0; i < dim; ++i){
      theta(i) = lambda(i);
      for(int j = 0; j <= i; ++j){
        theta(i) += lambda(10*(i+1) + j) * epsilon(j);
      }
    }
    T tau = exp(theta(0));
    T mu = theta(1);
    
    T logdetJ = theta(0);
    for(int i = 0; i < dim; ++i){
      logdetJ += log(fabs(lambda(10*(i+1) + i)));
    }
    
    T prior = 0;
    Matrix<T, Dynamic, 1> kernel(dim);
    kernel.fill(0);
    for(int i = 0; i < mean.n_elem; ++i){
      for(int j = 0; j <= i; ++j){
        kernel(i) += (theta(j) - mean(j)) * Linv(i, j);
      }
      prior += - 0.5 * pow(kernel(i), 2);
    }
    
    T loglik = 0;
    for(int i = oldJ; i < J; ++i){
      loglik += -log(tau) - 0.5 * (nu + 1) * log(1 + pow((theta(2 + i) - mu)/tau, 2) / nu);
      loglik += -log(sigmaj(i)) - pow(ybar(i) - theta(2 + i), 2) / (2 * sigmaj(i) * sigmaj(i));
    }
    
    
    return loglik + logdetJ + prior;
  }
};



