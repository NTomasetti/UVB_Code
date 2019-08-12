// [[Rcpp::depends(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <boost/math/distributions.hpp> 
#include <boost/math/special_functions.hpp> 


using namespace arma;
using namespace boost::math;
using namespace std;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// Autodiff for a single component MVN distribution

struct QlogdensSingle {
  const mat theta;
  QlogdensSingle(const vec& thetaIn) :
    theta(thetaIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow;
    
    int dim = theta.n_cols;
    T logMVN = 0;
    
    for(int i = 0; i < dim; ++i){
      
      // This autodiff library doesn't work well with matrix operations, so we will have to hard code all of the steps to get from L to Sigma^-1, and the matrix multiplication
      // It is all 2x2 so it isn't difficult to do
      T L11 = lambda(5*i + 2);
      T L21 = lambda(5*i + 3);
      T L22 = lambda(5*i + 4);
      
      // Variance matrix
      T sigma11 = L11 * L11;
      T sigma21 = L11 * L21;
      T sigma22 = L21 * L21 + L22 * L22
      
      // Inverse variance matrix
      T det = sigma11 * sigma22 - sigma21 * sigma21;
      T sigInv11 = sigma22 / det;
      T sigInv22 = sigma11 / det;
      T sigInv21 = - sigma21 / det;
      
      // Calculate bivariate normal log density - again manual matrix operation.
      
      logMVN += - log(2 * 3.14159) - 0.5 * (log(L11) + log(L22))  -  0.5 * 
        ((theta(0, i) - lambda(5*i)) * ((theta(0, i) - lambda(5*i)) * sigInv11 + (theta(1, i) - lambda(5*i+1)) * sigInv21) + 
        (theta(1, i) - lambda(5*i+1)) * ((theta(0, i) - lambda(5*i)) * sigInv21 + (theta(1, i) - lambda(5*i+1)) * sigInv22));
    }
    return logMVN;
  }
};

// Calculate the scroe gradient (for a given value log p)

// [[Rcpp::export]]
Rcpp::List mixtureNormalDPM(mat theta, Rcpp::NumericMatrix lambdaIn, double logp){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1> grad(dim);
  
  QlogdensSingle logQ(theta);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(logQ, lambda, qEval, grad);

  double elbo = logp - qEval;
  for(int i = 0; i < dim; ++i){
    grad(i) *= elbo;
  }
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = elbo);
}


