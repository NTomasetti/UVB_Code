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

using namespace arma;
using namespace boost::math;
using namespace std;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd; 
using Eigen::Map;

// Log Joint Density (ie log prior + log likelihood) Evaluations used

// Evaluate the log joint density for a single draw of theta (in MCMC)
// [[Rcpp::export]]
double ARjointDensSingle (vec x, vec theta, vec mu, mat varInv, int lags){
  
  double dens = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta - mu).t() * varInv * (theta - mu));
  for(int t = lags; t < x.n_elem; ++t){
    double mean = theta(1);
    for(int j = 1; j <= lags; ++j){
      mean += theta(1 + j) * (x(t - j) - theta(1));
    }
    dens += -0.5 * theta(0) - pow(x(t) - mean, 2) / (2 * exp(theta(0)));
  }
  return dens;
}

// As ARjointDensSingle, but for a mixture prior distribution, used as part of the score gradient
// [[Rcpp::export]]
double ARjointDensMix(vec x, vec theta, mat mean, cube SigInv, vec dets, vec weights, int lags){
  
  int mix = weights.n_rows;
  // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
  double prior = 0;
  for(int k = 0; k < mix; ++k){
    prior += weights(k) * pow(6.283185, -3) * dets(k) * exp(-0.5 * as_scalar((theta - mean.col(k)).t() * SigInv.slice(k) * (theta - mean.col(k))));
  }
  // Evaluate log likelihood
  double logLik = 0;
  for(int t = lags; t < x.n_elem; ++t){
    double mean = theta(1);
    for(int j = 1; j <= lags; ++j){
      mean += (x(t - j) - theta(1)) * theta(1 + j);
    }
    logLik += - 0.5 * theta(0) - pow(x(t) - mean, 2) / (2 * exp(theta(0)));
  }
  
  return std::log(prior) + logLik;
}

// Evalaute the log joint density for a matrix of theta draws, used for UVB-IS
// [[Rcpp::export]]
vec ARjointDensSingleMatrix (vec x, mat theta, vec mu, mat varInv, int lags){
  int N = theta.n_rows;
  
  vec dens(N);
  
  for(int i = 0; i < N; ++i){
    dens(i) = 0.5 * log(det(varInv)) - 0.5 * as_scalar((theta.row(i).t() - mu).t() * varInv * (theta.row(i).t() - mu));
    for(int t = lags; t < x.n_elem; ++t){
      double mean = theta(i, 1);
      for(int j = 1; j <= lags; ++j){
        mean += x(t - j) * theta(i, 1 + j);
      }
      dens(i) += -0.5 * theta(i, 0) - pow(x(t) - mean, 2) / (2 * exp(theta(i, 0)));
    }
  }
  dens -= max(dens);
  return exp(dens);
}

// As ARjointDens but setup for a mixture prior distribution, used for UVB-IS
// [[Rcpp::export]]
vec ARjointDensMixMatrix(vec x, mat theta, mat mean, cube SigInv, vec dets, vec weights, int lags){
  
  int mix = weights.n_rows;
  int N = theta.n_rows;
  vec dens(N);
  for(int i = 0; i < N; ++i){
    // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
    double prior = 0;
    for(int k = 0; k < mix; ++k){
      prior += weights(k) * pow(6.283185, -3) * dets(k) * exp(-0.5 * as_scalar((theta.row(i).t() - mean.col(k)).t() * SigInv.slice(k) * (theta.row(i).t() - mean.col(k))));
    }
    // Evaluate log likelihood
    double logLik = 0;
    for(int t = lags; t < x.n_elem; ++t){
      double mean = theta(i, 1);
      for(int j = 1; j <= lags; ++j){
        mean += (x(t - j) - theta(i, 1)) * theta(i, 1 + j);
      }
      logLik += - 0.5 * theta(i, 0) - pow(x(t) - mean, 2) / (2 * exp(theta(i, 0)));
    }
    dens(i) = log(prior) + logLik;
  }
  dens -= max(dens);
  return exp(dens);
}




// VB Score Gradient components, called from AR3VB and AR3Score functions defined later

// Some pieces for the mixture approximations are shared between both prior specifications
// This takes the derivative of log q when q is a mixture of multivariate normals
struct mixQlogdens {
  const vec theta;
  const int mix;
  mixQlogdens(const vec& thetaIn, const int& mixIn) :
    theta(thetaIn), mix(mixIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow; using std::sqrt;
    
    int dim = theta.n_rows;
    
    Matrix<T, Dynamic, 1> dets(mix);
    for(int k = 0; k < mix; ++k){
      dets(k) = exp(lambda(dim*k + dim*mix));
      for(int i = 1; i < dim; ++i){
        dets(k) *= exp(lambda(dim*mix + dim*k + i));
      }
      dets(k) = 1.0 / dets(k);
    }
    
    Matrix<T, Dynamic, 1> kernel(mix);
    kernel.fill(0);
    
    for(int k = 0; k < mix; ++k){
      for(int i = 0; i < dim; ++i){
        kernel(k) += pow((theta(i) - lambda(k*dim + i)) / exp(lambda(dim*mix + dim*k + i)), 2);
      }
    }
    
    Matrix<T, Dynamic, 1> pi(mix);
    T sumExpZ = 0;
    for(int k = 0; k < mix; ++k){
      pi(k) = exp(lambda(2*dim*mix + k));
      sumExpZ += pi(k);
    }
    pi /= sumExpZ;
    
    T density = 0;
    for(int k = 0; k < mix; ++k){
      density += pi(k) * dets(k) * pow(6.283185, -3) *  exp(-0.5 * kernel(k));
    }
    
    return log(density);
  }
};

// This takes the derivative of log q when q is a single component multivariate normal
struct QlogdensSingle {
  const vec theta;
  QlogdensSingle(const vec& thetaIn) :
    theta(thetaIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow;
    
    int dim = theta.n_elem;
    T logMVN = 0;
    
    for(int i = 0; i < dim; ++i){
      logMVN += - 0.5 * log(2 * 3.14159) - lambda(dim + i) -  0.5 * pow((theta(i) - lambda(i)) / exp(lambda(dim + i)), 2);
    }
    return logMVN;
  }
};

// Main function for score gradients (in SVB and UVB) for a single draw of theta, will:
// 1) Take the derivative of log q, when q is either a single or mixture of MVNs (controlled by mix argument)
// 2) Evalaute log(p(theta, y)), which works fine for either a mixture prior or single component - number of components controlled by length of weights 
// 3) Combine these to get the score gradient
// [[Rcpp::export]]
Rcpp::List AR3VB(vec data, Rcpp::NumericMatrix lambdaIn, vec theta, mat mean, cube SigInv, vec dets, vec weights, int mix = 1){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1> grad(dim);
  
  if(mix == 1){
    QlogdensSingle logQ(theta);
    stan::math::set_zero_all_adjoints();
    stan::math::gradient(logQ, lambda, qEval, grad);
    
  } else {
    mixQlogdens logQ(theta, mix);
    stan::math::set_zero_all_adjoints();
    stan::math::gradient(logQ, lambda, qEval, grad);
  }
  
  double logp = ARjointDensMix(data, theta, mean, SigInv, dets, weights, 3);
  double elbo = logp - qEval;
  for(int i = 0; i < dim; ++i){
    grad(i) *= elbo;
  }
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = elbo);
}

// As above but only deals with d logq / dlambda, used for UVB-IS
// [[Rcpp::export]]
Rcpp::List AR3Score(vec theta, Rcpp::NumericMatrix lambdaIn, int mix = 1){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1> grad(dim);
  
  if(mix == 1){
    QlogdensSingle logQ(theta);
    stan::math::set_zero_all_adjoints();
    stan::math::gradient(logQ, lambda, qEval, grad);
    
  } else {
    mixQlogdens logQ(theta, mix);
    stan::math::set_zero_all_adjoints();
    stan::math::gradient(logQ, lambda, qEval, grad);
    
  }

  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = qEval);
}

