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


// Pieces for the mixture approximations are shared between both prior specifications

// autodiff for a mixture normal q distribution

struct Qlogdens {
  const vec theta;
  const int mix;
  Qlogdens(const vec& thetaIn, const int& mixIn) :
    theta(thetaIn), mix(mixIn) {}
  template <typename T> //
  T operator ()(const Matrix<T, Dynamic, 1>& lambda)
    const{
    using std::log; using std::exp; using std::pow;
    
    int dim = theta.n_elem;
    
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

    Matrix<T, Dynamic, 1> weights(mix);
    T sumExpZ = 0;
    for(int k = 0; k < mix; ++k){
      weights(k) = exp(lambda(2*dim*mix + k));
      sumExpZ += weights(k);
    }
    weights /= sumExpZ;
    
    T density = 0;
    for(int k = 0; k < mix; ++k){
      density += weights(k) * dets(k) * pow(6.283185, -0.5 * dim) * exp(-0.5 * kernel(k));
    }
    
    T logMVN = log(density);
    
    return logMVN;
  }
};


// autodiff for a single MVN q distribution
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

// log joint distribution of theta and y, averaged over (prior) probabilities for each k
double pLogDens(mat y, vec theta, mat probK, mat mean, cube SigInv, vec weights, vec dets){
  int N = y.n_cols;
  int T = y.n_rows;
  int mix = probK.n_cols;
  
  int mixPrior = weights.n_rows;
  // Evaluate log(p(theta)), start by evaluative the quadratic in the MVN exponents
  double prior = 0;
  for(int k = 0; k < mixPrior; ++k){
    prior += weights(k) * pow(6.283185, -3) * dets(k) * exp(-0.5 * as_scalar((theta - mean.col(k)).t() * SigInv.slice(k) * (theta - mean.col(k))));
  }
  double logPrior = log(prior);
  
  double logLik = 0;
  for(int i = 0; i < N; ++i){
    for(int t = 0; t < T; ++t){
      double likelihood = 0;
      for(int m = 0; m < mix; ++m){
        likelihood += probK(i, m) *  pow(2 * 3.141598 * exp(theta(m)), -0.5) * exp(-pow(y(t, i) - theta(mix + m), 2) / (2 * exp(theta(m))));
      }
      logLik += log(likelihood);
    }
  }
  return logPrior + logLik;
}

// These models are parameterised by the mean and log standard deviations
// [[Rcpp::export]]
Rcpp::List mixtureNormal(mat data, Rcpp::NumericMatrix lambdaIn, mat probK, vec theta, mat mean, cube SigInv, vec weights, vec dets, int mix){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1> grad(dim);
  
  
  Qlogdens logQ(theta, mix);
  stan::math::set_zero_all_adjoints();
  stan::math::gradient(logQ, lambda, qEval, grad);
  double logp = pLogDens(data, theta, probK, mean, SigInv, weights, dets);
  double elbo = logp - qEval;
  for(int i = 0; i < dim; ++i){
    grad(i) *= elbo;
  }
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = elbo);
}

// Score derivative on its own, used in UVB-IS, either a single or mixture q distribution
// [[Rcpp::export]]
Rcpp::List mixNormScore(vec theta, Rcpp::NumericMatrix lambdaIn, int mix){
  Map<MatrixXd> lambda(Rcpp::as<Map<MatrixXd> >(lambdaIn));
  double qEval;
  int dim = lambda.rows();
  Matrix<double, Dynamic, 1> grad(dim);
  
  if(mix == 1){
    QlogdensSingle logQ(theta);
    stan::math::set_zero_all_adjoints();
    stan::math::gradient(logQ, lambda, qEval, grad);
    
  } else {
    Qlogdens logQ(theta, mix);
    stan::math::set_zero_all_adjoints();
    stan::math::gradient(logQ, lambda, qEval, grad);
    
  }
  
 
  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("val") = qEval);
}

// Metropolis Hastings Log Density Evaluation
// [[Rcpp::export]]
double MNLogDens(mat y, vec theta, vec K, vec mean, mat VarInv, int mix = 2){
  int N = y.n_cols;
  int T = y.n_rows;
  
  double logPrior = - as_scalar((theta - mean).t() * VarInv * (theta - mean));
  
  double logLik = 0;
  for(int i = 0; i < N; ++i){
    double mu = theta(mix + K(i));
    double var = exp(theta(K(i)));
    for(int t = 0; t < T; ++t){
      logLik += - 0.5 * log(var) - pow(y(t, i) - mu, 2) / (2 * var);
    }
  }
  return logLik + logPrior;
}

// log posterior of K
// [[Rcpp::export]]
vec postK(vec y, vec theta, vec kPrior){
  double p0 = kPrior(0);
  double p1 = kPrior(1);
  
  for(int t = 0; t < y.n_rows; ++t){
    p0 += - 0.5 * theta(0) - pow(y(t) - theta(2), 2) / (2 * exp(theta(0)));
    p1 += - 0.5 * theta(1) - pow(y(t) - theta(3), 2) / (2 * exp(theta(1)));
  }
  vec out {p0, p1};
  return out;
}


