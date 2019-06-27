#' Reparameterised VB gradient ascent.
#'
#' This function runs a Variational Bayes gradient ascent to find the parameters of an approximating distribution that minimises the KL divergence to the true posterior where the gradients are given by the reparameterised gradient estimator of Kingma and Ba 2014
#' Given a model function that calculates the derivative of the ELBO with resepct to lambda, it will return a list containing the optimal lambda,
#' history of ELBO values, and iterations required to converge for a given model function and dataset.
#' @param data Data passed directly into model, may also be used in conjunction with batch.
#' @param lambda Parameter vector or one column maxtrix of auxiliary parameters to be optimised.
#' @param model A function with first three arguments: data, lambda, epsilon, that returns a list with elements 'grad', the vector of the gradient of the ELBO with respect to lambda, and 'val', the value of the ELBO for this lambda.
#' @param dimTheta Integer. The dimension of the theta parameter vector
#' @param S Integer. The number of Monte Carlo estimates per iteration. Defaults to 25.
#' @param epsDist Either 'normal' or 'uniform'. Distribution of the auxiliary random variable. Defaults to 'normal'
#' @param batch Integer. If data is a matrix, calculate gradient based on batch many columns per iteration, if ncol(data) is divisible by batch. If batch = 0 use the whole data object per iteration.
#' @param maxIter Interger. Maximum number of gradient ascent iterations. Defaults to 5000.
#' @param alpha adam optimisation control parameter. Defaults to 0.01.
#' @param beta1 adam optimisation control parameter. Defaults to 0.9.
#' @param beta2 adam optimisation control parameter. Defaults to 0.99.
#' @param rollingWindowSize Integer. Take the mean of this many most recent iterations to assess convergence. Defaults to 5 if batch = 0, or ncol(data) / batch otherwise.
#' @param threshold Maximum difference in mean value of ELBO before convergence is achieved. Defaults to 0.01.
#' @param suppressProgress Boolean, if true the program will not print ELBO values to console during optimisation.
#' @param ... Extra arguments passed into model
#' @export
reparamVB <- function(data, lambda, model, dimTheta, S = 25, epsDist = 'normal', batch = 0, maxIter = 5000, alpha = 0.01, beta1 = 0.9, beta2 = 0.99, rollingWindowSize = 5, threshold = 0.01, suppressProgress = FALSE, ...){
  if(!is.matrix(lambda)){
    lambda <- matrix(lambda, ncol = 1)
  }
   dimLambda <- length(lambda)
  

  if(batch > 0 & is.matrix(data)){
    rollingWindowSize <- ncol(data) / batch
    if(!is.integer(rollingWindowSize)){
      stop('Number of columns of data must be divisible by batch')
    }
    subset <- c(0, seq(batch, ncol(data), batch))
    dataFull <- data
  } else if(batch > 0 & !is.matrix(data)){
    stop('batch > 0 is only available for matrix data')
  }

  if(!epsDist %in% c('normal', 'uniform')){
    stop('epsDist must be normal or uniform')
  }

  diff <- threshold + 1
  iter <- 1
  LB <- numeric(maxIter)
  M <- V <- numeric(dimLambda)
  e <- 1e-8
  meanLB <- 0
  oldMeanLB <- 0
  while(diff > threshold){
    if(iter > maxIter){
      break
    }
    if(batch > 0){
      set <- iter %% rollingWindowSize + 1
      data <- dataFull[,(subset[set]+1):subset[set+1]]
    }
    eval <- numeric(S)
    grad <- matrix(0, dimLambda, S)
    q <- numeric(S)

    unif <- matrix(runif(S*dimTheta), ncol = dimTheta)

    if(S == 1){
      epsilon <- unif[1,]
      q <- 0
      if(epsDist == 'normal'){
        epsilon <- qnorm(epsilon)
        q <- sum(dnorm(epsilon, log = TRUE))
      }
      logpj <- model(data, lambda, epsilon, ...)
      eval <- logpj$val
      grad <- logpj$grad
      gradient <- grad
      gradientSq <- grad^2
      LB[iter] <- eval - q
    } else {
      for(s in 1:S){
        epsilon <- unif[s,]
        q <- 0
        if(epsDist == 'normal'){
          epsilon <- qnorm(epsilon)
          q <- sum(dnorm(epsilon, log = TRUE))
        }
        logpj <- model(data, lambda, epsilon, ...)   
        eval[s] <- logpj$val
        grad[,s] <- logpj$grad
      }
      eval[eval == -Inf] = NA
      gradient <- rowMeans(grad, na.rm = TRUE)
      gradientSq <- rowMeans(grad^2, na.rm = TRUE)
      LB[iter] <- mean(eval - q, na.rm=TRUE) 
    }
   
    M <- beta1 * M + (1 - beta1) * gradient
    V <- beta2 * V + (1 - beta2) * gradientSq
    Mstar <- M / (1 - beta1^iter)
    Vstar <- V / (1 - beta2^iter)
    update <- alpha * Mstar / (sqrt(Vstar) + e)
    if(any(is.na(update))){
      print('Break - NA gradient value')
      break
    }
    lambda <- lambda + update
  
    if(iter %% rollingWindowSize == 0){
      oldMeanLB <- meanLB
      meanLB <- mean(LB[iter:(iter-(rollingWindowSize - 1))])
      diff <- abs(meanLB - oldMeanLB)
    }
    if(iter %% 100 == 0 & !suppressProgress){
        print(paste0('Iteration: ', iter, ' ELBO: ', meanLB))
    }
    iter <- iter + 1
  }
  iter <- iter - 1
  if(!suppressProgress){
    print(paste0('iter: ', iter, ' ELBO: ', meanLB))
  }
  return(list(lambda=lambda,
              LB = LB[1:iter], 
              iter = iter))
}
