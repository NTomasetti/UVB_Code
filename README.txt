This repository contains the R code used for Updating Variational Bayes.

The code to implement a single simulation for the AR3 model and Normal Mixture Model in Section 5, as well as the applications in Section 6 and Section 7 can be found in the relevant files named 'Main code for... .R'. 

Implementation of the DPM model requires the NGSIM dataset, which can be downloaded from a link in the references of the submitted paper. There is a data pre-processing step to estimate the vehicle midlines, which can be found in carDataTransform.R

Each of the four main files requires functions defined in code contained in the folder Functions. Some of these files are specific to a particular application, and named accordingly, though there is also five generic functions: gaussianDensity.R, metropolisHastings.R, reparamVB.R, scoreVB.R, and UVBIS.R.

gaussianDensity.R will input a vector of means, a vector of standard deviations, and a vector of transformations from 'identity', 'exp', 'sigmoid', 'stretched sigmoid' and calculate the density of y = f(x), where x is normal and f() is the given function:
'identity' - f(x) = x
'exp' - f(x) = exp(x)
'sigmoid' - f(x) = 1 / (1 + exp(-x))
'stretched sigmoid' - f(x) = 2 / (1 + exp(-x)) - 1
These can be used to calculate the density of distributions that are defined by transformations of the standard normal, as used in reparameterised gradient VB. 

metropolisHastings.R contains a generic framework to run a Metropolis Hastings MCMC algorithm with a Normal Random Walk proposal distribution. The user will need to supply a function to calculate the log joint density of the model - log(p(y, theta)).

reparamVB.R contains a generic framework to run a reparamaterised gradient VB optimisation for either a uniform or standard normal pre-transformation approximating distribution. The user will need to supply a function that calculates the actual gradient of the ELBO as a function of the data, parameter vector lambda, and auxiliary random variables epsilon.

scoreVB.R contains a similar framework to run VB optimisation with the score gradients. It supports use of a mixture of diagonal covariance multivariate normals as the approximating distribution, and adds control variates to reduce the gradient variance. The user will need to supply a function that calculates the actual gradient of the ELBO as a function of the data, parameter vector lambda, and simulated random variables theta

UVBIS.R is another generic framework to run UVB-IS optimisation. Most of the calculations for UVB-IS can be carried out before running the optimisation (i.e. simulating values of theta from the proposal, calculating the density of these values with respect to the proposal distribution, and calculating the log joint distribution of the data and simualted theta values, log(p(y, theta)). Given these values ('samples', 'dSamples', and 'logjoint'), UVB-IS just requires the score function of the new approximating distribution, which must be provided as an argument 'qDist'. See implementations of UVB-IS in `Main code for AR3 Simulation.R` and `Main code for Mixture Simulation.R` for examples.


