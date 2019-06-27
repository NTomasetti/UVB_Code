rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(rstan)
library(tidyverse)
source('Functions/DPM funs.R')
source('Functions/gaussianDensity.R')
source('Functions/DPM CAVI.R')
sourceCpp('Functions/DPM funs.cpp')
data <- readRDS('subCarsUVB.RDS') # This data is not on Github (too large) but can be re-created from carDataTransform.R if data is downloaded from the NGSIM link provided in that file.

# Set N / T
N <- 500
T <- 200

y <- data[1:T, 1:N]

# Control amount of data used per update

initialBatch <- 50
updateSize <- 25
batches <- c(0, seq(initialBatch, T, updateSize))

# Initialise params

lambda <- matrix(0, 5, N)
# Mean of log sd
lambda[1,] <- log(apply(y, 2, sd))
# Mean of mu
lambda[2,] <- colMeans(y)
# log sd of above parameters
lambda[c(3, 5), ] <- 1

# Prior params N(0, 10)
pMean <- matrix(0, 2, N)
pSig <-  array(rep(c(10 * diag(2)), N),c(2, 2, N))

# Save values
lamb <- list()
pk <- list()
set.seed(10)

# Initial UVB fit

UVBfit <- fitVBDPM(data = y[1:batches[2], 1:N],
                 lambda = lambda, 
                 S = 100,
                 priorMean = pMean,
                 priorSig = pSig,
                 threshold = 0.001 * N * initialBatch,
                 alphaDPM = 1,
                 alpha = 0.2,
                 maxIter = 5000)

lamb[[1]] <- UVBfit$lambda

# Marginal probabilities of each k
probs <- VBDPMprobs(y[1:batches[2], 1:N],
                    lambda = UVBfit$lambda, 
                    alphaDPM = 1,
                    S = 50)
pk[[1]] <- probs


# UVB updates

for(t in 3:length(batches)){
  print(paste(t, Sys.time()))
  updateMean <- UVBfit$lambda[1:2,]
  updateSig <- priorSig
  for(i in 1:N){
    l <- matrix(c(UVBfit$lambda[c(3, 4), i], 0, UVBfit$lambda[5, i]), 2)
    updateSig[,,i] <- l %*% t(l)
  }
 
  UVBfit <- fitVBDPP(data = y[(batches[t-1]+1):batches[t], 1:N],
                         lambda = UVBfit$lambda, 
                         S = 100,
                         priorMean = updateMean,
                         priorSig = updateSig,
                         logPrior = probs,
                         threshold = 0.001 * N * updateSize,
                         alpha = 0.025,
                         alphaDPM = 1,
                         maxIter = 5000)
  lamb[[t-1]] <- UVBfit$lambda
  
  probs <- VBDPPprobs(data = y[1:batches[t], ],
                      lambda = UVBfit$lambda, 
                      alphaDPM = 1,
                      S = 50)
  pk[[t-1]] <- probs

}

# Also add SVB, SVB at T_1 = UVB Init

lambSVB <- list()
lambSVB[[1]] <- lamb[[1]]

for(t in 3:length(batches)){
  print(paste(t, Sys.time()))
  lambda <- matrix(0, 5, N)
  # Mean of log sd
  lambda[1,] <- log(apply(y, 2, sd))
  # Mean of mu
  lambda[2,] <- colMeans(y)
  # log sd of above parameters
  lambda[c(3, 5), ] <- 1
  
  
  SVBfit <- fitVBDPM(data = y[1:batches[t], 1:N],
                     lambda = lambda, 
                     S = 100,
                     priorMean = pMean,
                     priorSig = pSig,
                     threshold = 0.001 * N * initialBatch,
                     alphaDPM = 1,
                     alpha = 0.2,
                     maxIter = 5000)
  lambSVB[[t-1]] <- SVBfit$lambda
}


##### Code for assorted plots of the UVB results


# Set parameters to UVB T_6 values
UVBfit$lambda <- lamb[[t-1]]

# Sample

sample <- VBDPMsample(y[1:batches[t], 1:N],
                  lambda = UVBfit$lambda, 
                  alphaDPM = 1,
                  S = 100)

# Check which values are sampled at least once - ignore clusters with zero samples values (and zero weighted density)
amount <- table(c(sample$k))
include <- as.numeric(names(amount[amount > 1]))

# For each sampled value we should calculate the density of the parameter
dens <- NULL
for(i in 1:length(include)){
  # Extract parameters for this value
  lam <- c(UVBfit$lambda[, include[i]])
  
  weight <- sum(sample == include[i]) / length(sample)
  
  sd <- sqrt(c(lam[3]^2, lam[4]^2 + lam[5]^2))
  
  gaussianDensity(mu = lam[1:2],
                  sd = sd,
                  transform = c('exp', 'identity'),
                  names = c('sigma^{2}', 'mu')) %>%
        mutate(dens = density * weightsP[1]) %>%
        select(-density) -> df
  
  df$group <- include[i]
  df$weight <- weight
  
  dens <- rbind(dens, df %>% select(var, support, dens, group, weight))
}

# Make the first panel of the figure
dens %>%
  filter(support > -2 & support < 2) %>%
  ggplot() + geom_line(aes(support, dens * weight, group = factor(group))) + facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
  theme(legend.position = 'none') + scale_x_continuous(limits = c(NA, 2)) + theme_bw() + 
  labs(x = NULL, y = NULL)

# Predictive densities for each cluster

predDens <- data.frame()
support <- seq(-6, 6, length.out = 1000)
for(i in 1:500){
  pdens <- rep(0, 1000)
  weight <- sum(sample$k == i) / length(sample$k)
  
  # Evaluate via the 100 samples from q
  for(s in 1:100){
    th <- sample$theta[, i, s]
    pdens <- pdens + dnorm(support, th[2], exp(th[1])) / 100
  }
  predDens <- rbind(
    predDens,
    data.frame(support = support,
               dens = pdens,
               group = i,
               weight = weight))
}

# Get the weights of each cluster
predDens %>%
  group_by(group) %>%
  filter(support == min(support)) %>% 
  arrange(desc(weight)) %>%
  .$weight -> weights

# Create the second panel, integrate the density to get means and variance
predDens %>%
  filter(weight >= weights[50]) %>%
  group_by(group, weight) %>%
  summarise(mean = sum(dens * support * (support - lag(support)), na.rm = TRUE),
            var = sum(dens * (support - mean)^2 * (support - lag(support)), na.rm = TRUE),
            sd = sqrt(var)) %>%
  ungroup() %>% 
  ggplot() + geom_point(aes(mean, sd, size = weight)) +
  labs(size = 'Probability', x = 'Predictive Mean', y = 'Predictive Standard Deviation') + theme_bw()



# Calculate the individual densities for the third panel
indivs <- DPMindivDens(y[1:batches[t], ], UVBfit$lambda, S = 50, support = support)

# Choose 50 verhicles randomly
sub <- sample(N, 50)

# Plot the third panel
indivs %>%
  as.data.frame() %>%
  mutate(support = support) %>%
  gather(car, dens, -support) %>%
  mutate(id = as.numeric(substr(car, 2, 10))) %>%
  filter(id %in% sub) %>%
  ggplot() + geom_line(aes(support, dens, group = id), alpha = 0.20) + 
  geom_line(data = predDens %>%
              group_by(support) %>%
              summarise(dens = sum(dens * weight)),
            aes(support, dens), size = 2, colour = '#202050') + 
  labs(x = 'Lateral Lane Position', y = 'Predictive Density') + 
  theme_bw()


# Fit MFVB for each time period

cavilist <- list()
for(t in 2:length(batches)){
  cavilist[[t-1]] <- CAVI(y[1:batches[t],], 0, 10, 1, ncol(y), 100)
}


### Various forecasting methods

# Add each forecast to the list to avoid memnory allocation issues
fcList <- list()
# Count which element of the list we are currently at
counter <- 1
# Forecast all 500 vehicles
fcN <- 500
time <- batches[2:length(batches)]
# Loop though time periods
for(t in 1:6){
  print(t)
  Tn <- time[t]
  # The individual t predicitve density
  for(i in 1:fcN){
    
    muhat <- mean(data[1:Tn, i])
    varhat <- var(data[1:Tn, i]) / Tn * (Tn + 1)
    fcDens <- dt((data[Tn + 1:50, i] - muhat) / sqrt(varhat), Tn-1, log = TRUE) - log(sqrt(varhat))
    fcList[[counter]] <- tibble(t = Tn + 1:50,
                       id = i,
                       ls = fcDens,
                       method = 'Independent',
                       maxT = Tn)
    counter <- counter + 1
  }
  
  # UVB Forecasts
  
  fcDens <- DPMindivDens(y[1:Tn, 1:N],
                         lamb[[t]])
  

  # Add each forecast to the list
  for(i in 1:fcN){
    sub <- tibble(support = seq(-6, 6, length.out = 1000),
                  dens = fcDens[,i])
    
    fun <- approxfun(sub$support, sub$dens)
    fcDens <- log(fun(data[Tn + 1:50, i]))

    fcList[[counter]] <- tibble(t = Tn + 1:50,
                                id = i,
                                ls = fcDens,
                                method = 'UVB',
                                maxT = Tn)
    counter <- counter + 1
    
  }
  
  # Repeat for SVB
  
  fcDens <- DPMindivDens(y[1:Tn, 1:N],
                         lambSVB[[t]])
  
  for(i in 1:fcN){
    sub <- tibble(support = seq(-6, 6, length.out = 1000),
                  dens = fcDens[,i])
    
    fun <- approxfun(sub$support, sub$dens)
    fcDens <- log(fun(data[Tn + 1:50, i]))
    
    fcList[[counter]] <- tibble(t = Tn + 1:50,
                                id = i,
                                ls = fcDens,
                                method = 'SVB',
                                maxT = Tn)
    counter <- counter + 1
    
  }
  
 
  # Repeat for MFVB
  
  
  cavifit <- cavilist[[t]]
  
  # We will get a density for each of the N clusters
  predDensCAVIT <- matrix(0, 1000, N)
  for(i in 1:500){
    mean <- rnorm(100, cavifit$muHat[i], sqrt(cavifit$lambdaHat[i]))
    var <- invgamma::rinvgamma(100, cavifit$alphaHat[i], cavifit$betaHat[i])
    for(m in 1:100){
      predDensCAVIT[,i] <- predDensCAVIT[,i] + dnorm(support, mean[m], sqrt(var[m])) / 100
    }
  }

  # Then we will weight the clusters for each vehicle according to the estimated value of rho
  carDensC <- NULL
  for(i in 1:fcN){
    probsT <- cavifit$rho[i,]
    predDensCAVIT %>%
      as.tibble() %>%
      mutate(support = support) %>%
      gather(group, dens, -support) %>%
      mutate(p = rep(probsT, each = 1000)) %>%
      group_by(support) %>%
      summarise(dens = sum(dens * p)) %>%
      mutate(vehicle = i) -> df
    carDensC <- rbind(carDensC, df)
  }
  
  
  for(i in 1:fcN){
    sub <- filter(carDensC, vehicle == i)
    fun <- approxfun(sub$support, sub$dens)
    fcDens <- log(fun(data[Tn+1:50, i]))

    fcList[[counter]] <- tibble(t = Tn+1:50,
                                id = i,
                                ls = fcDens,
                                method = 'MFVB',
                                maxT = Tn)
    counter <- counter + 1
  }
}

cbPalette <- c("#009E73", "#D55E00", "#E69F00", "#56B4E9",  "#CC79A7")

# Create the plot

fcList %>%
  bind_rows() %>%
  filter(is.finite(ls)) %>%
  mutate(method = factor(method, levels = c('Independent', 'MFVB', 'SVB',  'UVB')),
         t = t - maxT,
         maxT = paste0('T', (maxT - 25)/25, ' = ', maxT)) %>%
  group_by(id, method, maxT) %>%
  mutate(ls = cumsum(ls)) %>%
  ungroup() %>%
  group_by(method, t, maxT) %>%
  summarise(ls = mean(ls, na.rm = TRUE)) %>%
  ggplot() + geom_line(aes(t, ls, colour = method), size = 0.6) + 
  theme_bw() + 
  theme(axis.line = element_blank(), legend.position = 'bottom') +
  facet_wrap(~maxT)  +
  labs(x = expression('Tn + h'), y = 'Mean Cumulative Log Score (MCLS)', colour = 'Method') + 
  scale_colour_manual(values=cbPalette)
  

