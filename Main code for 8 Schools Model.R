library(rstan)
library(tidyverse)
source('Functions/gaussianDensity.R')
source('Functions/8schools funs.R')
Rcpp::sourceCpp('Functions/8schools.cpp')

# Main file containing the 8 Schools Analysis


# Data 
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

# Fit MCMC

fit <- stan(file = 'Functions/8schools.stan', data = schools_dat, 
            iter = 20000, chains = 1)

# Extract Samples

draws <- data.frame(mu = fit@sim$samples[[1]]$mu,
                    tau = fit@sim$samples[[1]]$tau,
                    t1 = fit@sim$samples[[1]]$`theta[1]`,
                    t2 = fit@sim$samples[[1]]$`theta[2]`,
                    t3 = fit@sim$samples[[1]]$`theta[3]`,
                    t4 = fit@sim$samples[[1]]$`theta[4]`,
                    t5 = fit@sim$samples[[1]]$`theta[5]`,
                    t6 = fit@sim$samples[[1]]$`theta[6]`,
                    t7 = fit@sim$samples[[1]]$`theta[7]`,
                    t8 = fit@sim$samples[[1]]$`theta[8]`)

# Save in new format

draws %>% 
  rename(`theta[1]` = t1, `theta[2]` = t2, `theta[3]` = t3, `theta[4]` = t4,
         `theta[5]` = t5, `theta[6]` = t6, `theta[7]` = t7, `theta[8]` = t8) %>%
  cbind(iter = 1:20000) %>%
  filter(iter > 10000) %>%
  gather(var, draw, -iter) -> MCMCdens

# A few plots

MCMCdens %>%
  filter(iter %% 20 == 0) %>%
  ggplot() + geom_line(aes(iter, draw)) + 
  facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
  theme_bw() + labs(x = NULL, y = NULL)

ggplot(MCMCdens) + geom_density(aes(draw)) + 
  facet_wrap(~var, scales = 'free', labeller = label_parsed) + 
  theme_bw() + labs(x = NULL, y = NULL)

# Calculate density function for each margin

MCMCdens %>%
  split(.$var) %>%
  map(calcDens) -> MCMClist

MCMClist <- lapply(MCMClist, function(x) data.frame(support = x$x, density = x$y))

# Result is the marginal posterior pdf at each value of support and variable: mu, tau, thetas
MCMCdens2 <- bind_rows(MCMClist)
MCMCdens2$var <- rep(c('mu', 'tau', paste0('theta[', 1:8, ']')), rep(1000, 10))
MCMCdens2$Inference <- 'MCMC'

# Extract the support used for each marginal value so we can calculate the VB margins on the same grid
mcmcsupport <- matrix(MCMCdens2$support[c(1001:2000, 1:1000, 2001:10000)], ncol = 10)



# Reformat data as a matrix to input to c++ functions used for VB

schoolsMat <- matrix(0, 8, 2)
schoolsMat[,1] <- schools_dat$y
schoolsMat[,2] <- schools_dat$sigma

# Starting values for lambda: mean is the zero vector, upper tri of variance is 10 * Identity Matrix
lambda <- c(rep(0, 10), diag(10, 10))

# Get all possible school order permutations

permutations <- gtools::permutations(8, 8)
nPermute <- 100

set.seed(4)

# Sample from possible permutations
selections <- sample(factorial(8), nPermute)
KLresults <- list()
densities <- list()

# Start with 1 schools, add 1 school per update
startJ <- 1
updateSize <- 1
nUpdates <- ceiling((8 - startJ) / updateSize)

# Fit VB to all schools (ie SVB)

VB <- reparamVB(schoolsMat, lambda, schoolsInit, J = 8, S = 50, dimTheta = 10, threshold = 0.005, alpha = 0.5)

# Save mean and variance matrix of the VB approximation
# The first 10 elemtns of lambda are the mean
# The next 100 is the upper triangle of the variance matrix, so convert that to Sigma
VBlist <- list(mean = VB$lambda[1:10],
               sigma = t(matrix(VB$lambda[11:110], 10)) %*% matrix(VB$lambda[11:110]))


# This function inputs the mean and standard deviation and calculate the density of the given marginal approximations
# So the first variable, tau, has an exponential transform (we fit an approximation to log tau!), and the rest 
# have an identity transform. 
# The grid is calcualted on the same support from the earlier MCMC margins
VBdens <- gaussianDensity(mu = VBlist$mean,
                                  sd = sqrt(diag(VBlist$sigma)),
                                  transform = c('exp', rep('identity', 9)),
                                  names = c('tau', 'mu', paste0('theta[', 1:8, ']')),
                                  support = mcmcsupport)

# Calcuate the KL divergence for each margin
VBdens$Inference <- 'SVB'
KLV <- KLdiv(VBdens, MCMCdens2)

# We will now calculate the vine copula used to estimate the joint pdf of the MCMC distribution
# THe pair of functions: margins and vine will allow us to estimate the pdf of any given value under the posterior
MCMCdens %>%
  spread(var, draw) %>%
  select(-iter) -> mcmcDraws

margins <- apply(mcmcDraws, 2, function(x){
  dens <- density(x)
  approxfun(dens)
}) 

cdf <- apply(mcmcDraws, 2, function(x){
  y <- seq(1/length(x), 1, 1/length(x))
  y[rank(x)]
})

vine <- VineCopula::RVineStructureSelect(cdf[seq(1, 90000, 200), ])

# The next step for the KL divergence is to sample from q as vbDraws,
# and calculate the density of these draws (under q) as vbDens
# We will sample 500 times 

vbDraws <- matrix(0, 500, 10)
vbDens <- rep(0, 500)

vbDraws <- mvtnorm::rmvnorm(500, VBlist$mean, VBlist$sigma)
vbDens <- mvtnorm::dmvnorm(vbDraws, VBlist$mean, VBlist$sigma, log = TRUE)

# MCMC results are mu then tau, vb is saved as log tau then mu. Switch the order of the first two columns, 
# Then take tau = exp(log(tau))
# Now VB columns will match MCMC columns
temp <- vbDraws[,1]
vbDraws[,1] <- vbDraws[,2]
vbDraws[,2] <- exp(temp)
vbDens <- vbDens - temp

# Calculate the density (under p) of these draws
vbMargins <- apply(vapply(1:10, function(x) margins[[x]](vbDraws[,x]), runif(500)), 1, prod)

vbcdf <- apply(vbDraws, 2, function(x){
  y <- seq(1/length(x), 1, 1/length(x))
  y[rank(x)]
})

vbCop <- VineCopula::RVinePDF(vbcdf, vine)
postDens <- vbCop * vbMargins

# We now have: 
# 1) Samples of theta, mu, tau ~ q
# 2) The density of these under q
# 3) The density of these under p
# Put them together to get the KL divergence
klVB <- vbDens - log(postDens)
klVB <- mean(klVB[is.finite(klVB)])
KLV <- rbind(KLV, tibble(var = 'joint', kl = klVB))
KLV

# Now we can run UVB / UVB-IS nPermute many times
for(j in 1:nPermute){
  # Print progress
  print(paste(j, Sys.time()))
  # Grab the j'th permutation of the data
  schoolsPermute <- schoolsMat[permutations[selections[j], ], ]
  
  # Initial fit of UVB (equiv to SVB with 1 school)
  J <- startJ
  UVBInit <- reparamVB(data = schoolsPermute,
                   lambda = lambda, 
                   model = schoolsInit,
                   J = J,
                   S = 50,
                   dimTheta = 2 + J,
                   threshold = 0.005,
                   alpha = 0.75)
  

  UVBlambda <- UVBInit$lambda
  
  # Update UVB nUpdates times
  
  for(i in 1:nUpdates){
    priorMean <- UVBlambda[1:(2 + J)]
    priorU <- matrix(UVBlambda[11:110], 10)
    priorLinv <- solve(t(priorU))
    
    # Update J
    
    J <- min(8, J + updateSize)
    
    UVBUpdate <- reparamVB(data = schoolsPermute,
                       lambda = UVBlambda, 
                       model = schoolsUVB,
                       J = J, 
                       S = 50, 
                       dimTheta = J + 2, 
                       threshold = 0.005,
                       alpha = 0.5,
                       mean = priorMean,
                       Linv = priorLinv)
    UVBlambda <- UVBUpdate$lambda
    
  }
  UVBlist <- list(mean = UVBlambda[1:10],
                  U = c(UVBlambda[11:110]))
  
  # Do this again for UVB-IS, starting at the UVBInit values
  ISsamples <- 500
  J <- startJ
  
  
  ISlambda <- UVBInit$lambda
  
  
  for(i in 1:nUpdates){
    # The number of schools in the previous update
    oldJ <- J
    # The number of schools after the current update
    J <- min(8, J + updateSize)
    
    # Grabbing the parameters of the previous approximation, slightly different 
    
    if(i == 1){
      priorMean <- ISlambda[1:(2 + J)]
      priorU <- matrix(ISlambda[11:110], 10)
      priorSigma <- (t(priorU) %*% priorU)[1:(2+J), 1:(2+J)]
      ISlambda <- c(priorMean, log(sqrt(diag(priorSigma))))
    } else {
      priorMean <- c(ISlambda[1:(2 + oldJ)], rep(0, updateSize))
      priorSd <- c(exp(ISlambda[(oldJ+3):(2*oldJ+4)]), rep(2, updateSize))
      priorSigma <- diag(priorSd^2)
      ISlambda <- c(priorMean, log(priorSd))
    }
    
    # Sample from previous lambda
    ISdraws <- mvtnorm::rmvnorm(ISsamples, priorMean, priorSigma)
    # Calcualate the density of these draws under the old version, log-likelihood etc.
    qDens <- mvtnorm::dmvnorm(ISdraws, priorMean, priorSigma)
    qDensPrior <- mvtnorm::dmvnorm(ISdraws[,1:(2 + oldJ)], priorMean[1:(2+oldJ)], priorSigma[1:(2+oldJ), 1:(2+oldJ)])
    loglik <- rep(0, ISsamples)
    
    for(a in 1:ISsamples){
      loglik[a] <- sum(dnorm(ISdraws[a, 2+(oldJ+1):J], ISdraws[a, 2], exp(ISdraws[a, 1]), log = TRUE)) + 
        sum(dnorm(schoolsPermute[(oldJ+1):J,1], ISdraws[a, 3:(2+J)], schoolsPermute[(oldJ+1):J, 2], log = TRUE))
    }
    logjoint <- loglik + log(qDensPrior)
    logjoint <- logjoint - max(logjoint)
    
    # Perform update
    
    ISUpdate <- ISschools(lambda = ISlambda, 
                      qScore = schoolsScore,
                      samples = ISdraws,
                      dSamples = qDens,
                      logjoint = logjoint,
                      alpha = 0.5,
                      threshold = 0.001)
    ISlambda <- ISUpdate$lambda
  }
  
  
  ISlist <- list(mean = ISlambda[1:10],
                 U = exp(ISlambda[11:20]))
  
  # Undo the permuations 
  ord <- order(permutations[selections[j], ])
  
  # Calculate the density of the approximating dist. 
  
  UVBmu <- UVBlist$mean[c(1:2, ord+2)]
  UVBvar <- (t(matrix(UVBlist$U, 10)) %*% matrix(UVBlist$U, 10))[c(1:2, ord+2), c(1:2, ord+2)]
  UVBsd <- sqrt(diag(UVBvar))
  
  ISmu <- ISlist$mean[c(1:2, ord+2)]
  ISsd <- ISlist$U[c(1:2, ord+2)]
  
  UVBdens <- gaussianDensity(mu = UVBmu,
                                     sd = UVBsd,
                                     transform = c('exp', rep('identity', 9)),
                                     names = c('tau', 'mu', paste0('theta[', 1:8, ']')),
                                     support = mcmcsupport)
  
  ISdens <- gaussianDensity(mu = ISmu,
                                    sd = ISsd,
                                    transform = c('exp', rep('identity', 9)),
                                    names = c('tau', 'mu', paste0('theta[', 1:8, ']')),
                                    support = mcmcsupport)
  
  
  jointDens <- rbind(UVBdens, ISdens)
  jointDens$Inference <- rep(c('UVB', 'UVB-IS'), rep(10000, 2))
  jointDens$per <- selections[j]
  densities[[i]] <- jointDens
  
  # Marginal KL divergences
  
  KLU <- KLdiv(UVBdens, MCMCdens2)
  KLI <- KLdiv(ISdens, MCMCdens2)
  
  # Sample from UVB / UVB-IS, calculate joint KL divergence as we did earlier
  # The posterior has not changed, so we do not need to re-estimate the vine.
  
  uvbDraws <- mvtnorm::rmvnorm(500, UVBmu, UVBvar)
  uvbDens <- mvtnorm::dmvnorm(uvbDraws,  UVBmu, UVBvar, log = TRUE)
  
  temp <- uvbDraws[,1]
  uvbDraws[,1] <- uvbDraws[,2]
  uvbDraws[,2] <- exp(temp)
  uvbDens <- uvbDens - temp
  
  uvbMargins <- apply(vapply(1:10, function(x) margins[[x]](uvbDraws[,x]), runif(500)), 1, prod)
  
  uvbcdf <- apply(uvbDraws, 2, function(x){
    y <- seq(1/length(x), 1, 1/length(x))
    y[rank(x)]
  })
  
  uvbCop <- VineCopula::RVinePDF(uvbcdf, vine)
  postDens <- uvbCop * uvbMargins
  
  klUVB <- uvbDens - log(postDens)
  klUVB <- mean(klUVB[is.finite(klUVB)])
  
  isDraws <- mvtnorm::rmvnorm(500, ISmu, diag(ISsd^2))
  isDens <- mvtnorm::dmvnorm(isDraws,ISmu, diag(ISsd^2), log = TRUE)
  
  temp <- isDraws[,1]
  isDraws[,1] <- isDraws[,2]
  isDraws[,2] <- exp(temp)
  isDens <- isDens - temp
  
  isMargins <- apply(vapply(1:10, function(x) margins[[x]](isDraws[,x]), runif(500)), 1, prod)
  
  iscdf <- apply(isDraws, 2, function(x){
    y <- seq(1/length(x), 1, 1/length(x))
    y[rank(x)]
  })
  
  isCop <- VineCopula::RVinePDF(iscdf, vine)
  postDens <- isCop * isMargins
  
  klIS <- isDens - log(postDens)
  klIS <- mean(klIS[is.finite(klIS)])
  
  tibble(var = KLV$var,
         UVB = c(KLU$kl, klUVB),
         `UVB-IS` = c(KLI$kl, klIS),
         per = selections[j]) -> KLresults[[j]]
  
  
}

# Make some plots
bind_rows(densities) %>%
  group_by(support, var, Inference) %>%
  summarise(density = mean(density)) %>%
  ungroup() %>%
  rbind(MCMCdens2) %>%
  rbind(VBdens) %>%
  ggplot() + geom_line(aes(support, density, colour = Inference)) + 
  facet_wrap(~var, scales = 'free', labeller = label_parsed)

bind_rows(KLresults) %>%
  gather(method, kl, -var, -per) %>%
  group_by(per, method) %>% 
  filter(!any(is.na(kl)) & !any(is.infinite(kl)) & all(kl > 0)) %>% 
  ungroup() %>%
  group_by(var, method) %>%
  summarise(kl = mean(kl)) %>%
  spread(method, kl) %>%
  left_join(KLV) %>%
  rename(SVB = kl) %>%
  mutate(UVB = round(UVB - SVB, 2),
         `UVB-IS` = round(`UVB-IS` - SVB, 2)) %>%
  select(var, UVB, `UVB-IS`) %>% 
  t() %>%
  knitr::kable(format = 'latex')



cbPalette <- c(MCMC = "#000000", SVB = "#E69F00", UVB = "#56B4E9", `IS-UVB` = "#009E73")

  
ggplot() + geom_line(data = jointDens, aes(support, density, colour = Inference)) +
  facet_wrap(~var, scales = 'free', labeller = label_parsed, ncol = 5) + 
  theme_bw() + labs(x = NULL, y = NULL) + 
  scale_colour_manual(values=cbPalette, breaks = c('MCMC', 'SVB', 'UVB', 'IS-UVB')) + 
  theme(legend.position = 'bottom')




