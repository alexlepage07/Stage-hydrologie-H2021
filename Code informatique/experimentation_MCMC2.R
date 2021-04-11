library(MHadaptive) # Permet d'utiliser la méthode du maximum de vraisemblance généralisée
library(qrmtools)   # Donne accès aux fonctions dGPD, pGPD, qGPD et rGPD.
library(tidyverse)  # Package pour être efficace en data science.

nll_Pareto <- function(param, data, t0, trend) {
   shape <- param[1]
   if (trend) {
      scale <- function(.t) {
         return(exp(param[2] + param[3] * as.numeric(.t)))
      }
   } else {
      scale <- function(.t) param[2]
   }
   
   if (min(shape, scale(t0)) <= 0) return(0)
   nll <- -sum(dPar(data, shape, scale(t0), log=T))
   return(nll)
} 

fit_W_GPD <- function(init_par, data, t0, trend, quiet=T, seed=2021) {
   
   prior_reg <- function(param, trend) {
      
      prior_shape <- dgamma(param[1], 5, 1, log = T)
      
      if (trend) {
         prior_a <- dnorm(param[2], init_par[2], 2, log=T)
         prior_b <- dnorm(param[3], 0, 1, log=T)
         prior_scale <- prior_a + prior_b
      } else {
         prior_scale <- dgamma(param[2], 15, 1/5, log=T)
      }
      return(prior_shape + prior_scale)
   }
   
   li_func <- function(param, data, t0, trend) {
      return(prior_reg(param, trend) - nll_Pareto(param, data, t0, trend))
   }
   
   if (trend) {
      par_names <- c('shape', 'a', 'b')
   } else {
      par_names <- c('shape', 'scale')
   }
   set.seed(seed)
   mcmc_r <- Metro_Hastings(
      li_func = li_func,
      pars = init_par,
      par_names = par_names,
      iterations = 2e+4,
      burn_in = 1e+4,
      quiet=T,
      data=data,
      t0=t0,
      trend=trend
   )
   init_par <- apply(mcmc_r$trace, 2, getmode)
   
   mcmc_r <- Metro_Hastings(
      li_func = li_func,
      pars = init_par,
      prop_sigma = mcmc_r$prop_sigma,
      par_names = par_names,
      iterations = 2e+4,
      burn_in = 1e+4,
      quiet=quiet,
      data=data,
      t0=t0,
      trend=trend
   )
   par <- apply(mcmc_r$trace, 2, getmode)
   names(par) <- par_names
   return(par)
}

# Avec stationnarité
.t0 <- as.numeric(t0)
.t0 <- (.t0 - min(.t0)) / (max(.t0) - min(.t0))

data_test <- sapply(.t0, function(.t){
   rPar(1, 2, exp(4 + 1e-2*as.numeric(.t)))})
fit_W_GPD(init_par = c(2, 5, 0), data=data_test, t0=.t0, trend=T)

# Sans stationnarité
data_test <- rPar(1000, 2, 50)
fit_W_GPD(init_par = c(10, 40), trend=F)

