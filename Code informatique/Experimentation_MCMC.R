library(MCMC4Extremes)
library(qrmtools)

setwd("C:/Users/alexl/Google Drive/1. Université/2. Stage hydrologie/Code informatique")
source('Rainfall_functions.R')

# Obtaining posterior distribution of a vector of simulated points
x=rGPD(100, shape=0.3, scale=1) # in this case beta is the scale parameter sigma

# Obtaning 1000 points of posterior distribution
ajuste = gpdp(x, 0, 500)
# plot.gpdp(ajuste, 'histogram')
# plot.gpdp(ajuste, 'predictive')

# bayesian method
params <- parametrise_GPD(x, method='bayesian')

ADGofTest::ad.test(x, pGPD,
                   shape=params[1],
                   scale=params[2])
ks.test(x, pGPD,
        shape=params[1],
        scale=params[2])

car::qqPlot(x, 'GPD', 
            shape=params[1],
            scale=params[2],
            lwd=0.5)

# PWM method
params <- parametrise_GPD(x, method='pwm')

ADGofTest::ad.test(x, pGPD,
                   shape=params[1],
                   scale=params[2])
ks.test(x, pGPD,
        shape=params[1],
        scale=params[2])

car::qqPlot(x, 'GPD', 
            shape=params[1],
            scale=params[2],
            lwd=0.5)

# ML method
params <- parametrise_GPD(x, method='ml')

ADGofTest::ad.test(x, pGPD,
                   shape=params[1],
                   scale=params[2])
ks.test(x, pGPD,
        shape=params[1],
        scale=params[2])

car::qqPlot(x, 'GPD', 
            shape=params[1],
            scale=params[2],
            lwd=0.5)
