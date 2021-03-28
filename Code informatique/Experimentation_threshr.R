rm(list=ls())
library(threshr)

# Set a vector of training thresholds
data("danish")
u_vec <- quantile(danish, probs = seq(0.1, 0.95, by = 0.01))
thresh_cv <- ithresh(data = danish, u_vec = u_vec, trans='none')
summary(thresh_cv)

plot(thresh_cv, lwd = 2, cex.axis = 0.8)
plot(thresh_cv, which_u='best')


x <- evir::rgpd(2e+3, xi=0.1, mu=0, beta=2)
u_vec <- quantile(x, probs = seq(0, 0.5, by = 0.01))
thresh_cv <- ithresh(data = x, u_vec = u_vec, trans='BC', n=1e+4)
summary(thresh_cv)[3]

plot(thresh_cv)
plot(thresh_cv, which_u='best')


# Note:
# 1. Smoother plots result from making n larger than the default n = 1000.
# 2. In some examples below validation thresholds rather higher than is
#    advisable have been used, with far fewer excesses than the minimum of
#    50 suggested by Jonathan and Ewans (2013).

## North Sea significant wave heights, default prior -----------------------
#' # A plot akin to the top left of Figure 7 in Northrop et al. (2017)
#' # ... but with fewer training thresholds

u_vec_ns <- quantile(ns, probs = seq(0.1, 0.9, by = 0.1))
ns_cv <- ithresh(data = ns, u_vec = u_vec_ns, n_v = 2)
plot(ns_cv, lwd = 2, add_legend = TRUE, legend_pos = "topright")
mtext("significant wave height / m", side = 3, line = 2.5)

## Gulf of Mexico significant wave heights, default prior ------------------

u_vec_gom <- quantile(gom, probs = seq(0.2, 0.9, by = 0.1))
# Setting a prior using its name and parameter value(s) --------------------
# This example gives the same prior as the default
gom_cv <- ithresh(data = gom, u_vec = u_vec_gom, n_v = 2, prior = "mdi",
                  h_prior = list(a = 0.6))

## Setting a user-defined (log-)prior R function ---------------------------
# This example also gives the same prior as the default
# (It will take longer to run than the example above because ithresh detects
#  that the prior is an R function and sets use_rcpp to FALSE.)

user_prior <- function(pars, a, min_xi = -1) {
   if (pars[1] <= 0 | pars[2] < min_xi) {
      return(-Inf)
   }
   return(-log(pars[1]) - a * pars[2])
}
user_bin_prior <- function(p, ab) {
   return(stats::dbeta(p, shape1 = ab[1], shape2 = ab[2], log = TRUE))
}
gom_cv <- ithresh(data = gom, u_vec = u_vec_gom, n_v = 2, prior = user_prior,
                  h_prior = list(a = 0.6), bin_prior = user_bin_prior,
                  h_bin_prior = list(ab = c(1 / 2, 1 / 2)))

