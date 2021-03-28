rm(list=ls())
library(tidyverse)
library(readr)
library(lubridate)
library(Kendall)
library(copula)
library(VineCopula)

# Importation des fonctions
setwd("C:/Users/alexl/Google Drive/1. Université/2. Stage hydrologie/Stage-hydrologie-H2021/Code informatique")
source("Rainfall_functions.R")


# Importation des données ----
setwd("C:/Users/alexl/Google Drive/1. Université/2. Stage hydrologie/data")
data <- readxl::read_xlsx('Clearwater river.xlsx', sheet = 'T089R09W4')


# ================== Prétraitement des données =================================
# Informations sur la station
data$Township %>% unique()

# Formater les données
data %>% glimpse()

data <- data %>% transmute(
   DATE = as.Date(Date),
   PRCP = `Precip. (mm)`
) %>% filter(month(DATE) %in% c(4,5,6))
data %>% summary()

# Vérifier l'intégrité de la BD.
data_intigrity_assertions(data)


# Clustering
data <- preprocess_data(  # Appuyer la méthode de NA imputing avec la littérature
   data, 
   NA.neighbour_days = 2,
   C.separeting_days = 1
   )

# ================== Analyse préliminaire ======================================
# Vérification de la dépendance séquentielle ----
verify_serial_dependence(data)


# Vérification de la stationnarité ----
MannKendall(data$prcp_agg)
data %>% ggplot(aes(x=start_date, y=prcp_agg)) + 
   geom_path() +
   geom_smooth()

data %>% filter(time_since_last < 52) %>% select(time_since_last) %>% unlist() %>% MannKendall()
MannKendall(data$time_since_last)
data %>% ggplot(aes(x=start_date, y=time_since_last)) + 
   geom_path() +
   geom_smooth()


# ================== Séparation en données extrêmes et données normales.========
set.seed(42)
p_vec <- rnorm(100, mean=0.86, sd=0.04) %>% sort()
u_vec <- quantile(data$prcp_agg, probs = p_vec) 

u <- treshold_selection(data, u_vec, method='bader', gf_test = 'ad')

data.extremes <- calculate_excedences(data, u, calculate_W=T)
data.normales <- data %>%
   filter(prcp_agg < u) %>%
   select(prcp_agg, start_date) %>%
   group_by(year(start_date)) %>%
   summarise(sum(prcp_agg), .groups='drop') %>%
   rename('year' = `year(start_date)`, total_prcp = `sum(prcp_agg)`)

# ================== Distribution des précipitations saisonnières ==============
Z <- data.normales$total_prcp
mktest <- MannKendall(Z)$sl

hist(Z, freq=FALSE)
shapiro.test(Z)  # Le fait que les données soient non stationnaires peut expliquer la faible p-value.

Z_lm <- lm(total_prcp ~ year, data=data.normales)
Z_sd <- summary(Z_lm)$sigma

if (mktest > 0.05) {
   Z_lm$coefficients <- c(mean(Z), 0)
   Z_sd <- sd(Z)
   Z_mu <- mean(Z)
   Z_pred <- qnorm(ecdf(Z)(Z), Z_mu, Z_sd)
   
   car::qqPlot(Z, distribution="norm", mean=Z_mu, sd=Z_sd, lwd=0.3, id=F)
   abline(0,1, col='red', lwd=1)
} else {
   Z_mu <- predict(Z_lm)
   Z_pred <- qnorm(ecdf(Z)(Z), Z_mu, Z_sd)
   car::qqPlot(Z_lm, lwd=0.3, id=F, ylab='studentized residuals')
   abline(0,1, col='red', lwd=1)
}

ks.test(Z, Z_pred)
ADGofTest::ad.test(Z, pnorm, mean=Z_mu, sd=Z_sd)

# ================== Distribution des précipitations extrêmes ==================
MannKendall(data.extremes$excedence)
data.extremes %>% ggplot(aes(x=start_date, y=excedence)) + 
   geom_path() +
   geom_smooth()

params_GPD <- fit_GPD(data.extremes$excedence, method = 'mm')

ADGofTest::ad.test(data.extremes$excedence, pGPD,
                   shape=params_GPD[1],
                   scale=params_GPD[2])
ks.test(data.extremes$excedence, pGPD,
        shape=params_GPD[1],
        scale=params_GPD[2])
car::qqPlot(data.extremes$excedence, 'GPD', 
            shape=params_GPD[1],
            scale=params_GPD[2],
            lwd=0.5, id=F, cex=0.7,
            ylab='Empirical quantiles')
abline(0,1, col='red', lwd=1)


# ================== Distribution des temps inter-occurence (extrêmes) =========
data.extremes %>%
   ggplot(aes(x=start_date, y=W)) + 
   geom_path() +
   geom_smooth()

W_fitted <- fit_W(data.extremes, show_qqplot=T)
dw <- W_fitted$density
pw <- W_fitted$distribution
qw <- W_fitted$quantiles


# ================== Distribution de la duration (extrêmes) ====================
data.extremes %>%
   ggplot(aes(x=start_date, y=duration)) + 
   geom_path() +
   geom_smooth()

D_fitted <- fit_duration(data.extremes, show_qqplot=T)
dD <- D_fitted$density
pD <- D_fitted$distribution
qD <- D_fitted$quantiles


# ================== Analyse de la dépendance empirique ========================
# Dépendence entre les variables D et W.
cat('rhos de Spearman:', fill=T)
DW <- data.extremes %>%
   select(excedence, duration, W) %>%
   mutate(duration = as.numeric(duration)) 
calculate_correlation(DW, alternative = 'two-sided')

# Dépendence entre W et X
data.extremes %>% select(excedence, W) %>% na.omit() %>% 
   calculate_correlation(alternative = 'greater')

# Dépendence entre les normales et les valeurs extrêmes.
data.extremes %>% 
   mutate(year = year(start_date)) %>%
   group_by(year) %>%
   summarise(annual_extremes = sum(u + excedence), .groups='drop') %>% 
   right_join(data.normales, by='year') %>%
   select(-year) %>%
   cor(use = "pairwise.complete.obs", method='spearman')


# ============== Modélisation de la dépendence entre X-u|X>u et W ==============
# eval_dependence_trend(data.extremes, variable='W', stride=nrow(data.extremes)/2)
extreme_value_copulas <- c(4, 5, 104)
copulas_except_elliptic <- c(-1,-2)

best_copula.W <- copula_selection(
   data.extremes, variable='W', family_set = copulas_except_elliptic,
   gofTest = 'white')

# Copule empirique
data.extremes %>% 
   select(W, excedence) %>%
   ggplot(aes(x = pobs(W),
              y = pobs(excedence))) +
   geom_density_2d_filled(show.legend=F) + 
   xlab('Ranks of W') + 
   ylab('Ranks of (X-u|X>u)')

# Copule théorique
U_sim <- BiCopSim(1e+4, obj = best_copula.W)
colnames(U_sim) <- c('u1', 'u2')
U_sim %>% as_tibble() %>%
   ggplot(aes(x = u1, y = u2)) +
   geom_density_2d_filled(show.legend=F)


# ============== Modélisation de la dépendence entre X-u|X>u et D ==============
data.extremes %>% select(excedence, duration) %>% na.omit() %>% 
   calculate_correlation(alternative = 'greater')

# eval_dependence_trend(data.extremes, variable='duration',
# stride=nrow(data.extremes)/2)

best_copula.D <- copula_selection(
   data.extremes, variable='duration', family_set = copulas_except_elliptic, 
   gofTest = 'white')

# Copule empirique
data.extremes %>% 
   select(duration, excedence) %>%
   ggplot(aes(x = pobs(duration),
              y = pobs(excedence))) +
   geom_density_2d_filled(show.legend=F) + 
   xlab('Ranks of D') + 
   ylab('Ranks of (X-u|X>u)')

# Copule théorique
U_sim <- BiCopSim(1e+4, obj = best_copula.D)
colnames(U_sim) <- c('u1', 'u2')
U_sim %>% as_tibble() %>%
   ggplot(aes(x = u1, y = u2)) +
   geom_density_2d_filled(show.legend=F)


# ============== Modélisation de la dépendence entre X et Z ====================
data_XZ <- data.extremes %>% 
   mutate(year = year(start_date)) %>%
   group_by(year) %>%
   summarise(annual_extremes = sum(u + excedence), .groups='drop') %>% 
   right_join(data.normales) %>%
   rename(X = annual_extremes,
          Z = total_prcp)

# eval_dependence_trend(data.extremes, variable='duration',
#                       stride=nrow(data.extremes)/2)

best_copula.XZ <- copula_selection.XZ(
   data_XZ, family_set = copulas_except_elliptic, gofTest = 'white')


# ==== Simulation pour tenter de reproduire la distribution des pluies totales ====
years <- data$start_date %>% year() %>% unique()

simulations <- simulate_annual_prcp(years, W_fitted, D_fitted, params_GPD,
                                    best_copula.W, best_copula.D, best_copula.XZ,
                                    Z_lm, treshold=u, nsim=3e+2)


ps <- function(x, year) {
   ind <- which(years == year)
   ecdf(simulations[ind,])(x)
}
qs <- function(p, year) {
   ind <- which(years == year)
   n <- ncol(simulations)
   sort(simulations[ind,])[ceiling(p*n)]
}
TVaR <- function(p, year) {
   ind <- which(years == year)
   n <- ncol(simulations)
   sorted_obs <- sort(simulations[ind,])
   VaR <- sorted_obs[ceiling(p*n)]
   TVaR <- mean(sorted_obs[sorted_obs > VaR])
   return(TVaR)
}
qglobal <- function(p) {
   vec_sim <- as.vector(simulations)
   n <- length(vec_sim)
   q <- sort(vec_sim)[ceiling(p*n)]
   return(q)
}
pglobal <- function(x) {
   vec_sim <- as.vector(simulations)
   ecdf(vec_sim)(x)
}
dglobal <- function(x) {
   vec_sim <- as.vector(simulations)
   ecdf(vec_sim)(x + 0.5) - ecdf(vec_sim)(x - 0.5)
}


annual_prcp <- data %>% 
   group_by(year(start_date)) %>%
   summarise(prcp_agg = sum(prcp_agg), .groups='drop') %>%
   select(prcp_agg) %>%
   unlist(use.names = F)

n <- length(years)
qqplot(annual_prcp, qglobal((1:n)/(n+1)),
       xlab='Empirical quantiles',
       ylab='Theorical quantiles')
abline(0, 1, col='blue')

car::qqPlot(annual_prcp, "global",
            xlab='Empirical quantiles',
            ylab='Theorical quantiles',
            lwd=0.5, id=F)
abline(0, 1, col='red')

ks.test(annual_prcp, pglobal)
ad.test(annual_prcp, pglobal)
# ------------------------------------------------------------------------------
