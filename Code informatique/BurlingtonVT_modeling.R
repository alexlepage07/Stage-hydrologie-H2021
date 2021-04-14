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


# Importation des données de https://www.ncdc.noaa.gov/cdo-web/search ----
setwd("C:/Users/alexl/Google Drive/1. Université/2. Stage hydrologie/data/BurlingtonRainfall_Data")
data <- list.files(path=".", full.names = TRUE) %>% 
   lapply(read_csv) %>% 
   bind_rows
#'  PRCP = Precipitation (mm or inches as per user preference, inches to
#'  hundredths on Daily Form pdf file)

# ================== Prétraitement des données =================================
# Informations sur la station
data$STATION %>% unique()
data$NAME %>% unique()

BurlingtonVT <- data %>% filter(NAME == "BURLINGTON, VT US")
BurlingtonVT$STATION %>% unique()
BurlingtonVT$DATE %>% year() %>% unique() %>% max()

AirPort <- data %>% filter(NAME == "BURLINGTON INTERNATIONAL AIRPORT, VT US")
AirPort$STATION %>% unique()
AirPort$DATE %>% year() %>% unique()


# Faire le ménage dans les données
data1 <- data %>% filter(
   NAME == "BURLINGTON, VT US", year(DATE) <= 1942)
data2 <- data %>% filter(
   NAME == "BURLINGTON INTERNATIONAL AIRPORT, VT US", year(DATE) > 1942)
data <- rbind(data1, data2)

data <- data %>% select(c(DATE, PRCP)) # Retirer les colonnes non nécessaires.
data <- data %>% filter(month(DATE) %in% c(4,5,6)) # Ne conserver que les mois
# d'avril, de mai et de juin.


# Conversion des pouces en mm.
data <- data %>% mutate(PRCP = PRCP * 25.4)


# Vérifier l'intégrité de la BD.
# data_intigrity_assertions(data)
data %>% filter(year(DATE) == 1943) %>% select(DATE) %>% c()
# Ajout de données manquantes
data <- rbind(
   data,
   list("1943-06-01", 0.03),
   list("1943-06-02", 0),
   list("1943-06-03", 0)
) %>% arrange(DATE)


# Clustering
data <- preprocess_data(  # Appuyer la méthode de NA imputing avec la littérature
   data, 
   NA.neighbour_days = 2,
   C.separeting_days = 1
   )

# ================== Analyse préliminaire ======================================
# Vérification de la dépendance séquentielle ----
# verify_serial_dependence(data)


# Vérification de la stationnarité ----
MannKendall(data$prcp_agg)
data %>% ggplot(aes(x=start_date, y=prcp_agg)) + 
   geom_path() +
   geom_smooth() + ylab("K") + xlab('time')

MannKendall(data$time_since_last)
data %>% 
   filter(time_since_last < 52) %>%
   ggplot(aes(x=start_date, y=time_since_last)) + 
   geom_path() +
   geom_smooth() + ylab("time inter-clusters") + xlab('time')


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
hist(data.normales$total_prcp, freq=FALSE)
shapiro.test(data.normales$total_prcp)

Z_fitted <- fit_Z(data.normales)


# ================== Distribution des précipitations extrêmes ==================
MannKendall(data.extremes$excedence)
data.extremes %>% ggplot(aes(x=start_date, y=excedence)) + 
   geom_path() +
   geom_smooth()

params_GPD <- fit_GPD(data.extremes$excedence, method = 'mm')

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
   geom_smooth() +
   xlab('Time')

hist(data.extremes$W, freq=F,
     xlab='W', ylab='Probability', main=NULL)

W_fitted <- fit_W(data.extremes)

# Vérifier l'adéquation
x = data.extremes$W
qW = W_fitted$quantiles
.max = max(data.extremes$W)
t0 <- data.extremes$t0 + data.extremes$duration
nsim = 1e+4

simuls <- sapply(t0, function(.t) {
   qW(runif(nsim), .t)
}) %>% as.vector()
simuls <- simuls[simuls <= .max] %>% sort()

pSimul <- ecdf(simuls)
qSimul <- function(p) simuls[ceiling(p*length(simuls))]
dSimul <- function(x) {
   sapply(x, function(.x) {
      pSimul(.x+1) - pSimul(.x)
   })
}

car::qqPlot(x, 'Simul',
            xlab='Theorical quantiles',
            ylab='Empirical quantiles',
            lwd = 0.5, id = F)
abline(0, 1, col='red')


# ================== Distribution de la duration (extrêmes) ====================
data.extremes %>%
   ggplot(aes(x=start_date, y=duration)) + 
   geom_path() +
   geom_smooth() + xlab('Time') + ylab('D')
hist(as.numeric(data.extremes$duration), freq=F,
     xlab='Duration', ylab='Probability', main=NULL)

D_fitted <- fit_D(data.extremes)

# Vérifier l'adéquation
x = data.extremes$duration %>% as.numeric()
qD = D_fitted$quantiles
.max = max(x)
t0 = data.extremes$t0
nsim = 1e+4

simuls <- sapply(t0, function(.t) {
   qD(runif(nsim), .t)
}) %>% as.vector()
simuls <- simuls[simuls <= .max] %>% sort()

qSimul <- function(p) simuls[ceiling(p*length(simuls))]
qqplot(x, qSimul(ecdf(x)(x)),
            xlab='Theorical quantiles',
            ylab='Empirical quantiles')
abline(0, 1, col='red')


# ================== Analyse de la dépendance empirique ========================
# Dépendence entre les variables D et W.
cat('rhos de Spearman:', fill=T)
DW <- data.extremes %>%
   select(excedence, duration, W) %>%
   mutate(duration = as.numeric(duration))
cor(DW, method = 'spearman') %>% xtable::xtable(digits=4)
cor(DW, method = 'kendall') %>% xtable::xtable(digits=4)

DW %>% select(duration, W) %>% mutate_all(as.numeric) %>%
   calculate_correlation(alternative = 'greater')

uu <- data.extremes %>% select(excedence, W, duration) %>%
   mutate_all(as.numeric) %>% pobs()
scatterplot_Chaoubi(uu, kernel_dist='norm')


# Dépendence entre les normales et les valeurs extrêmes.
data.extremes %>% 
   mutate(year = year(start_date)) %>%
   group_by(year) %>%
   summarise(annual_extremes = sum(u + excedence), .groups='drop') %>% 
   right_join(data.normales, by='year') %>%
   select(-year) %>%
   cor(use = "pairwise.complete.obs", method='spearman')


# ============== Modélisation de la dépendence =================================
pX <- function(x) pGPD(x, params_GPD[1], params_GPD[2])
pW <- W_fitted$distribution
pD <- D_fitted$distribution

uu <- data.extremes %>% select(excedence, W, duration) %>%
   mutate_all(as.numeric) %>% as.matrix()
t0 <- data.extremes$t0
D <- data.extremes$duration
uu[,1] <- uu[,1] %>% pX()
uu[,2] <- uu[,2] %>% pW(t0+D)
uu[,3] <- uu[,3] %>% pD(t0)

cor(uu)
calculate_correlation(uu[,c(2,3)], alternative = 'greater')
scatterplot_Chaoubi(uu[which(uu[,1]>0 & uu[,2]>0),], kernel_dist='norm')

Matrix = c(1, 0, 0,
           2, 2, 0,
           3, 3, 3) %>% matrix(ncol=3, byrow=T)
RVineMatrixCheck(Matrix)

# BiCopCompare(uu[,1], uu[,2], familyset = -c(1,2))
familyset_12 <- c(5, 114, 13, 4, 0)

# BiCopCompare(uu[,1], uu[,3], familyset = -c(1,2))
familyset_13 <- c(10, 13, 4, 17, 5)

cop12 <- BiCopEstList(uu[,1], uu[,2], familyset = familyset_12, rotation=F)
cop13 <- BiCopEstList(uu[,1], uu[,3], familyset = familyset_13, rotation=F)

cop12$models
cop12$summary %>% xtable::xtable() %>% print(include.rownames=F)

cop13$models
cop13$summary %>% xtable::xtable() %>% print(include.rownames=F)

RVM <- RVineCopSelect(uu, familyset = c(-1,-2), indeptest = T,
                      Matrix = Matrix, selectioncrit ='BIC')


# ============== Modélisation de la dépendence entre X et Z ====================
data_VZ <- data.extremes %>% 
   mutate(year = year(start_date)) %>%
   group_by(year) %>%
   summarise(annual_extremes = sum(u + excedence), .groups='drop') %>% 
   right_join(data.normales) %>%
   rename(V = annual_extremes,
          Z = total_prcp) 
 
uu <- data_VZ %>% select(-year) %>% na.omit() %>% pobs()
# BiCopCompare(uu[,1], uu[,2], familyset = -c(1,2))
familyset = c(24, 5, 33, 224, 26)
cop_VZ <- BiCopEstList(uu[,1], uu[,2], family=familyset, rotations = F)
cop_VZ$model
cop_VZ$summary %>% xtable::xtable() %>% print(include.rownames=F)

best_copula.XZ <- BiCopSelect( 
      uu[, 1], uu[, 2],
      family = familyset,
      indeptest = T,
      selectioncrit = 'BIC'
   )
best_copula.XZ <- BiCop(family = 0)
# eval_dependence_trend(data.extremes, variable='duration',
#                       stride=nrow(data.extremes)/2)

# best_copula.XZ <- copula_selection.XZ(
#    data_XZ, family_set = c(-1,-2), gofTest = 'white')


# ==== Simulation pour tenter de reproduire la distribution des pluies totales ====
years <- data$start_date %>% year() %>% unique()

simulations <- simulate_annual_prcp(years, Z_fitted, W_fitted, D_fitted, params_GPD,
                                    RVM=RVM, best_copula.XZ=best_copula.XZ, 
                                    treshold=u, nsim=1e+3)

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
   ecdf(vec_sim)(x) - ecdf(vec_sim)(x - 1)
}


annual_prcp <- data %>% 
   group_by(year(start_date)) %>%
   summarise(prcp_agg = sum(prcp_agg), .groups='drop') %>%
   select(prcp_agg) %>%
   unlist(use.names = F)

car::qqPlot(annual_prcp, "global",
            xlab='Theorical quantiles',
            ylab='Empirical quantiles',
            lwd=0.5, id=F, line='robust')
abline(0, 1, col='red')

ks.test(annual_prcp, pglobal)
ad.test(annual_prcp, pglobal)


# Pour aller plus loin ====
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
# ------------------------------------------------------------------------------
