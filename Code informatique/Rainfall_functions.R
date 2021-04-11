library(ADGofTest)  # Statistique d'Anderson-Darling
library(copula)     # Librairie générale pour travailler avec des copules.
library(eva)        # Un autre papier plus récent sur la détermination du seuil.
library(lubridate)  # Travailler efficacement avec les dates
library(modifiedmk) # Version modifiée du test de Mann-Kenndall pour des données auto-corrélées.
library(parallel)   # Permet de parréliser les simulations du modèle final.
library(qrmtools)   # Donne accès aux fonctions dGPD, pGPD, qGPD et rGPD.
library(threshr)    # Détermination du seuil de valeur extrême de façon automatique
library(tidyverse)  # Package pour être efficace en data science.
library(utils)      # Permet de faire des bars de progression.
library(VineCopula) # Librairie efficace pour faire la sélection de copule.

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}


data_intigrity_assertions <- function(data) {
   #' Fonction qui vérifie que chacune des années de la base de données
   #' contient exactement 91 jours. Dans le cas contraire, un message d'erreur 
   #' est retourné.
   #' 
   #' @param data (tibble): DataFrame tibble contenant les données sur lesquels 
   #'        les opérations seront appliquées. Les champs de cette base de 
   #'        données doivent au moins contenir:
   #'        - DATE (date): La date, en format as.Date, couvrant la période du 
   #'               1er avril au 30 juin de chaucune des années disponibles 
   #'               (91 jours par année).
   #' 
   years <- data$DATE %>% year() %>% unique() %>% sort()
   
   # Regarder s'il y a des années manquantes
   if (diff(years) %>% unique() %>% sum() != 1) {
      missing_years = years[which(diff(years) != 1) - 1] + 1
      if (length(missing_years) == 1)
         stop(paste("The year", missing_years, "is missing"))
      else
         stop(paste("The years", list(missing_years), "are missing"))
   }
   
   # Regarder si toutes les années contiennent 91 jours.
   if ((max(years) - min(years) + 1) * 91 != nrow(data)) {
      nb_of_days <- data$DATE %>% year() %>% table()
      print(nb_of_days[which(nb_of_days !=91 )])
      stop(paste("Some years doesn't include 91 days."))
   }
}


preprocess_data <- function(data, NA.neighbour_days=2, C.separeting_days=1, n_print=0) {
   #' Fonction qui impute les données manquantes et qui fait le regroupement des
   #' jours de pluie en périodes de précipitations continues.
   #' 
   #' @param data (tibble): DataFrame tibble contenant les données sur lesquels 
   #'        les opérations seront appliquées. Les champs de cette base de 
   #'        données doivent être:
   #'        - DATE (date): La date, en format as.Date, couvrant la période du 
   #'               1er avril au 30 juin de chaucune des années disponibles 
   #'               (91 jours par année).
   #'        - PRCP (float): La quantité de pluie tombée pour une journée
   #'        La colonne DATE 
   #' @param NA.neighbour_days (int): Nombre de jours à regarder avant et après
   #'        les jours où il y a des précipitations avec des NA afin de réaliser 
   #'        une inférence par moyenne.
   #' @param C.separeting_days (int): Nombre de jours d'ensoleillement qui 
   #'        séparent deux périodes de précipitations (clusters).
   #' @param n_print (int): Nombre d'observations à afficher pour valider que
   #'        le clustering et le calcul des temps inter-occurence se sont
   #'        bien fait. Si n_print=0, alors aucun tableau n'est affiché.
   #'       
   #' @return la base de données data dont les observations ont été regroupées.
   #'         de nouvelles colonnes sont insérées:
   #'         - cluster: l'indice de groupement utilisé;
   #'         - prcp_agg: les précipitation aggrégées pour une période donnée;
   #'         - duration: la durée de la période de pluie;
   #'         - time_since_last: le temps écoulé depuis la dernière période de 
   #'                    pluie. Si la période correspond à la première d'une 
   #'                    année, alors cette variable prend la valeur NA.
   #'         - start_date: La date de début de la période de pluie.
   
   na_imput <- function(data, neighbours_days=2) {
      #' Fonction qui fait de l'inférence sur les valeurs manquntes en considérant
      #' en utilisant la moyenne des jours avant et après.
      #' 
      #' @param data : Les données sur lesquelles imputer les données manquantes.
      #' @param neighbours_days : Nombre de jours avant et après la période avec 
      #'             des NA qu'il faut utiliser pour calculer une moyenne.
      #'             
      #' @return Les données avec les valeurs manquantes modifiées.
      na_mask <- is.na(data$PRCP)
      if (sum(na_mask) == 0) {
         return(data)
      }
      na_dates <- data[na_mask, 'DATE']
      na_dates_bornes <- cbind(
         na_dates - neighbours_days,
         na_dates + neighbours_days
      )
      colnames(na_dates_bornes) <- c('date_inf', 'date_sup')
      
      for (r in 1:nrow(na_dates)) {
         dates_to_check = seq(
            from = na_dates_bornes[r,1],
            to = na_dates_bornes[r,2], by=1)
         
         prcp <- data %>% 
            filter(DATE %in% dates_to_check) %>%
            select(PRCP) %>% unlist()
         
         data[data$DATE %in% na_dates[r,1], 'PRCP'] <- mean(prcp, na.rm = TRUE)
      }
      
      return(data)
   }
   
   rain_clustering <- function(data, nb_sunny_days=1, n_print=20) {
      #' Fonction qui regroupe les données journalières par périodes de 
      #' précipitations.
      #' 
      #' @param data : Une base de données tibble.
      #' @param nb_sunny_days : Nb de jours nécessaires pour séparer deux 
      #'        périodes (Par défaut=1).
      #'
      #'@return un dataframe tibble où les précipitations continues sont agrégées.
      years_ <- data$DATE %>% year() %>% unique()
      n <- nrow(data)
      cluster_indicator <- numeric(n)
      time_since_last <- numeric(n)
      C = 0
      W = 1
      
      for (i in 2:n) {
         
         if (data[i, 'PRCP'] == 0) {
            # Si on a un jour d'ensoleillement
            cluster_indicator[i]  <- 0
            W = W + 1
         } else {
            # Si on a un jour de pluie on lui attribue un cluster
            last_days <- max(i - nb_sunny_days, 1):(i-1)
            
            if (sum(cluster_indicator[last_days]) == 0) {
               #' Si les @nb_sunny_days derniers jours ont été ensoleillés,
               #' alors on change de cluster.
               C = C + 1
               time_since_last[i] <- W
               W <- 0
            }
            cluster_indicator[i] <- C
         }
         
         # Deux années distinctes ne peuvent appartenir à un même cluster.
         if (i %% 91 == 0)  {
            C = C + 1
            time_since_last[i] <- 0
            W <- 0
         }
      }
      
      if (n_print > 0) {
         print(rbind(
            "Cluster_indicator" = head(cluster_indicator, n_print),
            "Precipitations" = head(ceiling(data$PRCP), n_print),
            "Time_since_last" = head(ceiling(time_since_last), n_print)
         ), quotes=F)
      }
      
      data <- data %>% 
         mutate(cluster = cluster_indicator,
                time_since_last = time_since_last) %>% 
         filter(cluster > 0) %>% 
         group_by(cluster) %>% 
         summarise(
            prcp_agg = sum(PRCP),
            start_date = min(DATE),
            duration = as.numeric(max(DATE) - min(DATE) + 1),
            time_since_last = sum(time_since_last),
            .groups = 'drop'
         ) %>% select(-cluster)
      
      return(data)
   }
   
   # Preprocessing
   data <- data %>% select(c(DATE, PRCP)) # Retirer les colonnes non nécessaires.
   data <- data %>% filter(month(DATE) %in% c(4,5,6)) # Ne conserver que les mois
   
   # Vérifier que les données contiennent le bon nombre de jours par années et 
   # qu'il ne manque pas d'années.
   data_intigrity_assertions(data)
   
   # Missing imputation
   data <- na_imput(data, NA.neighbour_days)

   # clustering
   data <- rain_clustering(data, C.separeting_days, n_print)
   
   return(data)
}


fit_GPD <- function(x, method='bayesian') {
   #' Fonction qui paramétrise la loi GPD telle que
   #' 
   #' F(x) = 1 - (1 + \xi x/\beta)^{-1 / \xi} if \xi != 0
   #' F(x) = 1 - exp(-x / \beta) if \xi == 0.
   #' 
   #' @param x (vector): Vecteur d'observations i.i.d. supposées suivre une 
   #'        loi GPD.
   #' @param method (str): La méthode utilisée pour définir les paramètres de 
   #'        la loi. Correspond à l'une des options suivantes:
   #'        - 'pwm' : Probability Weighted Moments (Voir Hosking et al. 1985)
   #'                  Particulièrement efficace si le paramètre de forme \xi 
   #'                  est compris dans l'intervalle [-0.4, -0.2].
   #'        - 'mle' : Maximum Likelihood. Utile lorsque \xi appartient à 
   #'                  [0.2, 0.4].
   #'        - 'mm' : Moments Method. Efficace lorsque \xi est compris dans 
   #'                 [-0.2, 0.2]. Remarque: Divergera si \xi > 0 et \beta <= 2!
   #'        - 'vwl : Variance Weibull Least Squares (Makkonen & Tikanmäki 2019)       
   #'                 Méthode minimisant l'écart quadratique entre les quantiles
   #'                 empiriques et les quantiles théoriques en utilisant du 
   #'                 bootstrap.
   #'        - 'bayesian : voir Castellanos, M. A. and Cabras, S. (2007)
   #' @return les paramètres de forme \xi et d'échelle \beta de la famille de
   #'         lois GPD.
   
   if (method == 'pwm') {
      param <- fit_GPD_PWM(x)
      xi <- param[1]
      scale <- param[2]
   }
   else if (method == 'ml') {
      param <- fit_GPD_MLE(x)
      xi <- param[1]
      scale <- param[2]
   }
   else if (method == 'mm') {
      param <- fit_GPD_MOM(x)
      xi <- param[1]
      scale <- param[2]
   }
   else if (method == 'bayesian') {
      out <-  MCMC4Extremes::gpdp(x, 0, 500)
      xi <- out$postmedian[2]
      scale <- out$postmedian[1]
   }
   else stop("The specified method is not implemented in this function.")
   
   print(ks.test(x, pGPD, shape=xi, scale=scale))
   print(ad.test(x, pGPD, shape=xi, scale=scale))
   
   return(c(xi, scale))
}


fit_Z <- function(data, show_qqplot=T) {
   #' Fonction qui paramétrise les précipitations normales selon une loi Gamma.
   #' L'hypothèse de tendance est testée via le critère de l'AIC.
   #' 
   #' @param data (tibble): Les données contenant les précipitations
   #'        non-extrêmes totales pour chacune des saisons printanières. Les
   #'        colonnes de ce dataFrame doivent être "total_prcp" et "year."
   #' @param show_qqplot (bool): Est-ce qu'on affiche le qqplot lors de 
   #'        l'exécution de la fonction ?
   #' 
   #' @return (list): La fonction quantile, la fonction de densité et la 
   #'        fonction de répartition de la loi paramétrée selon le maximum de 
   #'        vraisemblance.
   Z <- data$total_prcp
   years <- data$year
   
   cat("Mann-Kendall's test:", fill=T)
   print(MannKendall(Z))

   nll_gamma <- function(param, trend) {
      if (trend) {
         rate <- function(.year) {
            param[1] / (param[2] + param[3] * .year)
         }
      } else {
         rate <- function(.year) {
            param[1] / param[2]
         }
      }
      - sum(log(dgamma(Z, param[1], rate(years))))
   }
   
   nll_normal <- function(param, trend) {
      if (trend) {
         mu <- function(.year) {
            (param[2] + param[3] * .year)
         }
      } else {
         mu <- function(.year) {
             param[2]
         }
      }
      - sum(log(dnorm(Z, mu(years), param[1])))
   }
   
   sd_Z <- sd(Z)
   mean_Z <- mean(Z)
   shape <- mean_Z^2 / var(Z)
   
   optimizations <- list(
      Z_optim_gamma_trend <- optim(c(shape, mean_Z, 0), nll_gamma, gr=NULL, trend=T),
      Z_optim_normal_trend <- optim(c(sd_Z, mean_Z, 0), nll_normal, gr=NULL, trend=T),
      Z_optim_gamma_notrend <- optim(c(shape, mean_Z), nll_gamma, gr=NULL, trend=F),
      Z_optim_normal_notrend <- optim(c(sd_Z, mean_Z), nll_normal, gr=NULL, trend=F)
   )
   best <- which.min(c(
      aic_gamma_trend = 2 * (3 + Z_optim_gamma_trend$value),
      aic_normal_trend = 2 * (3 + Z_optim_normal_trend$value),
      aic_gamma_notrend = 2 * (2 + Z_optim_gamma_notrend$value),
      aic_normal_notrend = 2 * (2 + Z_optim_normal_notrend$value)
   ))
   best_par <- optimizations[[best]]$par
   best_dist <- c("gamma", "normal", "gamma", "normal")[best]
   
   trend <- if (best < 3) "with trend" else "whitout trend"
   cat("\nBest distribution:", best_dist, trend, fill=T)
   cat("Parameters:", best_par, fill=T)
   
   if (best_dist=='normal') {
      Z_sd <- best_par[1]
      
      mu <- function(.year) {
         if (length(best_par) == 3){
            return(best_par[2] + best_par[3] * .year)
         }
         return(best_par[2])
      }
      
      pZ <- function(z, .year) pnorm(z, mu(.year), Z_sd)
      dZ <- function(z, .year) dnorm(z, mu(.year), Z_sd)
      qZ <- function(p, .year) qnorm(p, mu(.year), Z_sd)
      
      if (show_qqplot) {
         car::qqPlot(Z, 'norm', mean=mu(.year), sd=Z_sd, lwd=0.5, id=F)
         abline(0, 1, col = 'red')
      }
      
   } else {
      shape <- best_par[1]
      
      rate <- function(.year) {
         if (length(best_par) == 3){
            return(shape / (best_par[2] + best_par[3] * .year))
         }
         return(shape / best_par[2])
      }
      
      pZ <- function(z, .year) pgamma(z, shape, rate(.year))
      dZ <- function(z, .year) dgamma(z, shape, rate(.year))
      qZ <- function(p, .year) qgamma(p, shape, rate(.year))
      
      if (show_qqplot) {
         car::qqPlot(Z, 'gamma', shape=shape, rate=rate(years),
                     ylab='Empirical quantiles', lwd=0.5, id=F)
         abline(0, 1, col = 'red')
      }
   }
   
   print(ks.test(Z, pZ, .year=years))
   print(ad.test(Z, pZ, .year=years))
   
   return(list(
      "quantiles"=qZ, "density"=dZ, "cdf"=pZ
   ))
}


fit_W <- function(data, heavy_tailed=F) {
   #' Fonction qui paramétrise la loi des temps inter-occurences des périodes 
   #' de pluies extrêmes selon l'une des trois distributions suivantes:
   #'  - La loi Géométrique
   #'  - La loi Binomiale négative
   #'  - Weibull discrétisée
   #' 
   #' @param data (tibble): Base de données tibble contenant les variables 
   #'        suivantes:
   #'        - W (float): Les temps inter-occurences des périodes de pluies
   #'        extrêmes.
   #'        - start_date (Date): La date de début des périodes de pluies
   #'        extrêmes.
   #'        - Duration (int): La durée des périodes de pluies extrêmes en
   #'        nombre de jours.
   #' @param heavy_tailed (bool): Est-ce que la distribution présente une queue
   #'        très lourde ? Dans l'affirmative, une loi GP sera paramétrée.
   #'        
   #' @return Une liste:
   #'        - density (func): Fonction de densité de probabilités paramétrée.
   #'        - distribution (func) Fonction de répartition des probabilités.
   #'        - quantiles (func): Fonction quantile paramétrée.
   #'        
   #'        Les paramètres utilisés pour ces deux fonctions sont:
   #'        - w (float): Les temps inter-occurences d'événements extrêmes.
   #'        - t0 (Date ou int): Date du début
   
   data_time <- data.extremes %>% select(W, start_date, t0, duration)
   W <- data_time$W %>% unlist(use.names = F)
   n <- length(W)
   
   t0 <- (data_time$t0 + data_time$duration) %>% as.numeric()
   max_t0 <- max(t0)
   t0 <- t0 / max_t0
   
   residual_times <- function(data_time) {
      t <- data_time$start_date %>% year() %>% unlist(use.names = F)
      min_year <- min(t)
      n_years <- max(t) - min_year + 1
      res <- numeric(n_years)
      .t0 <- numeric(n_years)
      for (i in 1:n_years) {
         m <- min_year + i - 1
         annual_times <- data_time %>% filter(year(start_date) == m)
         .t0[i] <- as.Date(paste0(as.character(m), '-04-01'))
         cum_time <- sum(annual_times$W, annual_times$duration)
         .t0[i] <- .t0[i] + cum_time
         res[i] <- 91 - cum_time
      }
      .t0 <- (.t0 - min(.t0)) / max_t0
      return(list("res"=res, "t0"=.t0, 'years'=seq(min_year, max(t), by=1)))
   }
   res <- residual_times(data_time)
   
   cat("Mann-Kendall's test:", fill=T)
   print(MannKendall(W))
   
   geom_nll <- function(param, trend) {
      
      mu <- function(.t0) {
         if (trend) {
            return(param[1] + param[2] * as.numeric(.t0))
         }
         return(param[1])
      }
      
      nll <- -sum(log(dgeom(W, 1 / mu(t0))))
      nll <- nll - sum(log(1 - pgeom(res$res, 1 / mu(res$t0))))
      return(nll)
   }
   
   nbin_nll <- function(param, trend) {
      shape <- param[1]
      
      mu <- function(.t0) {
         if (trend) {
            return(param[2] + param[3] * as.numeric(.t0))
         }
         return(param[2])
      }
      
      nll <- -sum(log(dnbinom(W, shape, mu=mu(t0))))
      nll <- nll - sum(log(1 - pnbinom(res$res, shape, mu=mu(res$t0))))
      return(nll)
   }
   
   weibull_nll <- function(param, trend) {
      shape <- param[1]
      
      scale <- function(.t0) {
         if (trend) {
            return(exp(param[2] + param[3] * as.numeric(.t0)) / gamma(1 + 1 / shape))
         }
         return(param[2] / gamma(1 + 1 / shape))
      }
      
      cdf <- function(w, scales) {
         mapply(pweibull, q=w, scale=scales, shape=shape)
      }
      
      pmf <- function(w, scales) {
         f <- function(.w, .sig) {
            pweibull(.w+1, shape, .sig) - pweibull(.w, shape, .sig)
         }
         mapply(f, w, scales)
      }
      
      if (min(shape, scale(t0)) <= 0) return(NA)
      nll <- -sum(log(pmf(W, scale(t0))))
      nll <- nll - sum(log(1 - cdf(res$res, scale(res$t0))))
      return(nll)
   }
   
   pareto_nll <- function(param, trend) {
      shape <- param[1]
      if (trend) {
         scale <- function(.t) {
            return(exp(param[2] + param[3] * as.numeric(.t)))
         }
      } else {
         scale <- function(.t) param[2]
      }
      
      cdf <- function(w, scales) {
         mapply(pPar, q=w, scale=scales, shape=shape)
      }
      
      pmf <- function(w, scales) {
         f <- function(.w, .sig) {
            if (.w == 0) return(0)
            pPar(.w+1, shape, .sig) - pPar(.w, shape, .sig)
         }
         mapply(f, w, scales)
      }
      
      if (min(shape, scale(t0)) <= 0) return(0)
      
      nll <- -sum(log(pmf(W, scale(t0))))
      nll <- nll - sum(log(1 - cdf(res$res, scale(res$t0))))
      return(nll)
   }

   if(heavy_tailed) {
      
      fit_W_pareto <- function(init_par, trend, quiet=T, seed=2021) {
         require(MHadaptive)
         
         prior_reg <- function(param, trend) {
            
            prior_shape <- dgamma(param[1], 3, 1, log = T)
            
            if (trend) {
               prior_a <- dnorm(param[2], init_par[2], 2, log=T)
               prior_b <- dnorm(param[3], 0, 1, log=T)
               prior_scale <- prior_a + prior_b
            } else {
               prior_scale <- dgamma(param[2], 15, 1/5, log=T)
            }
            return(prior_shape + prior_scale)
         }
         
         li_func <- function(param, trend) {
            return(prior_reg(param, trend) - pareto_nll(param, trend))
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
            iterations = 5e+4,
            burn_in = 1e+4,
            quiet=T,
            trend=trend
         )
         init_par <- apply(mcmc_r$trace, 2, getmode)
         cat("prelim_pars:", init_par, fill=T)
         
         mcmc_r <- Metro_Hastings(
            li_func = li_func,
            pars = init_par,
            prop_sigma = mcmc_r$prop_sigma,
            par_names = par_names,
            iterations = 5e+4,
            burn_in = 1e+4,
            quiet=T,
            trend=trend
         )
         par <- apply(mcmc_r$trace, 2, getmode)
         names(par) <- par_names
         value <- pareto_nll(par, trend)
         cat('Pareto, Trend=', trend, ', value:', value, fill=T)
         cat('parameters', par, fill=T)
         return(list(
            "par" = par,
            "value" = value
         ))
         
         unloadNamespace("MHadaptive")
         unloadNamespace("MASS")
      }
      
      # No trend
      init_par <- c(2, mean(W))
      cat("Pareto initial parameters:", init_par, fill=T)
      pareto_notrend <- fit_W_pareto(init_par = init_par, trend=F)
      
      # With trend
      init_par <- pareto_notrend$par
      pareto_trend <- fit_W_pareto(
         init_par = c(init_par[1], log(init_par[2]), 0),
         trend=T)
      
   } else {
      pareto_notrend <- list(par = NA, value = NA)
      pareto_trend <- list(par = NA, value = NA)
   }

   optimizations <- list(
      # With trend
      optim_geom_trend <- optim(c(mean(W), 0), geom_nll, gr=NULL, trend=T),
      optim_nbin_trend <- optim(c(1, mean(W), 0), nbin_nll, gr=NULL, trend=T),
      optim_weib_trend <- optim(c(1, log(mean(W)), 0), weibull_nll, gr=NULL, trend=T),
      pareto_trend,
      # Without trend
      optim_geom_notrend <- list(
         "par" = mean(W),
         "value" = geom_nll(mean(W), trend=F)),
      optim_nbin_notrend <- optim(c(1, mean(W)), nbin_nll, gr=NULL, trend=F),
      optim_weib_notrend <- optim(c(1, mean(W)), weibull_nll, gr=NULL, trend=F),
      pareto_notrend
   )
   
   aics <- cbind(
      # With trend
      aic_geom_trend = 2 * (2 + optim_geom_trend$value),
      aic_nbin_trend = 2 * (3 + optim_nbin_trend$value),
      aic_weib_trend = 2 * (3 + optim_weib_trend$value),
      aic_pare_trend  = 2 * (3 + pareto_trend$value),
      # Without trend
      aic_geom_notrend = 2 * (1 + optim_geom_notrend$value),
      aic_nbin_notrend = 2 * (2 + optim_nbin_notrend$value),
      aic_weib_notrend = 2 * (2 + optim_weib_notrend$value),
      aic_pare_notrend  = 2 * (2 + pareto_notrend$value)
   )
   best <- which.min(aics)
   print(aics)
   
   best_par <- optimizations[[best]]$par
   best_dist <- rep(
      c('geometric', 'Negative binomial', 'Weibull', "Pareto"), 2)[best]
   
   trend <- if (best < 5) "with trend" else "whitout trend"
   cat("\nBest distribution:", best_dist, trend, fill=T)
   cat("Parameters:", best_par, fill=T)
   
   if (best_dist=='geometric') {
      
      mu <- function(.t0) {
         if (length(best_par) == 2){
            return(best_par[1] + best_par[2] / max_t0 * as.numeric(.t0))
         }
         return(best_par[1])
      }
      
      pW <- function(w, t0) pgeom(w, 1 / mu(t0))
      dW <- function(w, t0) dgeom(w, 1 / mu(t0))
      qW <- function(p, t0) qgeom(p, 1 / mu(t0))
      
   } else if (best_dist=='Negative binomial') {
      shape <- best_par[1]
      
      mu <- function(.t0) {
         if (length(best_par) == 3){
            return(best_par[2] + best_par[3] / max_t0 * as.numeric(.t0))
         }
         return(best_par[2])
      }
      
      pW <- function(w, t0) pnbinom(w, shape, mu=mu(t0))
      dW <- function(w, t0) dnbinom(w, shape, mu=mu(t0))
      qW <- function(p, t0) qnbinom(p, shape, mu=mu(t0))
      
   } else if (best_dist=='Weibull') {
      shape <- best_par[1]
      
      scale <- function(.t0) {
         if (length(best_par) == 3) {
            return(exp(best_par[2] + best_par[3] / max_t0 * as.numeric(.t0)) / gamma(1 + 1 / shape))
         }
         return(best_par[2] / gamma(1 + 1 / shape))
      }
      
      pW <- function(w, t0) pweibull(w, shape, scale(t0))
      qW <- function(p, t0) ceiling(qweibull(p, shape, scale(t0)))
      dW <- function(w, t0) {
         sapply(w, function(.w) {
            if (.w == 0) 0
            else pD(.w, t0) - pD(.w - 1, t0) 
            })
      }
      
   } else { # best_dist=='Pareto'
      
      if (length(best_par) == 3) {
         shape <- best_par["Shape (Intercept)"]
         scale <- function(.t) {
            a <- best_par["Scale (Intercept)"]
            b <- best_par["Scale time"]
            return(a + b / max_t0 * as.numeric(.t))
         }
      } else {
         shape <- best_par['shape']
         scale <- function(.t) best_par['scale']
      }
      
      pW <- function(w, t0) pPar(w, shape, scale(t0))
      qW <- function(p, t0) ceiling(qPar(p, shape, scale(t0)))
      dW <- function(w, t0) {
         sapply(w, function(.w) {
            if (.w == 0) 0
            else pD(.w, t0) - pD(.w - 1, t0) 
         })
      }
   }
   
   return(list(
      "density" = dW,
      "distribution" = pW,
      "quantiles" = qW,
      'params' = best_par
   ))
}


fit_D <- function(data) {
   #' Fonction qui paramétrise la loi de la durée des périodes 
   #' de pluies extrêmes selon l'une des trois distributions suivantes:
   #'  - La loi exponentielle
   #'  - La loi Weibull
   #'  - La loi Gamma
   #' 
   #' @param data (tibble): Base de données tibble contenant les variables 
   #'        suivantes:
   #'        - W (float): Les temps inter-occurences des périodes de pluies
   #'        extrêmes.
   #'        - start_date (Date): La date de début des périodes de pluies
   #'        extrêmes.
   #'        - Duration (int): La durée des périodes de pluies extrêmes en
   #'        nombre de jours.
   #'        
   #' @return Une liste:
   #'        - density (func): Fonction de densité de probabilités paramétrée.
   #'        - distribution (func) Fonction de répartition des probabilités.
   #'        - quantiles (func): Fonction quantile paramétrée.
   #'        
   #'        Les paramètres utilisés pour ces deux fonctions sont:
   #'        - x (float): La durée des événements extrêmes.
   #'        - t0 (Date ou int): Date du début

   data_time <- data.extremes %>% select(W, start_date, duration, t0)
   D <- data_time$duration %>% as.numeric() %>% unlist(use.names = F)
   t0 <- data_time$t0
   max_t0 <- max(t0)
   t0 <- t0 / max_t0

   cat("Mann-Kendall's test:", fill=T)
   print(MannKendall(D))
   
   geom_nll <- function(param, trend) {
      
      mu <- function(.t0) {
         if (trend) {
            return(exp(param[1] + param[2] * as.numeric(.t0)))
         }
         return(param[1])
      }
      
      nll <- -sum(log(dgeom(D, 1 / mu(t0))))
      return(nll)
   }
   
   nbin_nll <- function(param, trend) {
      shape <- param[1]
      
      mu <- function(.t0) {
         if (trend) {
            return(exp(param[2] + param[3] * as.numeric(.t0)))
         }
         return(param[2])
      }
      
      nll <- -sum(log(dnbinom(D, shape, mu=mu(t0))))
      return(nll)
   }
   
   pois_nll <- function(param, trend) {
      mu <- function(.t0) {
         if (trend) {
            return(exp(param[1] + param[2] * as.numeric(.t0)))
         }
         return(param[1])
      }
      
      nll <- -sum(log(dpois(D, mu(t0))))
      return(nll)
   }
   
   weibull_nll <- function(param, trend) {
      shape <- param[1]
      
      scale <- function(.t0) {
         if (trend) {
            return(exp(param[2] + param[3] * as.numeric(.t0)) / gamma(1 + 1 / shape))
         }
         return(param[2] / gamma(1 + 1 / shape))
      }
      
      pmf <- function(d, scales) {
         f <- function(.d, .sig) {
            if (.d == 0) return(0)
            pweibull(.d, shape, .sig) - pweibull(.d-1, shape, .sig)
         }
         mapply(f, d, scales)
      }
      
      if (min(shape, scale(t0)) <= 0) return(NA)
      nll <- -sum(log(dweibull(D, shape, scale(t0))))
      return(nll)
   }
   
   optimizations <- list(   
      # With trend
      optim_geom_trend <- optim(c(log(mean(D)), 0), geom_nll, gr=NULL, trend=T),
      optim_nbin_trend <- optim(c(1, log(mean(D)), 0), nbin_nll, gr=NULL, trend=T),
      optim_pois_trend <- optim(c(log(mean(D)), 0), pois_nll, gr=NULL, trend=T),
      optim_weib_trend <- optim(c(1, log(mean(D)), 0), weibull_nll, gr=NULL, trend=T),
      # Without trend
      optim_geom_notrend <- list(
         "par" = mean(D),
         "value" = geom_nll(mean(D), trend=F)),
      optim_nbin_notrend <- optim(c(1, mean(D)), nbin_nll, gr=NULL, trend=F),
      optim_pois_notrend <- list(
         "par" = mean(D),
         "value" = pois_nll(mean(D), trend=F)),
      optim_weib_notrend <- optim(c(1, mean(D)), weibull_nll, gr=NULL, trend=F)
   )
   
   aics <- cbind(
      # With trend
      aic_geom_trend = 2 * (2 + optim_geom_trend$value),
      aic_nbin_trend = 2 * (3 + optim_nbin_trend$value),
      aic_pois_trend = 2 * (2 + optim_pois_trend$value),
      aic_weib_trend = 2 * (3 + optim_weib_trend$value),
      # Without trend
      aic_geom_notrend = 2 * (1 + optim_geom_notrend$value),
      aic_nbin_notrend = 2 * (2 + optim_nbin_notrend$value),
      aic_pois_notrend = 2 * (1 + optim_pois_notrend$value),
      aic_weib_notrend = 2 * (2 + optim_weib_notrend$value)
   )
   print(aics)
   best <- which.min(aics)

   best_par <- optimizations[[best]]$par
   best_dist <- rep(
      c('geometric', 'Negative binomial', 'Poisson', 'Weibull'), 2)[best]
   
   trend <- if (best < 5) "with trend" else "whitout trend"
   cat("\nBest distribution:", best_dist, trend, fill=T)
   cat("Parameters:", best_par, fill=T)
   
   if (best_dist=='geometric') {
      
      mu <- function(.t0) {
         if (length(best_par) == 2){
            return(exp(best_par[1] + best_par[2] / max_t0 * as.numeric(.t0)))
         }
         return(best_par[1])
      }
      
      pD <- function(d, t0) pgeom(d, 1 / mu(t0))
      dD <- function(d, t0) dgeom(d, 1 / mu(t0))
      qD <- function(p, t0) qgeom(p, 1 / mu(t0))
      
   } else if (best_dist=='Negative binomial') {
      shape <- best_par[1]
      
      mu <- function(.t0) {
         if (length(best_par) == 3){
            return(exp(best_par[2] + best_par[3] / max_t0 * as.numeric(.t0)))
         }
         return(best_par[2])
      }
      
      pD <- function(d, t0) pnbinom(d, shape, mu=mu(t0))
      dD <- function(d, t0) dnbinom(d, shape, mu=mu(t0))
      qD <- function(p, t0) qnbinom(p, shape, mu=mu(t0))
      
   } else if (best_dist=='Poisson') {
      
      mu <- function(.t0) {
         if (length(best_par) == 2){
            return(exp(best_par[1] + best_par[2] / max_t0 * as.numeric(.t0)))
         }
         return(best_par[1])
      }
      
      pD <- function(d, t0) ppois(d, mu(t0))
      dD <- function(d, t0) dpois(d, mu(t0))
      qD <- function(p, t0) qpois(p, mu(t0))
      
   } else {  # best_dist=='Weibull'
      shape <- best_par[1]
      
      scale <- function(.t0) {
         if (length(best_par) == 3) {
            return(exp(best_par[2] + best_par[3] / max_t0 * as.numeric(.t0)) / gamma(1 + 1 / shape))
         }
         return(best_par[2] / gamma(1 + 1 / shape))
      }
      
      pD <- function(d, t0) pweibull(d, shape, scale(t0))
      qD <- function(p, t0) ceiling(qweibull(p, shape, scale(t0)))
      dD <- function(d, t0) {
         sapply(d, function(.d) {
            if (.d == 0) 0
            else pD(.d, t0) - pD(.d - 1, t0) 
         })
      }
   }
   return(list(
      "density" = dD,
      "distribution" = pD,
      "quantiles" = qD,
      'params' = best_par
   ))
}


calculate_interoccurence_times <- function(data, n_print=20) {
   #' Fonction qui calcule les temps séparant les événements extrêmes.
   #' 
   #' @param data (tibble): DataFrame contenant les excès de seuil.
   #'        Nécessite les colonnes suivantes:
   #'        - excedence (float): Les excédents de seuil
   #'        - time_since_last (int): Nombre de jours depuis la dernière 
   #'          période de pluie.
   #'        - duration (int): Nombre de jours qu'a duré la période de pluie.
   #'        - start_date (date): Date du début de caque période de pluie.
   #' @param n_print (int): Nombre d'observations à afficher pour valider que
   #'        le clustering et le calcul des temps inter-occurence se sont
   #'        bien fait. Si n_print=0, alors aucun tableau n'est affiché.
   #'        
   #' @return data: le DataFrame tibble filtré pour ne conserver que les
   #'         précipitations extrêmes et dont on a ajouté une colonne pour les
   #'         temps inter-occurences.
   
   nb_extremes <- data %>% filter(excedence > 0) %>% nrow()
   n <- nrow(data)
   pb <- txtProgressBar(min=1, max=nb_extremes, style=3)
   W <- numeric(nb_extremes)
   j <- 1
   
   if (data$excedence[1] == 0) {
      W[j] <-
         sum(W[j], data$time_since_last[1], data$duration[1], na.rm = T)
   } else {
      W[j] <-  sum(W[j], data$time_since_last[1], na.rm = T)
      j <- j + 1
   }
   
   for (i in 2:n) {
      if (year(data$start_date[i]) == year(data$start_date[i - 1])) {
         if (data$excedence[i] == 0) {
            W[j] <-
               sum(W[j], data$time_since_last[i], data$duration[i], na.rm = T)
         } else {
            W[j] <- sum(W[j], data$time_since_last[i], na.rm = T)
            j <- j + 1
         }
      } else {
         W[j] <- data$time_since_last[i]
         if (data$excedence[i] > 0) {
            j <- j + 1
         }
      }
      setTxtProgressBar(pb, j)
      if (j > nb_extremes) break
   }
   
   
   if (n_print > 0) {
      cat("\n Premiers W:", head(W, 3), fill=T)
      print(
         rbind(
            "Year" = head(year(data$start_date), n_print),
            "Time_since_last" = head(data$time_since_last, n_print),
            "Durations" = head(data$duration, n_print),
            "Excedence" = head(data$excedence > 0, n_print)
         ))
   }
   
   data <- data %>%
      filter(excedence>0) %>%
      select(start_date, duration, excedence) %>%
      mutate(W = W, .after = duration)
   
   return(data)
}


calculate_excedences <- function(data, treshold, calculate_W=T, n_print=0) {
   #' Fonction qui calcule les excédent de seuil pour les précipitations
   #' agrégées d'une base de données.
   #' 
   #' @param data (tibble): DataFrame tibble contenant les données de 
   #'        précipitation. Elle doit contenir une colonne nommée prcp_agg qui
   #'        correspond aux précipitations aggrégées par période de pluie 
   #'        continue.
   #' @param treshold (float): Le seuil qui sépare les valeurs extrêmes des 
   #'        valeurs normales. Si treshold=0, celui-ci est déterminé 
   #'        automatiquement.
   #' @param calculate_W (bool): Est-ce qu'on désire calculer le temps 
   #'        séparant les valeurs extrêmes ?
   #' @param n_print (int): Nombre d'observations à afficher pour valider que
   #'        le clustering et le calcul des temps inter-occurence se sont
   #'        bien fait. Si n_print=0, alors aucun tableau n'est affiché.
   #'        
   #' @return data (tibble): le DataFrame tibble offert en entrée auquel on 
   #'        ajoute deux colonnes:
   #'        - excedence : Excendent au dessus du seuil (X-u|X>u)
   #'        - normal : la normale (min(X, u)),
   #'        où u désigne le seuil.
   if (treshold == 0) {
      set.seed(42)
      p_vec <- rnorm(100, mean=0.86, sd=0.04) %>% sort()
      u_vec <- quantile(data$prcp_agg, probs = p_vec) 
      u <- treshold_selection(data, u_vec, method='bader', gf_test = 'ad', verbose=2)
   }
   data <- data %>% mutate(excedence = pmax(prcp_agg - treshold, 0))
   
   nb_extremes <- data %>% filter(excedence > 0) %>% nrow()
   
   if (calculate_W) {
      data <- calculate_interoccurence_times(data, n_print)
   }
   data['t0'] <- (data$start_date - min(data$start_date)) %>% as.numeric()
   
   return(data)
}


treshold_selection <- function(data, u_vec, method='bader', gf_test='ad', 
                               nsim=NULL, verbose=2, seed=2021) {
   #' Fonction qui optimise le seuil u qui sépare les valeurs de précipitations
   #' extrêmes de celles qui sont normales. S'appuie sur Northrop et al. 2017.
   #' 
   #' @param data (tibble): DataFrame tibble contenant les données de 
   #'        précipitation. Elle doit contenir une colonne nommée prcp_agg qui
   #'        correspond aux précipitations aggrégées par période de pluie 
   #'        continue.
   #' @param u_vec (float): Vecteur des seuil à essayer.
   #' @param method (str) : L'une des options suivantes:
   #'        - 'northrop': Méthode proposée par Northrop et al. 2017. 
   #'        Voir ithresh du package threshr.
   #'        - 'bader': méthode proposée par Bader et al. 2018.
   #'        Voir gpdSeqTests du package eva.
   #' @param gf_test (str): Argument utilisé avec la méthode de sélection de Bader.
   #'        Test d'adéquation utilisé pour sélectionner le meilleur seuil.
   #'        L'une des options suivantes:
   #'        - 'ad' : Anderson-Darling test
   #'        - 'cvm' : Cramer-von Mises Test
   #'        - 'pbscore' : Parametric Bootstrap Score Test
   #'        - 'multscore': Fast weighted bootstrap alternative to the 
   #'          parametric bootstrap procedure for the Generalized Pareto score 
   #'          test.
   #'        - 'imasym' : Asymptotic Adjusted Information Matrix (IM) Test. 
   #'          Runs the IM Test using bootstrap estimated covariance matrix. 
   #'          Asymptotically (in sample size) follows the F(3, bootnum - 3) 
   #'          distribution.
   #'        - 'impb': Bootstrapped Information Matrix (IM) Test. Runs the IM
   #'          Test using a two-step iterative procedure, to boostrap the 
   #'          covariance estimate and critical values.
   #' @param nsim (int): Le nombre d'itération du bootstrap
   #' @param verbose (int): 
   #'        - Si verbose = 0, aucune statistique d'adéquation ni aucun graphique
   #'         n'est affiché.
   #'         - Si verbose = 1, affiche le résultat des tests d'adéquation de
   #'         Kolmogorov-Smirnov et d'Anderson-Darling.
   #'         - si verbose = 2, affiche en plus le graphique quantiles-à-
   #'         quantiles
   #'         - si verbose > 2, affiche en plus le graphique de seuil optimal 
   #'         présenté avec la méthode plot.ithresh du package threshr.
   #' @param seed (int): L'ancrage utilisé pour initialiser l'algorithme MCMC 
   #'         utilisé pour entraîner le modèle servant à optimiser le seuil.
   #' 
   #' @return float: retourne le seuil optimal pour maximiser l'adéquation de 
   #'         la queue de distribution au modèle GPD.  
   set.seed(seed)
   
   X <- data$prcp_agg %>% unlist(use.names = F)
   X <- X + rnorm(length(X), 0, 1e-6)
   
   if (method=='northrop'){
      u_cv <- ithresh(data = X, u_vec = u_vec, n=nsim, n_v=5)
      u <- summary(u_cv)[3]
      
      params <- fit_GPD(X[X>u]-u, method='pwm')
      
      num_above <- data %>% filter(prcp_agg > u) %>% nrow()
      cat(paste("Threshold selection method:", method), fill=T)
      
      if (verbose > 2) {
         plot(u_cv)
      }
   } else if (method=='bader'){
      u_cv_balder <- gpdSeqTests(
         data = X,
         thresholds = u_vec,
         method = gf_test,
         nsim = nsim,
         allowParallel = T,
         numCores = 6
      )
      best <- u_cv_balder %>% top_n(1, ForwardStop)
      params <- best %>% select(est.shape, est.scale) %>% unlist()
      u <- best$threshold
      
      num_above <- best$num.above
      cat(paste("Threshold selection method:", method), fill=T)
      cat(paste("Threshold selection criterion:", gf_test), fill=T)
      
   } else {
      stop(paste0("The method", method, "is not implemented",
                  "Select either 'northrop' or 'bader'."))
   }

   if (verbose > 0){
      cat(paste("Associated percentile:", 
                  round(ecdf(data$prcp_agg)(u) * 100, 2)), fill=T)
      cat(paste("Number of observations over threshold:", num_above, "\n"))
      
      print(paste("Anderson-Darling p-value:", gpdAd(X[X>u]-u)$p.value), quote=F)
      print(paste("Cramer-von Mises p-value:", gpdCvm(X[X>u]-u)$p.value), quote=F)

      if (verbose > 1){
         cat("GPD parameters:", params, fill=T)
         car::qqPlot(X[X>u]-u, 'GPD',
                     shape=params[1],
                     scale=params[2],
                     lwd=0.5, id=F, cex=0.8,
                     ylab='Empirical quantiles')
         abline(0,1,col='red')
      }
   }
   
   return(u)
}

verify_serial_dependence <- function(data, max_lag=1, seed=2021, return_xtable=F) {
   #' Fonction qui évalue l'hypothèse nulle d'indépendence séquentielle selon 
   #' l'approche proposée dans Genest et Rémillard 2004.
   #' 
   #' @param data (tibble): Un DataFrame contenant les colonnes suivantes:
   #'        - start_date (date): Dates de début des périodes de pluie
   #'        - prcp_agg (float): La quantité de pluie totale tombée par période
   #' @param max_lag (int): Le nombre d'années de décalage (lag) pour calculer 
   #'        la corrélation.
   #' @param seed (int): Ancrage pour le générateur de nombre pseudo-aléatoires.
   #'        Celui-ci est utilisé pour prendre des années aléatoirement afin de
   #'        réaliser le test su chacune de ces années indépendamment.
   #' @param return_xtable (bool): Est-ce que la fonction retourne le tableau
   #'        des p-values sous forme de tableau xtable afin de le présenter
   #'        dans un rapport LaTeX ?
   #'        
   #' @return None: Affiche simplement à la console les résultats calculés avec 
   #'         les fonctions serialIndepTest et multSerialIndepTest du package 
   #'         copula.
   years <- data$start_date %>% year() %>% unique()
   
   cat(" Number of observations by year:", fill=T)
   n_obs_per_year <- data %>% 
      mutate(year=year(start_date)) %>%
      group_by(year) %>%
      summarise(nobs=length(prcp_agg), .groups='drop') %>%
      select(nobs) %>% 
      unlist() %>% summary() %>% print()
   
   n_years <- length(years)
   
   data <- data %>% mutate(duration = as.numeric(duration))
   
   pvalues <- array(dim=c(n_years, 6))
   
   i <- 1
   pb <- txtProgressBar(min=1, max=n_years, style = 3)
   for (year in years) {
      .data <- data %>% filter(year(start_date) == year)
      
      SITS <- serialIndepTestSim(nrow(.data), max_lag, verbose = F)
      
      X <- .data %>% select(prcp_agg) %>% unlist()
      pvalues[i, 1] <- serialIndepTest(X, SITS)$pvalues
      
      W <- .data %>% select(time_since_last) %>% unlist() %>% na.omit()
      pvalues[i, 2] <- serialIndepTest(W, SITS)$pvalues
      
      D <- .data %>% select(duration) %>% unlist()
      pvalues[i, 3] <- serialIndepTest(D, SITS)$pvalues
      
      XW <- .data %>% select(prcp_agg, time_since_last) %>% na.omit()
      pvalues[i, 4] <- multSerialIndepTest(XW, lag.max = max_lag, verbose = F)$pvalues
      
      XD <- .data %>% select(prcp_agg, duration) %>% na.omit()
      pvalues[i, 5] <- multSerialIndepTest(XD, lag.max = max_lag, verbose = F)$pvalues
      
      WD <- .data %>% select(time_since_last, duration) %>% na.omit()
      pvalues[i, 6] <- multSerialIndepTest(WD, lag.max = max_lag, verbose = F)$pvalues
      
      setTxtProgressBar(pb, i)
      i <- i+1
   }
   
   pvalues_stats <- apply(pvalues, 2, summary) %>% t()
   rownames(pvalues_stats) <- list("X", "W", "D", "(X,W)", "(X,D)", "(W,D)")

   cat("\n Serial dependence test:", fill=T)
   print(pvalues_stats)
   print(quantile(pvalues, 0.05))
   
   if (return_xtable) {
      cat("\n", fill = T)
      return(xtable::xtable(pvalues_stats, digits=3))
   }
}


calculate_correlation <- function(data, alternative='two-sided'){
   #' Fonction qui prend la base de données et qui retourne les mesures de
   #' corrélation de Spearman et de Kendall entre les variables de l'excédent
   #' de seuil et du temps séparant deux événements.
   #' 
   #' @param data (tibble): Un dataFrame tibble contenant deux colonnes dont on
   #'        voudrait vérifier l'indépendance.
   #' @param alternative (str): correspond à la forme de l'hypothèse alternative.
   #'        peut prendre l'une des valeurs suivantes:
   #'        - "two-sided"  H1 : \rho != 0
   #'        - "greater"  H1 : \rho > 0
   #'        - "lower"  H1 : \rho < 0
   #'        
   #' @return list: le rho de spearman et le tau de Kendall empirique.       
   correlations <- lapply(
      c('spearman', 'kendall'),
      cor,
      x = data[,1] %>% unlist(use.names = F),
      y = data[,2] %>% unlist(use.names = F),
      use = "pairwise.complete.obs"
   )
   names(correlations) <- c('spearman', 'kendall')
   
   # Test d'indépendance
   mh_stat <- sqrt(nrow(data) - 1) * correlations$spearman
   if (alternative=='lower') {
      mh_pvalue <- pnorm(mh_stat)
   } else if (alternative == 'greater') {
      mh_pvalue <- 1 - pnorm(mh_stat)
   } else {
      mh_pvalue <- 1 - pchisq(mh_stat^2, 1)
   }
   
   cat("Mantel-Haensel independence test", fill=T)
   cat("statistic:", mh_stat, ', p-value', mh_pvalue, fill=T)
   cat("alternative:", alternative, fill=T)
   cat("\n")
   
   return(correlations)
}


eval_dependence_trend <- function(data, variable='W', stride=50) {
   #' Fonction qui permet de voir s'il existe une tendance dans le niveau de 
   #' dépendance liant les variables des temps inter-occurence et de la 
   #' sévérité des pluies extrêmes. 
   #' La méthode se base sur {Chebana 2021 - Multivariate non-stationary 
   #' hydrological frequency analysis}.
   #' 
   #' @param data (tibble): DataFrame contenant les colonnes suivantes:
   #'        - W (int): Le temps séparant deux périodes de pluies extrêmes (en jours).
   #'        - excedence (float): l'excédent de pluie tombée par rapport au 
   #'          seuil de la méthode POT.
   #' @param variable (str): Le nom de la colonne utilisée pour vérifier une 
   #'        tendance dans la dépendance avec la variable d'excédent de seuil.
   #' @param stride (int): Le nombre d'observations contenues dans un pas de
   #'        convolution.
   #' 
   #' @return Les rhos de spearman calculés       
   W <- data %>% select(all_of(variable)) %>% unlist(use.names = F)
   X <- data$excedence
   n <- nrow(data)
   n_conv <- n - stride + 1
   rhos <- numeric(n_conv)
   for (i in 1:n_conv) {
      plage <- i:(i + stride - 1)
      rhos[i] <- cor(W[plage], X[plage], method = 'spearman')
   }
   plot(rhos, main='Dependence trend')
   print(bbsmk(rhos))
}


copula_selection.XZ <- function(data, dep_trend=F, family_set=NA,
                                gofTest='None', seed=42) {
   #' Fonction qui sélectionne la copule la plus adéquate parmi toutes les 
   #' copules programmées dans la librairie VineCopula.
   #' 
   #' @param data (tibble): DataFrame tibble contenant les colonnes suivantes:
   #'        - X : Les précipitations extrêmes totales pour chaque année.
   #'        - Z : Les précipitations normales totales pour chaque année.
   #' @param dep_trend (bool): Y-a-t'il des signes de tendances dans la dépendance ?
   #' @param family_set (int): Vecteur des choix possibles. Chaque famille 
   #'        implantée dans la librairie VineCopula comporte un chiffre de 
   #'        référence comme suit:
   #'        0 = independence copula 
   #         1 = Gaussian copula 
   #         2 = Student t copula (t-copula) 
   #         3 = Clayton copula 
   #         4 = Gumbel copula 
   #         5 = Frank copula 
   #         6 = Joe copula 
   #         7 = BB1 copula 
   #         8 = BB6 copula 
   #         9 = BB7 copula 
   #         10 = BB8 copula 
   #         13 = rotated Clayton copula (180 degrees;  survival Clayton'') \cr `14` = rotated Gumbel copula (180 degrees; survival Gumbel”) 
   #         16 = rotated Joe copula (180 degrees;  survival Joe'') \cr `17` = rotated BB1 copula (180 degrees; survival BB1”)
   #         18 = rotated BB6 copula (180 degrees;  survival BB6'')\cr `19` = rotated BB7 copula (180 degrees; survival BB7”)
   #         20 = rotated BB8 copula (180 degrees; “survival BB8”)
   #         23 = rotated Clayton copula (90 degrees) 
   #         '24' = rotated Gumbel copula (90 degrees) 
   #         '26' = rotated Joe copula (90 degrees) 
   #         '27' = rotated BB1 copula (90 degrees) 
   #         '28' = rotated BB6 copula (90 degrees) 
   #         '29' = rotated BB7 copula (90 degrees) 
   #         '30' = rotated BB8 copula (90 degrees) 
   #         '33' = rotated Clayton copula (270 degrees) 
   #         '34' = rotated Gumbel copula (270 degrees) 
   #         '36' = rotated Joe copula (270 degrees) 
   #         '37' = rotated BB1 copula (270 degrees) 
   #         '38' = rotated BB6 copula (270 degrees) 
   #         '39' = rotated BB7 copula (270 degrees) 
   #         '40' = rotated BB8 copula (270 degrees) 
   #         '104' = Tawn type 1 copula 
   #         '114' = rotated Tawn type 1 copula (180 degrees) 
   #         '124' = rotated Tawn type 1 copula (90 degrees) 
   #         '134' = rotated Tawn type 1 copula (270 degrees) 
   #         '204' = Tawn type 2 copula 
   #         '214' = rotated Tawn type 2 copula (180 degrees) 
   #         '224' = rotated Tawn type 2 copula (90 degrees) 
   #         '234' = rotated Tawn type 2 copula (270 degrees)
   #' 
   #' @param gofTest (str): L'une de ces options suivantes:
   #'        - "white" = Test d'adéquation basé sur la matrice d'information de 
   #'        White (White, 1982;  Huang and Prokhorov (2011))
   #'        - "Kendall" = Test d'adéquation basé sur le processus de Kendall  
   #'        (Wang and Wells, 2000; Genest et al., 2006)
   #'        - "None" = Aucun test d'adéquation n'est effectué.
   #' @param seed (int): Ancrage de simulation pour pouvoir reproduire les 
   #'        résultats.
   #' 
   #' @return un objet de la classe BiCop (voir le package VineCopula).
   set.seed(seed)
   
   XZ <- data %>% select(-year) %>% na.omit()
   years <- data$year
   # white_noise <- rnorm(2 * nrow(XZ), mean = 0, sd=1e-8) %>% matrix(ncol = 2)
   # XZ <- XZ + white_noise
   uu <- pobs(XZ)
   
   # BiCopCompare(uu[,1], uu[,2])
   bestCopula <- BiCopSelect(uu[,1], uu[,2], familyset = family_set, indeptest = T)
   
   if (dep_trend){
      copula_nll.1 <- function(params) {
         nll_copula <- function(year) {
            u <- uu[which(years == year), ]
            par <- params[1] + params[2] * year
            return(-sum(log(
               BiCopPDF(uu[, 1], uu[, 2], bestCopula$family, par = par)
            )))
         }
         return(sum(sapply(years, nll_copula)))
      }   
      copula_nll.2 <- function(params) {
         nll_copula <- function(year) {
            u <- uu[which(years == year), ]
            par <- params[1] + params[2] * year**2
            return(-sum(log(
               BiCopPDF(uu[, 1], uu[, 2], bestCopula$family, par = par)
            )))
         }
         return(sum(sapply(years, nll_copula)))
      }   
      copula_nll.3 <- function(params) {
         nll_copula <- function(year) {
            u <- uu[which(years == year),]
            par <- params[1] + params[2] * year + params[3] * year**2
            return(-sum(log(
               BiCopPDF(uu[, 1], uu[, 2], bestCopula$family, par = par)
            )))
         }
         return(sum(sapply(years, nll_copula)))
      }
      
      optim.1 <- optim(c(bestCopula$par, 0), copula_nll.1)
      optim.2 <- optim(c(bestCopula$par, 0), copula_nll.2)
      optim.3 <- optim(c(optim.1$par, 0), copula_nll.3)
      
      .aic <- c(
         bestCopula$AIC,
         aic.1 <- 2 * (2 - optim.1$value),
         aic.2 <- 2 * (2 - optim.2$value),
         aic.2 <- 2 * (3 - optim.3$value)
      )
      best <- which.min(.aic)
      params <- list(bestCopula$par, optim.1$par, optim.2$par, optim.3$par)
      models <- c('a', "a+bt", "a+bt^2", "a+bt+ct^2")
      cat("best model:", models[[best]], fill=T)
      cat("parameters:", params[[best]], fill=T)
   }
   
   print(bestCopula, quote = F)
   cat("Kendall's tau", fill=T)
   cat("Empirical:", corKendall(as.matrix(XZ))[2], fill=T)
   cat("Theorical:", BiCopPar2Tau(bestCopula$family, par=bestCopula$par), fill=T)
   
   if (gofTest == 'white') {
      cat("Goodness-of-fit test of White", fill=T)
      gf_test <- BiCopGofTest(uu[,1], uu[,2], obj = bestCopula, method='white')
      cat("Statistic:", round(gf_test$statistic, 4), ', ')
      cat("p-value:", gf_test$p.value)
   } else if (gofTest == 'kendall'){
      cat("Goodness-of-fit test of Kendall", fill=T)
      gf_test <- BiCopGofTest(uu[,1], uu[,2], obj = bestCopula, method='kendall')
      cat("Statistic:", round(gf_test$statistic, 4), ', ')
      cat("p-value:", gf_test$p.value)
   }
   
   return(bestCopula)
}


simulate_annual_prcp <- function(years, Z_fitted, W_fitted, D_fitted, params_GPD,
                                 RVM, best_copula.XZ, treshold, nsim=1e+3, seed=2021) {
   #' Fonction qui simule les précipitations annuelles totales pour les années
   #' spécifiées.
   #' 
   #' @param years (int): Vecteur donnant les années pour lesquelles on désire
   #'        faire les simulations.
   #' @param Z_fitted (list): Le résultat de la fonction fit_Z.
   #' @param W_fitted (list): Le résultat de la fonction fit_W.
   #' @param D_fitted (list): Le résultat de la fonction fit_duration
   #' @param params_GPD (float): Le vecteur des paramètres de la loi GPD, tels
   #'        que sortis par la fonction fit_GPD.
   #' @param RVM (RVineMatrix): Un objet de la classe RVineMatrix 
   #'        du package VineCopula.
   #' @param best_copula.XZ (list): Le résultat de la fonction copula_selection.XZ
   #'        pour modéliser la dépendance entre les variables V et Z.
   #' @param nsim (int): Le nombre de simulations désiré.
   #' @param seed (int): L'ancrage de simulations.
   eps <- 1e-8
   
   data_XZ <- data.extremes %>% 
      mutate(year = year(start_date)) %>%
      group_by(year) %>%
      summarise(annual_extremes = sum(u + excedence), .groups='drop') %>% 
      right_join(data.normales, by='year') %>%
      rename(X = annual_extremes,
             Z = total_prcp)
   
   T0 <- data.extremes$start_date %>% min()
   
   Fn_X <- ecdf(data_XZ$X)
   
   
   simul_annual_prcp <- function(.year) {
      qw <- W_fitted$quantiles
      qD <- D_fitted$quantiles
      qZ <- Z_fitted$quantiles
      
      Ti <- as.Date(paste0(as.character(.year), '-04-01'))
      Ti <- (Ti - T0) %>% as.numeric()
      t_max <- Ti + 91
      N_max <- 30
      
      u.XWD <- RVineSim(N_max, RVM)

      Wi <- list()
      i <- 0
      while(Ti <= t_max) {
         i <- i + 1
         Wi <- append(Wi, qw(u.XWD[i, 2], Ti))
         Ti <- Ti + as.numeric(tail(Wi, 1))
         Di <- qD(u.XWD[i, 3], Ti)
         Ti <- Ti + Di
      }
      Nt <- i - 1
      if (Nt <= 0) {
         X <- 0
      } else {
         # Total des précipitations extrêmes
         X <- sum(treshold + qGPD(u.XWD[1:Nt, 1], params_GPD[1], params_GPD[2]))
      }
      # Total des précipitations saisonnières
      u_Z <- BiCopCondSim(
         N = 1,
         cond.val = min(1-eps, max(eps, Fn_X(X))),
         cond.var = 1,
         obj = best_copula.XZ
      )
      Z <- qZ(u_Z, .year)
      # Total des précipitations de l'année
      St <- Z + X
      return(St)
   }
   
   f <- function(i) {
      S <- sapply(years, simul_annual_prcp)
      setTxtProgressBar(pb, i)
      return(S)
   }
   
   set.seed(seed)
   pb <- txtProgressBar(min=1, max=nsim, style=3)
   simuls <- sapply(1:nsim, f)
   return(simuls)
}


scatterplot_Chaoubi <- function(U, kernel_dist='norm') {
   #' Fonction implémentée par Ihsan Chaoubi, reprenant la méthode du package 
   #' VineCopula pour produire le graphique de la structure de dépendence 
   #' bivarié de chaque paire de v.a. d'un vecteur de pseudo-observations.
   #' 
   #' @param U (matrix): Matrice de pseudos observations où chaque colonne 
   #'        représente une variable aléatoire. Pour produire les 
   #'        pseudo-observations, voir la fonction pobs de la librairie copulas.
   #'        @param kernel_dist (string): Le triangle du bas de la matrice de 
   #'               graphique produite utilise une méthode par noyeau pour 
   #'               approximer la forme de la copule. Par défaut, un noyeau 
   #'               gaussien est utilisé pour effectuer le lissage de la 
   #'               distribution empirique. Les autres optios possibles sont
   #'               'uniforme' (Aucun lissage n'est effectué), 'exp' et 'flexp'.
   require(VineCopula)
   require(kdecopula)
   
   panel_up <- function(x, y) {
      points(
         x, y,
         cex = 0.3,
         col = rgb(red = 0, green = 0, blue = 0, alpha = 0.4),
         pch = 19
      )
   }
   
   panel_lower <-  function(x, y, cc=topo.colors(18), mar=kernel_dist, ...) {
      # set number of levels and if contour labels are drawn
      if (length(cc) == 1) {
         nlvls <- 11
         dl <- TRUE
      } else {
         nlvls <- length(cc)
         dl <- FALSE
      }
      # set levels according to margins
      if (mar %in% c("norm", "exp", "flexp")) {
         lvls <- seq(0.0, 0.25, length.out = nlvls)
      } else if (mar == "unif") {
         lvls <- seq(0.3, 4, length.out = nlvls)
      }
      # set default parameters
      pars <- list(u1 = x,
                   u2 = y,
                   size = 100,
                   levels = lvls,
                   margins = mar,
                   axes = FALSE,
                   drawlabels = dl)
      # get non-default parameters
      #pars <- modifyList(pars, list(...))
      pars <- modifyList(modifyList(pars, list(col = NULL)), list(col = cc))
      op <- par(usr = c(-3, 3, -3, 3), new = TRUE)
      # call BiCopMetaContour
      do.call(BiCopKDE, pars)
      r <- round(cor(x,y,method="kendall"), digits=3)
      text(1, -2.5, bquote( tau~" = "~.(r)), cex = 1.2,col="blue")
      on.exit(par(op))
   }

   pairs(
      as.copuladata(U),
      labels = c(
         expression(paste("u" ^ "(X)")),
         expression(paste("u" ^ "(W)")),
         expression(paste("u" ^ "(D)"))
      ),
      cex.labels = 1.5,
      diag.panel = NULL,
      upper.panel = panel_up ,
      lower.panel = panel_lower,
      label.pos = 0.5
   )
}
