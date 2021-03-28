library(ADGofTest)  # Statistique d'Anderson-Darling
library(copula)     # Librairie générale pour travailler avec des copules.
library(eva)        # Un autre papier plus récent sur la détermination du seuil.
library(lubridate)  # Travailler efficacement avec les dates
library(modifiedmk) # Version modifiée du test de Mann-Kenndall pour des données auto-corrélées.
library(qrmtools)   # Donne accès aux fonctions dGPD, pGPD, qGPD et rGPD.
library(threshr)    # Détermination du seuil de valeur extrême de façon automatique
library(tidyverse)  # Package pour être efficace en data science.
library(utils)      # Permet de faire des bars de progression.
library(VineCopula) # Librairie efficace pour faire la sélection de copule.
library(parallel)   # Permet de parréliser les simulations du modèle final.

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
            duration = max(DATE) - min(DATE) + 1,
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
      Fn <- ecdf(x)
      M100 <- mean(x)
      M101 <- mean(x * (1 - Fn(x)))
      xi <- -(M100 - 4*M101) / (2*M101 - M100)
      scale <- (1-xi) * M100
   }
   else if (method == 'ml') {
      nll <- function(param) {
         -sum(log(dGPD(x, param[1], param[2])))
      }
      param <- optim(par = fit_GPD(x, 'pwm'), fn=nll)$par
      xi <- param[1]
      scale <- param[2]
   }
   else if (method == 'mm') {
      mean_x <- mean(x)
      var_x <- var(x)
      xi <- 0.5 * (1 - mean_x^2/var_x)
      scale <- (1 - xi) * mean_x
   }
   else if (method == 'vwl') {
      # N <- trunc(length(x) * 0.75)
      
      func_to_minimize <- function(param) {
         N <- length(x)
         qm <- qGPD((1:N) / N, param[1], param[2])  # Quantiles théoriques
         xm <- sort(x)  # Quantiles empiriques
         return(mean((xm - qm) ^ 2))
      }
      
      param <- optim(par = fit_GPD(x, 'pwm'),
                     func_to_minimize)
      xi <- param[1]
      scale <- param[2]
   }
   else if (method == 'bayesian') {
      out <-  MCMC4Extremes::gpdp(x, 0, 500)
      xi <- out$postmedian[2]
      scale <- out$postmedian[1]
   }
   else stop("The specified method is not implemented in this function.")
   
   return(c(xi, scale))
}


fit_W <- function(data, show_qqplot=T) {
   #' Fonction qui paramétrise la loi des temps inter-occurences des périodes 
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
   #' @param show_qqplot (bool) : Est-ce que l'on désire afficher les qqplots 
   #'        pour valider l'adéquation de la loi paramétrée.
   #'        
   #' @return Une liste:
   #'        - density (func): Fonction de densité de probabilités paramétrée.
   #'        - distribution (func) Fonction de répartition des probabilités.
   #'        - quantiles (func): Fonction quantile paramétrée.
   #'        
   #'        Les paramètres utilisés pour ces deux fonctions sont:
   #'        - w (float): Les temps inter-occurences d'événements extrêmes.
   #'        - t0 (Date ou int): Date du début

   data_time <- data %>% select(W, start_date, duration)
   W <- data_time$W %>% unlist(use.names = F)
   t0 <- data_time$start_date + data_time$duration %>% unlist(use.names = F)
   years <- data_time$start_date %>% year() %>% unlist(use.names = F)
   n <- length(W)

   residual_times <- function(data_time) {
      t <- data_time$start_date %>% year() %>% unlist(use.names = F)
      min_year <- min(t)
      n_years <- max(t) - min_year + 1
      res <- numeric(n_years)
      t0 <- numeric(n_years)
      for (i in 1:n_years) {
         y <- min_year + i - 1
         annual_times <- data_time %>% filter(year(start_date) == y)
         t0[i] <- as.Date(paste0(as.character(y), '-04-01'))
         cum_time <- sum(annual_times$W, annual_times$duration)
         t0[i] <- t0[i] + cum_time
         res[i] <- 91 - cum_time
      }
      return(list("res"=res, "t0"=t0, 'years'=seq(min_year, max(t), by=1)))
   }
   res <- residual_times(data_time)
   
   mktest <- MannKendall(W)$sl[1]
   cat("Mann-Kendall's p-value:", mktest, fill = T)
   
   if (mktest <= 0.05) {
      
      mu <- function(t, s, param) {
         a <- param[1]
         b <- param[2]
         s <- as.numeric(s)
         (a + b * s) * t + 0.5 * b * t ^ 2
      }
      mu_inv <- function(u, s, param) {
         a <- param[1]
         b <- param[2]
         s <- as.numeric(s)
         (-(a+b*s) + sqrt((a+b*s)^2 + 2*b*u)) / b
      }
      
      dexp_ns <- function(t, s, param) {
         a <- param[1]
         b <- param[2]
         dmu <- function(t,s) {
            s <- as.numeric(s)
            a + b * (s + t)
         }
         exp(-mu(t, s, param)) * dmu(t, s)
      }
      exp_nll <- function(param) {
         nll <- -sum(log(dexp_ns(W, t0, param)))
         nll <- nll - sum(log(exp(-mu(res$res, res$t0, param))))
         return(nll)
      }
      
      gamma_nll <- function(param) {
         shape <- param[1]
         
         rate <- function(.year) {
            shape / (param[2] + param[3] * as.numeric(.year))
         }
         nll <- -sum(log(dgamma(W, shape, rate(years))))
         nll <- nll - sum(log(1-pgamma(res$res, shape, rate(res$years))))
         
         return(nll)
      }
      
      weibull_nll <- function(param) {
         shape <- param[1]
         
         scale <- function(.year) {
            (param[2] + param[3] * as.numeric(.year)) / gamma(1 + 1 / shape)
         }
         nll <- -sum(log(dweibull(W, shape, scale(years))))
         nll <- nll - sum(log(1-pweibull(res$res, shape, scale(res$years))))
         return(nll)
      }
      
      
      optim_exp <- optim(c(1 / mean(W), 0), exp_nll)
      optim_gamma <- optim(c(1, mean(W), 0), gamma_nll)
      optim_weibull <- optim(c(1, mean(W), 0), weibull_nll)
      
      aic_exp <- 2 * (2 + optim_exp$value)
      aic_gamma <-  2 * (3 + optim_gamma$value)
      aic_weibull <-  2 * (3 + optim_weibull$value)
      
   } 
   else {
      
      exp_nll <- function(param) {
         nll <- -sum(log(dexp(W, t0, param)))
         nll <- nll - sum(log(1 - pexp(res$res, res$t0, param)))
         return(nll)
      }
      
      gamma_nll <- function(param) {
         nll <- -sum(log(dgamma(W, param[1], param[2])))
         nll <- nll - sum(log(1 - pgamma(res$res, res$t0, param)))
         return(nll)
      }
      
      weibull_nll <- function(param) {
         nll <- -sum(log(dweibull(W, param[1], param[2])))
         nll <- nll - sum(log(1 - pweibull(res$res, res$t0, param)))
         return(nll)
      }
      
      optim_exp <- list()
      optim_exp$par <- 1 / mean(duration)
      optim_exp$value <- exp_nll(optim_exp$par)
      
      optim_gamma <- optim(c(1, optim_exp$par), gamma_nll)
      optim_weibull <- optim(c(1, 1 / optim_exp$par), weibull_nll)
      
      aic_exp <- 2 * (1 + optim_exp$value)
      aic_gamma <-  2 * (2 + optim_gamma$value)
      aic_weibull <-  2 * (2 + optim_weibull$value)
   }

   best <- which.min(c(aic_exp, aic_gamma, aic_weibull))
   cat('\n')
   print(cbind(aic_exp, aic_gamma, aic_weibull))
   cat('\n')

   distributions <- c("exp", "gamma", "weibull")
   params <- list(optim_exp$par,
                  optim_gamma$par,
                  optim_weibull$par)

   best_distribution <- distributions[best]
   best_param <- params[[best]] %>% unlist()

   cat("best distribution:", best_distribution, ', ')
   cat("Parameters:", best_param, fill = T)
   
   if (best_distribution == 'exp') {
      if (mktest <= 0.05) {
         dW <- function(w, t0) dexp_ns(w, t0, best_param)
         pW <- function(w, t0) 1 - exp(-mu(w, t0, best_param))
         qW <- function(u, t0) mu_inv(qexp(u), t0, best_param)
      } else {
         dW <- function(w, t0) dexp(w, best_param[1], best_param[2])
         pW <- function(w, t0) pexp(w, best_param[1], best_param[2])
         qW <- function(u, t0) qexp(u, best_param[1], best_param[2])
      }
      
   } else if (best_distribution == 'gamma') {
      if (mktest <= 0.05) {
         shape <- best_param[1]
         rate <- function(t0) {
            .year <- year(t0)
            shape / (best_param[2] + best_param[3] * .year)
         }
         dW <- function(w, t0) dgamma(w, shape, rate(t0))
         pW <- function(w, t0) pgamma(w, shape, rate(t0))
         qW <- function(u, t0) qgamma(u, shape, rate(t0))
      } else {
         dW <- function(w, t0) dgamma(w, best_param[1], best_param[2])
         pW <- function(w, t0) pgamma(w, best_param[1], best_param[2])
         qW <- function(u, t0) qgamma(u, best_param[1], best_param[2])
      }
      
   } else { # Weibull
      if (mktest <= 0.05) {
         shape <- best_param[1]
         
         scale <- function(t0) {
            .year <- year(t0)
            (best_param[2] + best_param[3] * .year) / gamma(1 + 1/shape)
         }
         dW <- function(w, t0) dweibull(w, shape, scale(t0))
         pW <- function(w, t0) pweibull(w, shape, scale(t0))
         qW <- function(u, t0) qweibull(u, shape, scale(t0))
      } else {
         dW <- function(w, t0) dweibull(w, best_param[1], best_param[2])
         pW <- function(w, t0) pweibull(w, best_param[1], best_param[2])
         qW <- function(u, t0) qweibull(u, best_param[1], best_param[2])
      }
   }
   
   if (show_qqplot) {
      p <- pW(91, t0)
      qqplot(W, qW(p * ecdf(W)(W), t0),
             xlab = 'Empirical observations',
             ylab = 'Theorical quantiles')
      abline(0, 1, col = 'blue')
   }
   
   return(list(
      "density" = dW,
      "distribution" = pW,
      "quantiles" = qW,
      'params' = best_param
   ))
}


fit_duration <- function(data, show_qqplot=T) {
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
   #' @param show_qqplot (bool) : Est-ce que l'on désire afficher les qqplots 
   #'        pour valider l'adéquation de la loi paramétrée.
   #'        
   #' @return Une liste:
   #'        - density (func): Fonction de densité de probabilités paramétrée.
   #'        - distribution (func) Fonction de répartition des probabilités.
   #'        - quantiles (func): Fonction quantile paramétrée.
   #'        
   #'        Les paramètres utilisés pour ces deux fonctions sont:
   #'        - x (float): La durée des événements extrêmes.
   #'        - t0 (Date ou int): Date du début

   data_time <- data %>% select(W, start_date, duration)
   duration <- data_time$duration %>% as.numeric() %>% unlist(use.names = F)
   t0 <- data_time$start_date %>% unlist(use.names = F)
   years <- data_time$start_date %>% year() %>% unlist(use.names = F)

   mktest <- MannKendall(duration)$sl[1]
   cat("Mann-Kendall's p-value:", mktest, fill = T)
   
   if (mktest <= 0.05) {
      
      mu <- function(t, s, param) {
         a <- param[1]
         b <- param[2]
         s <- as.numeric(s)
         (a + b * s) * t + 0.5 * b * t ^ 2
      }
      mu_inv <- function(u, s, param) {
         a <- param[1]
         b <- param[2]
         s <- as.numeric(s)
         (-(a+b*s) + sqrt((a+b*s)^2 + 2*b*u)) / b
      }
      
      dexp_ns <- function(t, s, param) {
         a <- param[1]
         b <- param[2]
         dmu <- function(t,s) {
            s <- as.numeric(s)
            a + b * (s + t)
         }
         exp(-mu(t, s, param)) * dmu(t, s)
      }
      exp_nll <- function(param) {
         nll <- -sum(log(dexp_ns(duration, t0, param)))
         return(nll)
      }
      exp_nll2 <- function(param) {
         mu2 <- function(.year) {
            (param[1] + param[2] * as.numeric(.year))
         }
         nll <- -sum(log(dexp(duration, 1/mu2(years))))
         return(nll)
      }
      
      gamma_nll <- function(param) {
         shape <- param[1]
         
         rate <- function(.year) {
            shape / (param[2] + param[3] * as.numeric(.year))
         }
         nll <- -sum(log(dgamma(duration, shape, rate(years))))
         return(nll)
      }
      
      weibull_nll <- function(param) {
         shape <- param[1]
         
         scale <- function(.year) {
            (param[2] + param[3] * as.numeric(.year)) / gamma(1 + 1 / shape)
         }
         nll <- -sum(log(dweibull(duration, shape, scale(years))))
         return(nll)
      }
      
      
      optim_exp <- optim(c(1/mean(duration), 0), exp_nll)
      optim_exp2 <- optim(c(mean(duration), 0), exp_nll2)
      optim_gamma <- optim(c(1, mean(duration), 0), gamma_nll)
      optim_weibull <- optim(c(1, mean(duration), 0), weibull_nll)
      
      aic_exp <- 2 * (2 + optim_exp$value)
      aic_exp2 <- 2 * (2 + optim_exp2$value)
      aic_gamma <-  2 * (3 + optim_gamma$value)
      aic_weibull <-  2 * (3 + optim_weibull$value)
      
   } 
   else {
      
      exp_nll <- function(param) {
         nll <- -sum(log(dexp(duration, param)))
         return(nll)
      }
      
      gamma_nll <- function(param) {
         nll <- -sum(log(dgamma(duration, param[1], param[2])))
         return(nll)
      }
      
      weibull_nll <- function(param) {
         nll <- -sum(log(dweibull(duration, param[1], param[2])))
         return(nll)
      }
      
      optim_exp <- list()
      optim_exp$par <- 1 / mean(duration)
      optim_exp$value <- exp_nll(optim_exp$par)
      optim_exp2 <- list("par"=0, "value"=1e+10)
      
      optim_gamma <- optim(c(1, optim_exp$par), gamma_nll)
      optim_weibull <- optim(c(1, 1 / optim_exp$par), weibull_nll)
      
      aic_exp <- 2 * (1 + optim_exp$value)
      aic_exp2 <- 2 * (1 + optim_exp2$value)
      aic_gamma <-  2 * (2 + optim_gamma$value)
      aic_weibull <-  2 * (2 + optim_weibull$value)
   }
   
   best <- which.min(c(aic_exp, aic_exp2, aic_gamma, aic_weibull))
   cat('\n')
   print(cbind(aic_exp, aic_exp2, aic_gamma, aic_weibull))
   cat('\n')
   
   distributions <- c("exp", "exp", "gamma", "weibull")
   params <- list(optim_exp$par,
                  optim_exp2$par,
                  optim_gamma$par,
                  optim_weibull$par)
   
   best_distribution <- distributions[best]
   best_param <- params[[best]] %>% unlist()
   
   cat("best distribution:", best_distribution, ', ')
   cat("Parameters:", best_param, fill = T)
   
   if (best_distribution == 'exp') {
      if (mktest <= 0.05) {
         dD <- function(d, t0) dexp_ns(d, t0, best_param)
         pD <- function(d, t0) 1 - exp(-mu(d, t0, best_param))
         qD <- function(u, t0) mu_inv(qexp(u), t0, best_param)
      } else {
         dD <- function(d, t0) dexp(d, best_param[1], best_param[2])
         pD <- function(d, t0) pexp(d, best_param[1], best_param[2])
         qD <- function(u, t0) qexp(u, best_param[1], best_param[2])
      }
      
   } 
   else if (best_distribution == 'gamma') {
      if (mktest <= 0.05) {
         shape <- best_param[1]
         rate <- function(t0) {
            .year <- year(t0)
            shape / (best_param[2] + best_param[3] * .year)
         }
         dD <- function(d, t0) dgamma(d, shape, rate(t0))
         pD <- function(d, t0) pgamma(d, shape, rate(t0))
         qD <- function(u, t0) qgamma(u, shape, rate(t0))
      } 
      else {
         dD <- function(d, t0) dgamma(d, best_param[1], best_param[2])
         pD <- function(d, t0) pgamma(d, best_param[1], best_param[2])
         qD <- function(u, t0) qgamma(u, best_param[1], best_param[2])
      }
      
   } 
   else { # Weibull
      if (mktest <= 0.05) {
         shape <- best_param[1]
         
         scale <- function(t0) {
            .year <- year(t0)
            (best_param[2] + best_param[3] * .year) / gamma(1 + 1/shape)
         }
         dD <- function(d, t0) dweibull(d, shape, scale(t0))
         pD <- function(d, t0) pweibull(d, shape, scale(t0))
         qD <- function(u, t0) qweibull(u, shape, scale(t0))
      } else {
         dD <- function(d, t0) dweibull(d, best_param[1], best_param[2])
         pD <- function(d, t0) pweibull(d, best_param[1], best_param[2])
         qD <- function(u, t0) qweibull(u, best_param[1], best_param[2])
      }
   }
   
   if (show_qqplot) {
      qqplot(duration, qD(ecdf(duration)(duration), t0),
             xlab = 'Empirical observations',
             ylab = 'Theorical quantiles')
      abline(0, 1, col = 'blue')
   }
   
   return(list(
      "density" = dD,
      "distribution" = pD,
      "quantiles" = qD,
      'params' = best_param
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


copula_selection <- function(data, variable='W', dep_trend=F, family_set=NA,
                             gofTest='None', seed=42) {
   #' Fonction qui sélectionne la copule la plus adéquate parmi toutes les 
   #' copules programmées dans la librairie VineCopula.
   #' 
   #' @param data (tibble): DataFrame tibble contenant les colonnes suivantes:
   #'        - excedence : L'excédent de seuil pour les précipitations.
   #'        - W : les temps inter-périodes de précipitations extrêmes.
   #' @param variable (str): Le nom de la colonne utilisée pour vérifier une 
   #'        tendance dans la dépendance avec la variable d'excédent de seuil.
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
   
   xw <- data %>% select(excedence, all_of(variable)) %>%
      na.omit()
   xw[variable] <- xw[variable] %>% unlist() %>% as.numeric()
   years <- data$start_date %>% year()
   white_noise <- rnorm(2 * nrow(xw), mean = 0, sd=1e-8) %>% matrix(ncol = 2)
   xw <- xw + white_noise
   uu <- pobs(xw)
   
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
   cat("Empirical:", corKendall(as.matrix(xw))[2], fill=T)
   cat("Theorical:", BiCopPar2Tau(
      bestCopula$family, par=bestCopula$par, par2=bestCopula$par2),
      fill=T)
   
   if (grepl('BB', bestCopula$familyname)) {
      gofTest = 'kendall'
   }
   
   if (gofTest == 'white') {
      cat("Goodness-of-fit test of White", fill=T)
      gf_test <- BiCopGofTest(uu[,1], uu[,2], obj = bestCopula, method='white')
      cat("Statistic:", round(gf_test$statistic, 4), ', ')
      cat("p-value:", gf_test$p.value)
   } else if (gofTest == 'kendall'){
      cat("Goodness-of-fit test of Kendall", fill=T)
      gf_test <- BiCopGofTest(uu[,1], uu[,2], obj = bestCopula, method='kendall')
      cat("Cramer-von Mises statistic:", round(gf_test$statistic.CvM, 4), ', ')
      cat("p-value:", gf_test$p.value.CvM, fill=T)
      cat("Kolmogorv-Smirnov statistic:", round(gf_test$statistic.KS, 4), ', ')
      cat("p-value:", gf_test$p.value.KS, fill=T)
   }
   
   return(bestCopula)
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
   white_noise <- rnorm(2 * nrow(XZ), mean = 0, sd=1e-8) %>% matrix(ncol = 2)
   XZ <- XZ + white_noise
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

simulate_annual_prcp <- function(years, W_fitted, D_fitted, params_GPD,
                                 best_copula.W, best_copula.D, best_copula.XZ,
                                 Z_lm, treshold, nsim=1e+3, seed=42) {
   #' Fonction qui simule les précipitations annuelles totales pour les années
   #' spécifiées.
   #' 
   #' @param years (int): Vecteur donnant les années pour lesquelles on désire
   #'        faire les simulations.
   #' @param W_fitted (list): Le résultat de la fonction fit_W.
   #' @param D_fitted (list): Le résultat de la fonction fit_duration
   #' @param params_GPD (float): Le vecteur des paramètres de la loi GPD, tels
   #'        que sortis par la fonction fit_GPD.
   #' @param best_copula.W (list): Le résultat de la fonction copula_selection 
   #'        pour modéliser la dépendance entre les variables X-u|X>u et W.
   #' @param best_copula.D (list): Le résultat de la fonction copula_selection 
   #'        pour modéliser la dépendance entre les variables X-u|X>u et D.
   #' @param Z_lm (lm): Modèle linéaire entraîné sur les précipitations 
   #'        saisonnières
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
   Fn_X <- ecdf(data_XZ$X)
   
   
   simul_annual_prcp <- function(.year) {
      qw <- W_fitted$quantiles
      qD <- D_fitted$quantiles
      
      Ti <- as.Date(paste0(as.character(.year), '-04-01'))
      t_max <- as.Date(paste0(as.character(.year), '-06-30'))
      N_max <- 30
      
      u.WX <- BiCopSim(N_max, obj = best_copula.W)
      u.DX <- BiCopCondSim(
         N = N_max,
         cond.val = u.WX[, 1],
         cond.var = 1,
         obj = best_copula.D
         )
      Wi <- list()
      i <- 0
      while(Ti <= t_max) {
         i <- i + 1
         Wi <- append(Wi, qw(u.WX[i, 2], Ti))
         Ti <- Ti + as.numeric(tail(Wi, 1))
         Di <- qD(u.DX[i], Ti)
         Ti <- Ti + Di
      }
      Nt <- i - 1
      if (Nt <= 0) {
         X <- 0
      } else {
         # Total des précipitations extrêmes
         X <- sum(treshold + qGPD(u.WX[1:Nt, 1], params_GPD[1], params_GPD[2]))
      }
      # Total des précipitations saisonnières
      Z_mu <- predict.lm(Z_lm, as.data.frame(list("year"=.year)))
      Z_sd <- summary(Z_lm)$sigma
      u_Z <- BiCopCondSim(
         N = 1,
         cond.val = min(1-eps, max(eps, Fn_X(X))),
         cond.var = 1,
         obj = best_copula.XZ
      )
      Z <- qnorm(u_Z, Z_mu, Z_sd)
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


fit_model <- function(data, nsim=1e+2) {
   #' Pipeline qui prend en entrée le dataset et qui retourne les fonctions
   #' nécessaires pour faire une analyse complête du risque et pour vérifier
   #' l'adéquation du modèle.
   #' 
   #' @param data (tibble): Les données de précipitations quotidiennes pour 
   #'        plusieurs années. Ces données doivent inclure des observations pour
   #'        tous les jours du printemps, soit du 1er avril au 30 juin de 
   #'        chacune des années. Les champs obligatoires sont:
   #'        - DATE (Date): La date de l'observation en format as.Date.
   #'        - PRCP (float): Les précipitations (en mm de pluie)
   #' @param nsim (int): Nombre de simulations pour étudier les prédictions du
   #'        modèle.
   #' 
   #' @return une liste de fonctions utiles à l'étude complète du risque:
   #'        - cdf_S: Fonction de répartition conditionnellement à l'année
   #'        - Var_S: VaR conditionnelle selon l'année
   #'        - TVaR_S: TVaR conditionnelle selon l'année
   #'        - global_quantiles: VaR globale (toutes années confondues)
   #'        - global_cdf: Fonction de répartition globale (toutes années 
   #'          confondues).
   
   cat('preprocessing:', fill=T)
   data <- preprocess_data(
      data, 
      NA.neighbour_days = 2,
      C.separeting_days = 1
   )
   
   verify_serial_dependence(data, max_lag=1)
   cat('--------------------------------------------------------', fill = T)
   
   # Sélection du seuil de valeurs extrêmes
   set.seed(42)
   p_vec <- rnorm(100, mean=0.86, sd=0.04) %>% sort()
   u_vec <- quantile(data$prcp_agg, probs = p_vec) 
   u <- treshold_selection(data, u_vec, method='bader', gf_test = 'ad')
   
   # Séparation des données en valeurs extrêmes et valeurs saisonnières
   data.extremes <- calculate_excedences(data, u, calculate_W=T)

   # Entraînement d'un modèle linéaire pour prédire les précipitations 
   # saisonnières.
   data.normales <- data %>%
      filter(prcp_agg < u) %>%
      select(prcp_agg, start_date) %>%
      group_by(year(start_date)) %>%
      summarise(sum(prcp_agg), .groups='drop') %>%
      rename('year' = `year(start_date)`, total_prcp = `sum(prcp_agg)`)

   Z_lm <- lm(total_prcp ~ year, data=data.normales)
   car::qqPlot(Z_lm, lwd=0.3, id=F, main="QQ-plot of seasonal prcp")
   
   # Entraînement des valeurs extrêmes
   cat("Mann-Kendall test for extremes:")
   print(MannKendall(data.extremes$excedence))
   params_GPD <- fit_GPD(data.extremes$excedence, method = 'mm')
   
   # Entraînement de la loi pour les temps inter-événements extrêmes
   cat("Mann-Kendall test for inter-occurence times:")
   print(MannKendall(data.extremes$W))
   W_fitted <- fit_W(data.extremes, show_qqplot=T)
   
   # Modélisation de la dépendance
   data.extremes %>% select(excedence, W) %>% na.omit() %>% 
      calculate_correlation_WX(alternative = 'greater')
   
   eval_dependence_trend(data.extremes, stride=150)
   
   best_copula <- copula_selection(
      data.extremes,
      dep_trend=F,
      plot_copula=T,
      gofTest = 'white')
   
   years <- data$start_date %>% year() %>% unique()
   simulations <- simulate_annual_prcp(years, W_fitted, best_copula, params_GPD,
                                       Z_lm, treshold = u, nsim = nsim)
   
   # Conditionnal functions
   cdf_S <- function(x, year) {
      ind <- which(years == year)
      ecdf(simulations[ind, ])(x)
   }
   VaR_S <- function(p, year) {
      ind <- which(years == year)
      n <- ncol(simulations)
      sort(simulations[ind, ])[ceiling(p * n)]
   }
   TVaR_S <- function(p, year) {
      ind <- which(years == year)
      n <- ncol(simulations)
      sorted_obs <- sort(simulations[ind, ])
      VaR <- sorted_obs[ceiling(p * n)]
      TVaR <- mean(sorted_obs[sorted_obs > VaR])
      return(TVaR)
   }
   
   # Global functions
   global_quantiles <- function(p) {
      simulations <- as.vector(simulations)
      n <- ncol(simulations)
      sort(simulations)[ceiling(p * n)]
   }
   globap_cdf <- function(x) {
      vec_sim <- as.vector(simulations)
      ecdf(vec_sim)(x)
   }
   
   return(
      list(
         "distribution" = cdf_S,
         "quantiles" = VaR_S,
         "TVaR" = TVaR_S,
         "global.quantiles" = global_quantiles,
         "global.cdf" = globap_cdf
      )
   )
   
}
