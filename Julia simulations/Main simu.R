# This simulation format should be preferred over the previous one.
rm(list = ls())


# Simulation:
# packages loading ---------------------------------
list.of.packages <- c("tidyverse", "MASS", "stats", "matrixcalc", "stargazer", "foreach", "AER", "doParallel", "ivmodel", "np", "latex2exp", "SteinIV", "magrittr", "weights", "hdm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("tidyverse")
library("MASS")
library("stargazer")
library("foreach")
library("doParallel")
library("AER")
library("matrixcalc")
library("stats")
library("ivmodel")
library("np")
library("latex2exp")
library("SteinIV")
library("magrittr")
library("DEoptim")
library("hdm")     # package for high dimensional models, with the estimator of Belloni, Chen, Chernozukov etc (2012)

# source("Lee_and_Shin.R")   # Donald and Newey (2001) estimator and Lee and Shin (2018) estimator
source("auxiliary functions.R")
source("DGP parameters.R")
# source("P2SLS.R")
# Setup for parallel looping ----------------------------------

parallel:::detectCores()   # detect the number of cores in the CPU of da computa
clust <- makeCluster(parallel:::detectCores()-1)
registerDoParallel(clust )

# Main loop -----------------------------------------
ptm <- proc.time()  
listu <- foreach(s = 1:nsim, .combine = rbind, .packages = c("AER", "MASS", "doParallel", "tidyverse", "np", "hdm"))%:%
  foreach(i = 1:length(VCOV_list) , .combine = rbind)%:%
  # foreach(f = fs_fun, .combine = rbind)%:%
  # foreach( t = type, .combine = rbind)%:%
  foreach( n = N, .combine = rbind) %do%  
  {
    cat(paste("\n","Starting iteration",s,"\n"),
        file="log.txt", append=FALSE)
    
    VCOV <- VCOV_list[[i]]
    # data <- BWdgp( n, beta = param, varz = varz, VCOV, con, polyorder = d )
    data <- dgp(n, parameters = param, varz = 1, VCOV, g = fs_fun, DGP = "normal")
    y <- data$y
    x <- data$x
    z <- data$z1
    
    varu <-  VCOV[1,1]
    varv <-  VCOV[2,2]
    covuv <- VCOV[1,2] 
    varnames <- "x"
    
    
    ols <- lm(y ~ -1 + x )
    betaols <- t(ols$coefficients)
    betaols_var <- rbind(diag(vcov(ols))) 
    tsls <- ivreg(y ~ -1 + x | z )
    beta2sls <- t(tsls$coefficients)
    beta2sls_var <- rbind(diag(vcov(tsls))) 
    
    # map the bandwidth with estimators
    mybeta_ll <- bandrange %>% map_dfr(myest_np, expl = y, endo = x, exo = NULL, instru = z, type = "ll", intercept = FALSE)  # faster than myest
    mybeta_lc <- bandrange %>% map_dfr(myest_np, expl = y, endo = x, exo = NULL, instru = z, type = "lc", intercept = FALSE)
    
    bandwidth <- bandrange
    
    ll <- dplyr::select(mybeta_ll, !ends_with("var"))   # the estimated coefficients
    ll_var <- dplyr::select(mybeta_ll, ends_with("var")) # The estimated variances
    
    # ll_var <- tibble(ll_var = mybeta_ll$varhat) 
    
    lc <- dplyr::select(mybeta_lc, !ends_with("var"))  
    lc_var <- dplyr::select(mybeta_lc, ends_with("var")) # The estimated variances
    
    np <- np_2SLS(Y = y, endo = x, exo = NULL, Z = z, type = "ll", intercept = FALSE)
    np_bhat <- t(np$beta_hat)
    np_varhat <- np$varhat_homoskedasticity#"endo", "endo"]
    opti_h_ll <- np$opti_h
    
    WMD <- wmd( x = x , y = y , z = z, exo = NULL, intercept = FALSE )
    WMDF <- t(WMD$WMDF)
    WMD <- t(WMD$WMD)
    
    iv <- ivmodel::ivmodel(Y = y, D = x , Z = z, beta0 = 1)
    LIML <- as.numeric(ivmodel::LIML(iv, beta0 = 1, alpha = 0.05, heteroSE = FALSE, clusterID = NULL)$point.est)  

    # P2SLS
    # P2SLS <- P2SLS(Y = as.matrix(y), X = dat$dmexfb_alt, Z = dat$dins, K = K, W = as.matrix(W), option = "MSE", weights = dat$popweight, int = TRUE)
  
    # The estimators from Lee and Shin (2019) & Donald and Newey (2001)
    psiz <- poly(data$z1, degree = 10, raw = F, simple = TRUE)
    # # The estimator from Donald and Newey
    DN <- tsls_DN(y = y, x.end = x, z.excl = psiz, method.pre = "one.step", z.pre = psiz)
    DN <- t(DN$bt.hat)
    # # The estimator from Lee and Shin (2018)
    CSA <- tsls_CSA(y = y, x.end = x, z.excl = psiz, method.pre = "one.step", z.pre = psiz)$bt.hat
    CSA <- t(CSA)
    
    # Add rlassoIV by Belloni et al. (2012)
    lassoiv <- rlassoIV(x = psiz, d = x, z = psiz, y = y, select.Z = TRUE)  # Not sure it is correct. Double check
    lasso <- t(lassoiv$coefficients)
    lasso_se <- t(lassoiv$se)
    
    # Change the names of all the estimates
    colnames(betaols) <- paste("ols", varnames, sep = "_")
    colnames(beta2sls) <- paste("tsls", varnames, sep = "_")
    
    colnames(betaols_var) <- paste("ols", varnames, "var", sep = "_")
    colnames(beta2sls_var) <- paste("tsls", varnames, "var", sep = "_")
    
    colnames(ll) <- paste("ll", varnames, sep = "_")
    colnames(lc) <- paste("lc", varnames, sep = "_")
    
    colnames(ll_var) <- paste("ll", varnames, "var", sep = "_")
    colnames(lc_var) <- paste("lc", varnames, "var", sep = "_")    
    
    colnames(np_bhat) <- paste("np", varnames, sep = "_")
    colnames(np_varhat) <- paste("np", varnames, "var", sep = "_")
    
    colnames(WMD) <- paste("WMD", varnames, sep = "_")
    colnames(WMDF) <- paste("WMDF", varnames, sep = "_")    
    
    colnames(DN) <- paste("DN", varnames, sep = "_")
    colnames(CSA) <- paste("CSA", varnames, sep = "_")    
    
    colnames(lasso) <-  paste("Lasso", varnames, sep = "_") 
    colnames(lasso_se) <-  paste("Lasso", varnames, "var", sep = "_")
    
    
    # loo <- loo_choice(x = x, y = y, z = z, h = band, regtype = "ll")
    # bhat_loo <- loo$betahat.betahat
    # bhat_loo_var <- loo$betahat.varhat
    
    # loo <- loo_opti(x, y, z, regtype = "ll")
    # bhat_loo <- loo$betahat.betahat
    # bhat_loo_var <- loo$betahat.varhat
    # loo_h <- loo$opti_h
    # 
    # h_aic <- min_AIC_iv(expl = y, endo = x, instru = z, lower = 0.001, upper = 10, type = "ll", loo = FALSE)
    # aic <-  myest_np(expl = y, endo = x, instru = z, h = h_aic, type = "ll") 
    # bhat_aic <- aic$betahat
    # bhat_aic_var <- aic$varhat
    # psiz <- poly(z, degree = 10, raw = F, simple = TRUE)
    # DN <- tsls_DN(y = y, x.end = x, z.excl = psiz,  method.pre = "one.step", z.pre = psiz)$bt.hat 
    # CSA <- tsls_CSA(y = y, x.end = x, z.excl = psiz, method.pre = "one.step", z.pre = psiz)$bt.hat  
    
    # crit <- bandrange %>% map_dbl(iv_AIC, expl = y, endo = x, instru = z, type = t, loo = FALSE)
    # h_opti <- bandrange[which.min(crit)]
    # AIC_est <- myest(expl = y, endo = x, instru = z, h = h_opti, type = t, cv = FALSE)
    # 
    # # crit_loo <- bandrange %>% map_dbl(iv_AIC, expl = y, endo = x, instru = z, type = t, loo = TRUE)
    # # h_opti_loo <- bandrange[which.min(crit_loo)]
    # # AIC_est_loo <- myest(expl = y, endo = x, instru = z, h = h_opti_loo, type = t, cv = TRUE)
    # 
    # iv_AIC(expl = y, endo = x, instru = z, bandwidth = 0.01, type = "ll", loo = FALSE)
    
    # sol_ll <- min_AIC_iv(expl = y, endo = x, instru = z, lower = n^(-1), upper = 100, type = "ll", loo = FALSE)
    # hstar_ll <- sol_ll$minimizer
    # minimum_ll <- sol_ll$minimum    
    # AIC_opti_est_ll <- myest_np(expl = y, endo = x, instru = z, h = hstar_ll, type = "ll")$betahat
    # 
    # sol_lc <- min_AIC_iv(expl = y, endo = x, instru = z, lower = n^(-1), upper = 100, type = "lc", loo = FALSE)
    # hstar_lc <- sol_lc$minimizer
    # minimum_lc <- sol_lc$minimum    
    # AIC_opti_est_lc <- myest_np(expl = y, endo = x, instru = z, h = hstar_lc, type = "lc")$betahat
    
    print(s)
    return(data.frame(simu = s, sample = n, varu = varu, varv = varv, covuv = covuv, 
                      fs_fun = fs_fun,
                      betaols,
                      betaols_var,
                      beta2sls,
                      beta2sls_var,
                      bandwidth = bandrange,
                      ll,
                      ll_var,
                      lc,
                      lc_var,
                      np_bhat,
                      np_varhat,
                      WMD,
                      WMDF,
                      LIML,   #,
                      CSA,
                      DN,
                      lasso,
                      lasso_se
                      
    ))
    print(s)
  }
proc.time() - ptm 
save.image() 
message <- geterrmessage()
cat(paste("\n","Last error", message,"\n"),
    file="errors.txt", append=FALSE)
name_file <- paste("kernel 2SLS, ", fs_fun, " design", ".RData", sep = "")
# name_file <- paste("kernel 2SLS simulation, ", "all fs_fun", ".RData", sep = "")
save(listu , file = name_file)

# Results ----------------------------

# The most recent results (nsim = 10,000 and all designs)
# load("kernel 2SLS simulation, all fs_fun.RData")
# 
# # The most recent results with separate designs
# load("kernel 2SLS, A&L DGP.RData")
# load("kernel 2SLS, linear DGP.RData")
# load("kernel 2SLS, log DGP.RData")
# 
# # The results in the respective folders
# load("Simulation results/linear design/Alternative kernel 2SLS consistency, linear DGP.RData")
# load("Simulation results/A&L design/Alternative kernel 2SLS consistency, A&L DGP.RData")
# load("Simulation results/doppler + dnorm design/Alternative kernel 2SLS consistency, dnorm, doppler DGP.RData")
# load("Simulation results/binary x design/Alternative kernel 2SLS consistency, binary_x DGP.RData")
# load("Simulation results/binary x design/Alternative kernel 2SLS consistency, binary_x DGP (newer).RData")   # for power graphs
# load("Simulation results/log design/Alternative kernel 2SLS consistency, log DGP.RData")
# 
# # The results from the Julia simulations
# listu <- read_csv("first results (most recent).csv")
# listu <- read_csv("first results (most recent, Dec 6).csv")
listu <- read_csv("first results (Dec 13).csv")



MC_b <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "bias" ) 

MC_v <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "variance" )

MC_MSE <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "MSE" )

MC_medAD <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "medAD" )

MC_medbias <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
  group_by( sample, covuv, fs_fun, bandwidth ) %>% 
  summarise_all (MSE, parameter = param, type = "medbias" )

# MC_b <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
#   group_by( sample, varu, varv, covuv, fs_fun, bandwidth ) %>% 
#   summarise_all (MSE, parameter = param, type = "bias" ) 
# 
# MC_v <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
#   group_by( sample, varu, varv, covuv, fs_fun, bandwidth ) %>% 
#   summarise_all (MSE, parameter = param, type = "variance" )
# 
# MC_MSE <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
#   group_by( sample, varu, varv, covuv, fs_fun, bandwidth ) %>% 
#   summarise_all (MSE, parameter = param, type = "MSE" )
# 
# MC_medAD <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
#   group_by( sample, varu, varv, covuv, fs_fun, bandwidth ) %>% 
#   summarise_all (MSE, parameter = param, type = "medAD" )
# 
# MC_medbias <- listu[, -1] %>% dplyr::select(!ends_with("var")) %>%
#   group_by( sample, varu, varv, covuv, fs_fun, bandwidth ) %>% 
#   summarise_all (MSE, parameter = param, type = "medbias" )

# Graphs ---------------------
# MC_b <- MC_b %>% dplyr::select(!contains("W1") & !contains("W2"))
# MC_v <- MC_v %>% dplyr::select(!contains("W1") & !contains("W2"))
# MC_MSE <- MC_MSE %>% dplyr::select(!contains("W1") & !contains("W2"))


# Not needed if simulations come from Julia
MC_b <- MC_b %>% rename(ols = ols_x,
                        tsls = tsls_x,
                        ll = ll_x,
                        lc = lc_x,
                        np_tsls = np_x,
                        WMD = WMD_x,
                        WMDF = WMDF_x,
                        DN = DN_x,
                        CSA = CSA_x
                        #,
                        # ptsls = ptsls     
)
MC_v <- MC_v %>% rename(ols = ols_x,
                        tsls = tsls_x,
                        ll = ll_x,
                        lc = lc_x,
                        np_tsls = np_x,
                        WMD = WMD_x,
                        WMDF = WMDF_x,
                        DN = DN_x,
                        CSA = CSA_x#,
                        # ptsls = ptsls     
)
MC_MSE <- MC_MSE %>% rename(ols = ols_x,
                            tsls = tsls_x,
                            ll = ll_x,
                            lc = lc_x,
                            np_tsls = np_x,
                            WMD = WMD_x,
                            WMDF = WMDF_x,
                            DN = DN_x,
                            CSA = CSA_x#,
                            # ptsls = ptsls     
)

MC <- list(MC_b, MC_v, MC_MSE)
# If looking at different bias measures
# MC <- list(MC_b, MC_medbias, MC_medAD)
# unique_designs <- unique(MC_MSE[, 1:5])    # As the same designs show up multiple times due to the different values of the bandwidth, I report the unique designs
unique_designs <- unique(MC_MSE[, 1:3])    # The Julia simulation does not have varu and varv

1:nrow(unique_designs) %>% map(MC_graph, MC = MC, exo = FALSE)

# Relative bias, variance etc... -------------------

relative_b <- MC_b %>% dplyr::select(!ends_with("var") & !bandwidth & !contains("_W") & !ll_x & !lc_x) %>%
  group_by( sample, varu, varv, covuv, fs_fun) %>%
  summarise_all(mean) %>%   # all the numbers are constant for a same design, so the mean does not change anything. It allows to regroup everything to avoid repeating design because of previous bandwidths values
  mutate(ols = ols_x/tsls_x,
         liml = LIML/tsls_x,
         WMDF = WMDF/tsls_x,
         np_x = np_x/tsls_x,
         # ll_x = ll_x/tsls_x,
         # lc_x = lc_x/tsls_x,
         tsls_x = tsls_x/tsls_x  ) 


relative_v <- MC_v %>% dplyr::select(!ends_with("var") & !bandwidth & !contains("_W")) %>%
  group_by( sample, varu, varv, covuv, fs_fun) %>%
  summarise_all(mean) %>%   # all the numbers are constant for a same design, so the mean does not change anything. It allows to regroup everything to avoid repeating design because of previous bandwidths values
  mutate(ols = ols_x/tsls_x,
         liml = LIML/tsls_x,
         WMDF = WMDF/tsls_x,
         np_x = np_x/tsls_x,
         ll_x = ll_x/tsls_x,
         lc_x = lc_x/tsls_x,
         tsls_x = tsls_x/tsls_x  ) 


relative_MSE <- MC_MSE %>% dplyr::select(!ends_with("var") & !bandwidth & !contains("_W")) %>%
  group_by( sample, varu, varv, covuv, fs_fun) %>%
  summarise_all(mean) %>%   # all the numbers are constant for a same design, so the mean does not change anything. It allows to regroup everything to avoid repeating design because of previous bandwidths values
  mutate(ols = ols_x/tsls_x,
         liml = LIML/tsls_x,
         WMDF = WMDF/tsls_x,
         np_x = np_x/tsls_x,
         ll_x = ll_x/tsls_x,
         lc_x = lc_x/tsls_x,
         tsls_x = tsls_x/tsls_x  ) 

relative_medbias <- MC_medbias %>% dplyr::select(!ends_with("var") & !bandwidth & !contains("_W")) %>%
  group_by( sample, varu, varv, covuv, fs_fun) %>%
  summarise_all(mean) %>%   # all the numbers are constant for a same design, so the mean does not change anything. It allows to regroup everything to avoid repeating design because of previous bandwidths values
  mutate(ols = ols_x/tsls_x,
         liml = LIML/tsls_x,
         WMDF = WMDF/tsls_x,
         np_x = np_x/tsls_x,
         ll_x = ll_x/tsls_x,
         lc_x = lc_x/tsls_x,
         tsls_x = tsls_x/tsls_x  ) 

relative_medAD <- MC_medAD %>% dplyr::select(!ends_with("var") & !bandwidth & !contains("_W")) %>%
  group_by( sample, varu, varv, covuv, fs_fun) %>%
  summarise_all(mean) %>%   # all the numbers are constant for a same design, so the mean does not change anything. It allows to regroup everything to avoid repeating design because of previous bandwidths values
  mutate(ols = ols_x/tsls_x,
         liml = LIML/tsls_x,
         WMDF = WMDF/tsls_x,
         np_x = np_x/tsls_x,
         ll_x = ll_x/tsls_x,
         lc_x = lc_x/tsls_x,
         tsls_x = tsls_x/tsls_x  ) 

# Save the relative stats tables in .csv format 
design <- "log"
stats <- c("relative_b", "relative_v", "relative_MSE", "relative_medbias", "relative_medAD")
set <- list(relative_b, relative_v, relative_MSE, relative_medbias, relative_medAD)  
for (i in 1:length(set))
{
  name <- paste(design, "_", stats[i], ".csv", sep = "")
  write_csv(set[[i]], name)
}

# Power graphs-------------------------------------

# unique_designs <- unique(listu[, c(2:6)])    # As the same designs show up multiple times due to the different values of the bandwidth, I report the unique designs
unique_designs <- unique(listu[, c(2:4)])    # In Julia

grid <- seq(-1, 3, 0.1) # grid centered around the true parameter beta = 1
# Saves a bunch of power graphs according to each design
1:nrow(unique_designs) %>% map( power_graph, data = listu, H0_range = grid, level = 0.05, by = "sample")
1:nrow(unique_designs) %>% map( power_graph, data = listu, H0_range = grid, level = 0.05, by = "bandwidth")

