library("np")
library(tidyverse)
library(DEoptim)
data <- read_csv("try data.csv")
source("auxiliary functions.R")

# OLS
summary(lm(y ~ x, data = data))
# 2SLS
summary(AER::ivreg(y ~ x|z, data = data))
# WMDF
wmd( x = data$x , y = data$y , z = data$z, exo = NULL, intercept = FALSE )

pred <- lcll( dependent = data$x, explanatory = data$z, bandwidth = 0.001, type = "ll", kertype = "gaussian", loo = FALSE )  # Ridging has been implemented
np_pred <- npreg(bws = 0.001, txdat = data$z, tydat = data$x, regtype = "ll")$mean
np_pred - pred


# With the Epanechnikov kernel. Does not seem to be the same thing
pred <- lcll( dependent = data$x, explanatory = data$z, bandwidth = 0.1, type = "ll", kertype = "epanechnikov", loo = FALSE )  # Ridging has been implemented
np_pred <- npreg(bws = 0.1, txdat = data$z, tydat = data$x, regtype = "ll", ckertype = "epanechnikov")$mean
np_pred - pred

# Check leave-one-out CV criterion function
loocv(dep = data$x, expl = data$z, bw = 0.001, type = "ll")
loocv(dep = data$x, expl = data$z, bw = 1, type = "ll")


# Bandwidth computation
np_bw <- npregbw( xdat = data$z, ydat = data$x, regtype = "ll")
my_bw <- opti_h(dep = data$x, expl = data$z, type = "ll" )
np_bw
my_bw

# Things match pretty well between np and my functon opti_h, although the optimal bandwidth is different
loocv( dep = data$x, expl = data$z, bw = my_bw$par, type = "ll" ) # loo CV criterion at the optimal "opti_h" bandwidth
loocv( dep = data$x, expl = data$z, bw = np_bw$bw, type = "ll" ) # loo CV criterion at the optimal np bandwidth
lcll( dependent = data$x, explanatory = data$z, bandwidth = my_bw$par, type = "ll", kertype = "gaussian", loo = FALSE )  # Ridging has been implemented
lcll( dependent = data$x, explanatory = data$z, bandwidth = np_bw$bw, type = "ll", kertype = "gaussian", loo = FALSE )  # Ridging has been implemented



# On a bigger data set
n <- 1000
data <- dgp(n = n, g = "A&L")
# Bandwidth computation
np_bw <- npregbw( xdat = data$z1, ydat = data$x, regtype = "ll")
my_bw <- opti_h(dep = data$x, expl = data$z1, type = "ll" )
np_bw
my_bw  

loocv( dep = data$x, expl = data$z, bw = my_bw$par, type = "ll" ) # loo CV criterion at the optimal "opti_h" bandwidth
loocv( dep = data$x, expl = data$z, bw = np_bw$bw, type = "ll" ) # loo CV criterion at the optimal np bandwidth
head(lcll( dependent = data$x, explanatory = data$z, bandwidth = my_bw$par, type = "ll", kertype = "gaussian", loo = FALSE ))  # Ridging has been implemented
head(lcll( dependent = data$x, explanatory = data$z, bandwidth = np_bw$bw, type = "ll", kertype = "gaussian", loo = FALSE ) ) # Ridging has been implemented

# Bandwidth computation with AIC criterion (Hurvic, Simonov and Tsai)
np_bw <- npregbw( xdat = data$z, ydat = data$x, regtype = "ll", bwmethod = "cv.aic")
my_bw <- optimr::optimr(par = c(1), fn = hurvic_AIC, method = "Nelder-Mead", explanatory = data$z, dependent = data$x, type = "ll", loo = FALSE)
# They match!! 
np_bw
my_bw

# Try the iv version now. Works bof so far (find negative bandwidth sometimes)
min_AIC_iv(expl = data$y, endo = data$x, instru = data$z1, type = "ll", loo = FALSE )
  
# Try the loo CV for betahat. It works!
loo_opti(x = data$x, y = data$y, z = data$z1, regtype = "ll")
np_2SLS(Y = data$y, endo = data$x, exo = NULL, Z = data$z1, type = "ll", intercept = FALSE)
  
  
