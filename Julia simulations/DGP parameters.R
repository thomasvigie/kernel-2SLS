# This file lists the parameters used in the simulations --------------------------
nsim <- 10

param <- 1
#varz <- c(1,2)
varz <- 1
# N <- c(100, 500, 1000, 5000)
N <- c(100, 1000)
# K <- 10
bandrange <- c(.001, .005, .01, .05, .1, .5, 1, 10, 100 ) # for the bandwidths
bandrange <-  c(0.01, 0.1, 1, 10, 100)   # The Julia simulations only have 5 bandwidth values

# band <- c(0.01, 0.1, 0.5, 1, 1.5, 2, 3, 5, 10, 100) # for the leave-one-out selector
# type <- c("ll","lc") # we consider the local constant (Nadaraya-Watson estimator), and the local linear one.
# type <- "ll"
#cv <- c(TRUE, FALSE)
# fs_fun <- c("log", "linear", "A_L", "dnorm", "binary")
fs_fun <- "linear"
VCOV_list <- list()
VCOV_list[[1]] <- matrix(c( 1 , 0.1 ,
                            0.1 , 1 ) , nrow = 2, ncol = 2)
VCOV_list[[2]] <- matrix(c( 1 , 0.5 ,
                            0.5 , 1 ) , nrow = 2, ncol = 2)
VCOV_list[[3]] <- matrix(c( 1 , 0.9 ,
                            0.9 , 1 ) , nrow = 2, ncol = 2)
# VCOV_list[[4]] <- matrix(c( 1 , 0.5 ,
#                             0.5 , 0.5 ) , nrow = 2, ncol = 2)
# VCOV_list[[5]] <- matrix(c( 0.25 , 0.1 ,
#                             0.1 , 0.25 ) , nrow = 2, ncol = 2)
# VCOV_list[[6]] <- matrix(c( 0.5 , 0.5 ,
#                             0.5 , 1 ) , nrow = 2, ncol = 2)
# VCOV_list[[7]] <- matrix(c( 0.5 , 0.3 ,
#                             0.3 , 0.5 ) , nrow = 2, ncol = 2)

