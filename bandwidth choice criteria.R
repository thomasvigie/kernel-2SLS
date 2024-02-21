# Information criteria ---------------------------------

hurvic_AIC <- function(dependent, explanatory, bandwidth, type = c("lc", "ll"), loo = FALSE)
{
  # sets the default values for the type to local constant
  type <- match.arg(type)
  
  H <- kernel_hat_matrix(dependent, explanatory, bandwidth, type, loo)
  n <- nrow(as.matrix(dependent))  
  sig2 <- ( t(dependent)%*%t( diag(n) - H )  %*% ( diag(n) - H ) %*% dependent ) / n
  crit <- log(sig2) + (1 +  sum(diag(H))/ n)/( 1 - (sum(diag(H)) + 2) / n)
  return(crit)
}
# Example:
# u <- rnorm(100)
# x <- rnorm(100)
# y <- sin(x) + u
# hurvic_AIC(y, x, bandwidth = 0.1)
# hurvic_AIC(y, x, bandwidth = 0.1, loo = TRUE)   # The leave-one-out version



iv_AIC <- function(expl, endo, instru, bandwidth, type = c("lc", "ll"), loo = FALSE )
{
  # sets the default values for the type to local constant
  
  type <- match.arg(type)
  
  n <- nrow(as.matrix(expl))
  H <- kernel_hat_matrix(dependent = endo, explanatory = instru, bandwidth = bandwidth, type, loo)
  # In the IV case, the effective hat matrix is (X'H'X)^(-1) X'H'
  H <- endo%*%solve(t(endo)%*%t(H)%*%endo)%*%t(endo)%*%t(H)
  sig2 <- ( t(expl)%*%t( diag(n) - H )  %*% ( diag(n) - H ) %*% expl ) / n
  crit <- log(sig2) + (1 +  sum(diag(H))/ n)/( 1 - (sum(diag(H)) + 2) / n)
  
  return(crit)
}
# Example:
# mu <- c(0,0)
# VCOV = matrix(c( 1 , 0.5 ,
#                  0.5 , 1 ), nrow = 2, ncol = 2)
# errors <- MASS::mvrnorm(100, mu, Sigma = VCOV )
# z <- rnorm(100)
# x <- z + errors[, 2]
# y <- 2*x + errors[, 1]
# 
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1)
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1, loo = TRUE)
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1, type = "ll")
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1, type = "ll", loo = TRUE)

# Try minimizing the criteria wrt the bandwidth

# This functions minimize the criterion function

min_AIC_iv <- function(expl, endo, instru, type = c("lc", "ll"), loo = FALSE )
{
  type <- match.arg(type)
  
  # using optim
  # sol <- optim(par = h0, fn = iv_AIC, method = "Brent", expl = expl, endo = endo, instru = instru, type = type, loo = loo, lower = lower, upper = upper)
  # minimum <- sol$value
  # minimizer <- sol$par
  # convergence <- sol$convergence
  
  # using DEoptim
  # sol <- DEoptim(fn = iv_AIC, lower = lower, upper = upper, expl = expl, endo = endo, instru = instru, type = type, loo = loo, control = DEoptim.control(NP = 20, itermax = 50, trace = FALSE))
  # minimum <- sol$bestval
  # minimizer <- sol$`optim`$`bestmem`
  # convergence <- sol$iter
  
  # Using optimr
  sol <- optimr::optimr(par = c(1), fn = iv_AIC, method = "Nelder-Mead", expl = expl, endo = endo, instru = instru, type = "ll", loo = FALSE)
  
  return(sol)
}

# Example:
# mu <- c(0,0)
# VCOV = matrix(c( 1 , 0.5 ,
#                  0.5 , 1 ), nrow = 2, ncol = 2)
# errors <- MASS::mvrnorm(100, mu, Sigma = VCOV )
# z <- rnorm(100)
# x <- z + errors[, 2]
# y <- 2*x + errors[, 1]
# 
# min_AIC_iv(expl = y, endo = x, instru = z, lower = 0.001, upper = 10, type = "ll", loo = FALSE)

# With a differential evolution algorithm 
# library(DEoptim)
# DEoptim(fn = iv_AIC, lower = 0.0001, upper = 1, expl = y, endo = x, instru = z, type = "lc", loo = FALSE)

# k-fold cross validation -------------------------------------------
# this function takes the ith fold of a sample out, and computes the estimate with teh remaining sample
ith_fold_est <- function(y, x, z, h, type = c("lc", "ll"), i, k)
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  z <- as.matrix(z)
  n <- nrow(as.matrix(x))
  validation_size <- round(n/k)
  
  start <- (i - 1)*validation_size + 1
  end <- start + validation_size - 1
  x <- x[-c(start:end), ]
  y <- y[-c(start:end), ]
  z <- z[-c(start:end), ]
  res <- myest_np(y, x, z, h, type)
  return(res)
}

# this function loops over ith_fold to get the estimates with all the folds, and computes the resulting CV criterion using pre_est 
kfold_CV <- function( y, x, z, h, k, type)
{
  pre_est <- myest_np(y, x, z, h, type)
  # n <- nrow(as.matrix(x))
  # validation_size <- round(n/k)   
  # we map ith_fold so that it goes over each sample without the ith fold, and computes an estimate. So there will be k generated estimates.
  estimates <- 1:k %>% map_dbl(ith_fold_est, y = y, x = x, z = z, h = h, type = type, k = k) 
  
  crit <- mean( (estimates - c(pre_est))^2 )
  return(crit)
}

min_CV <- function( y, x, z, h, k, type)    # super slow though...
{
  sol <- optim(par = nrow(as.matrix(y))^(-1/5), kfold_CV, y = y, x = x, z = z, k = k, type = type)
  # minimum <- sol$par
  # argmin <- sol$value  
  return(list(minimum = sol$value, h_opti = sol$par, convergence = sol$convergence))
}
# Example
# n <- 100
# u <- rnorm(n)
# v <- 0.5*u + rnorm(n)
# z <- rnorm(n)
# x <- z + v
# y <- x + u
# 
# ith_fold_est(y, x, z, h = 1, type = "ll", i = 60, k = 100)
# fold_10 <- kfold_CV(y, x, z, h = 1, k = 10, type = "ll")
# kfold_CV(y, x, z, h = 1, k = n, type = "ll")
# min_CV(y = y, x = x, z = z, h = nrow(as.matrix(y))^(-1/5), k = 10, type = "ll")

