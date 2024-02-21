# Auxiliary functions for the simulations for kernel 2SLS

# Classical kernel functions for kernel regression
epanechnikov <- function (u)
{
  ker <- (3/4)*(1-u^2)*(as.numeric(I(abs(u)<=1)))
  return (ker)
}

gaussian <- function (u)
{
  ker <- (1/sqrt(2*pi))*exp( -u^2 /2)
  return (ker)
}

# MSE function for simulation purposes

MSE <- function(x , parameter, type = c("bias", "variance", "MSE", "medAD", "medbias"))
{
  if(missing(type)){type <- "MSE"}
  if(type == "bias")
  {y <- abs(mean(x) - parameter)  }
  else if(type == "variance")
  {y <- var(x)}
  else if(type == "medAD")
  {y <- median(abs(x - parameter))}
  else if(type == "medbias")
  {y <- median(x - parameter)}   # not sure about that, but makes sense
  else
  {y <- mean(( x - parameter )^2)}
  return(y)
}

# Function used for relative statistics (relative bias, variance, MSE) in the simulation study
relative_stat <- function(x, y) 
  {
    z <- x/y
    return(z)
  }

# gives names according to some character and the numeric variable h. Used for the simulation study
name_h <- function(h, char)    
{
  paste(char, h, sep = "_")
} 
 
# functions to repeat a vector in the repmat (Matlab) way. Used for the kernel estimators
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

rep.row<-function(x,n)    
{
  matrix(rep(x,each=n),nrow=n)
}

# Function that names the instruments Z for the simulation study
namingz <- function(z)
{
  numbers <- c(1:ncol(data.matrix(z)))
  names <- character(length = 0)
  for(i in 1:ncol(data.matrix(z)))
  {
    names[i] <- paste("z", as.character(numbers[i]), sep='')
  }
  return( names )
}









myest_np <- function (expl, endo, exo = NULL, instru, h, type = c("lc", "ll"), intercept = c(TRUE, FALSE)) # Same function as above, but uses the np package of Racine and Hayfield
{
  type <- match.arg(type)
  
  fitted <- np::npreg(bws = as.vector(h), txdat = instru, tydat = endo, ckertype = "gaussian", regtype = type)$mean
  # fitted <- lcll(endo, instru, h, type, cv)
  
  xhatkernel <- fitted
  
  X <- cbind(endo, exo)
  xhat <- cbind(xhatkernel, exo)
  colnames(xhat)[1] <- paste(type, "endo", sep = "_")
  if (intercept == TRUE)
  {
    int <- 1
    X <- cbind(int, X)
    xhat <- cbind(1, xhat)
    colnames(xhat)[1:2] <- paste(type, c("intercept", "endo"), sep = "_")
  }
  
  betatilde <- t(solve( t(xhat) %*% X) %*% t(xhat) %*% expl)
  colnames(betatilde) <- paste(type, colnames(X), sep = "_")
  sig2hat <- mean((expl - X%*%c(betatilde))^2)
  var <- sig2hat*solve(( t(xhat)%*%xhat ))
  var <- rbind(diag(var))
  colnames(var ) <- paste(type, colnames(X), "var", sep = "_")
  
  return( data.frame(betatilde, var ) ) 
}

# Example:
# n <- 100
# VCOV <- matrix(c(1,0.5,0.5,1), 2, 2)
# errors <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = VCOV  )
# u <- errors[, 1]
# v <- errors[, 2]
# z <- rnorm(n)
# x2 <- rnorm(n)
# x3 <- rnorm(n)
# x1 <- z + v
# y <- 1+ x1 + x2 + x3 + u
# # # # Check if they are the same. They should be for a Gaussian kernel, and the lc and ll options.
# model <- np_2SLS(y, endo = x1, exo = cbind(x2, x3), Z = z, type = "ll", intercept = TRUE)
# band <- model$opti_h
# myest_np(expl = y, endo = x, instru = z, h = band, type = "ll")
# myest(expl = y, endo = x, instru = z, h = band, type = "ll", cv = FALSE)

# myest_np(expl = y, endo = x, instru = z, h = 0.1, type = "ll")
# myest(expl = y, endo = x, instru = z, h = 0.1, type = "ll", cv = FALSE)


# Leave-one-out estimator ---------------------------------
# This function computes the hat matrix for a nonparametric kernel regression (local constant or local linear)
# There is also an option of having the leave-one-out hat matrix in each case.
kernel_hat_matrix <- function(dependent, explanatory, bandwidth, type = c("ll","lc"), loo = c(FALSE, TRUE) )
{
  type <- match.arg(type)
  
  n <- nrow(as.matrix(explanatory))
  
  firsterm <- apply(as.matrix(explanatory), 1, rep.col, n)
  secterm <- t(firsterm)
  U <- (firsterm-secterm)/bandwidth
  # K <- epanechnikov(U) # matrix [n,n]
  K <- gaussian(U)
  if (type == "ll")
  {
    Sn1 <- matrix( apply( K*(firsterm-secterm), 1, sum) , n, 1)
    Sn1 <- rep.col(Sn1 , n )
    Sn2 <- matrix( apply( K*(firsterm-secterm)^2 , 1, sum ) , n, 1) # should be column vectors. Both
    Sn2 <- rep.col(Sn2 , n)
    
    B <- K*(Sn2 - (firsterm-secterm)*Sn1)    # top term
    SB <- matrix(apply(B, 1, sum), nrow = n, ncol =  1)   # bottom term
    eta <- 1/n    # Ridging parameter
    Ksum <- matrix(apply(K, 1, sum) , nrow = n , ncol =  1)
    if(any(SB==0)==TRUE){
      ind <- which(SB==0)  # locate where the 0 are in the bottom term
      print(paste("Careful mate ! We applied some ridging because of a small bandwidth value at coordinates", ind))
      SB[ind]<- SB[ind] + eta*Ksum[ind]   # apply ridging to the bottom term
      B[ind, ] <- B[ind, ] + eta*K[ind, ] # apply ridging to the top term
    }  
    
    L <- B/rep.col(SB, n)
    
  }
  else 
  {
    SB <- matrix(apply(K, 1, sum) , nrow = n , ncol =  1)
    if(any(SB==0)==TRUE){
      print("Careful mate ! An adjustment was made to avoid Nan!")
      SB[which(SB==0)]<- 0.0001 
    } # if there is a sum that is equal to 0
    
    L <- K/rep.col(SB , n)
  }
  
  if (loo == TRUE)     # Implementing the leave one out version from the regular version (ridging implemented in the regular version above)
  {
    for (i in 1:n)
    {
      L[i,] <- L[i,]/sum(L[i,-i])
    }
    diag(L) <- 0     
  }
  
  return(L)
}

loocv <- function( dep, expl, bw, type = c("ll", "lc") )
{
n <- nrow(as.matrix(dep))
# L <- kernel_hat_matrix(dependent = dep, explanatory = expl, bandwidth = bw, type = type, loo = TRUE )
L <- kernel_hat_matrix(dependent = dep, explanatory = expl, bandwidth = bw, type = type, loo = FALSE )

# ind <- which(diag(L) == 0)
# eta <- 1/n
# if a diagonal element is equal to 1, it will screw up cross validation...
# L[ind] <- eta  # some more ridging


# LOO = zeros(n, n)
# for i in 1:n
#  LOO[i, :] = L[i, :]/sum(L[i, 1:end .!= i]) # Fix that issue here
#  LOO[i, i] = 0
# end
fit <- L%*%dep
# objective = sum((fit - dep).^2)
objective <- (1/n)*sum( ( (dep - fit)/ (1 - diag(L) ) ) ^2 )
# objective <- (1/n)*sum( (dep - fit)^2 )

return(objective)
}

# CV minimization
opti_h <- function(dep, expl, bw, type = c("ll", "lc") )
{
  # Different algorithms are proposed: golden search (optimize), DEoptim and Powell's method
  # hstar <- newuoa(par = c(0.001, 1, 10, 100, 1000, 1000000), fn = loocv, dep = dep, expl = expl, type = type)
  # hstar <- optim(par = 1, fn = loocv, method = "Brent", lower = 0.0001, upper = 10^3, dep = dep, expl = expl, type = type)
  # hstar <- DEoptim(fn = loocv, lower = 0.001, upper = 1000, dep = dep, expl = expl, type = type, control = DEoptim.control(NP = 20, itermax = 100, trace = FALSE))
  # hstar <- optimize(f = loocv, interval = c(0.0001, 10^6), expl = expl, dep = dep, type = type)
  hstar <- optimr::optimr(par = c(1), fn = loocv, method = "Nelder-Mead", expl = expl, dep = dep, type = type)
  # hstar <- 
  return(hstar)
}

# This function proposes a leave-one-out criterion for betahat 
# Formula taken from Raffaele Saggio (econometrica)
loo_MSE <- function(x, y, z, h, regtype = c("lc", "ll"))
{
  xhat <- lcll( dependent = x, explanatory = z, bandwidth = h, type = regtype, kertype = "gaussian", loo = FALSE )  # Ridging has been implemented
  # xhat <- np::npreg(bws = as.vector(h), txdat = z, tydat = x, ckertype = "gaussian", regtype = regtype)$mean
  HAT <- x%*%solve(t(xhat)%*%x)%*%t(xhat)
  betahat <-  myest(expl = y, endo = x, instru = z, h = h, type = regtype)   # add the variance etc like in myest_np
  betas_loo <- (x*betahat - y*diag(HAT) ) / ( (1 - diag(HAT))*x )
  MSE <- mean((betas_loo - c(betahat) )^2)
  
  # Equivalent of:
  # Sk <- t(xhat)%*%x
  # betas_loo[i] <- solve(Sk - xhat[i]*x[i])%*%t(xhat[-i])%*%y[-i]  
  return(MSE)
}

# This function minimizes the leave-one-out function loo_MSE
loo_choice <- function(x, y, z, h, regtype = c("lc", "ll"))
{
  loo <- h %>% map_dbl(loo_MSE, x = x, y = y, z = z, regtype = regtype)
  opt_h <- h[which.min(loo)]
  betahat <-  myest_np(expl = y, endo = x, instru = z, h = opt_h, type = regtype)
  return(data.frame(betahat = betahat, opti_h = opt_h))
}

loo_opti <- function(x, y, z, regtype = c("lc", "ll"))
{
  # res <- optimize(loo_MSE, interval = c(0.0001, 50), x = x, y = y, z = z, regtype = regtype)
  # res <- DEoptim(loo_MSE, lower = 0.001, upper = 10, DEoptim.control(trace = FALSE, itermax = 50), x = x, y = y, z = z, regtype = regtype )
  # opt_h <- res$`optim`$`bestmem`
  sol <- optimr::optimr(par = c(1), fn = loo_MSE, method = "Nelder-Mead", x = x, y = y, z = z, regtype = regtype)
  opt_h <- sol$par
  betahat <-  myest(expl = y, endo = x, instru = z, h = opt_h, type = regtype) 
  return(data.frame(betahat = betahat, opti_h = opt_h))
}


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


# MC graphs----------------------------------------

# Needed for straight lines in log scale
fa <- function(x, y, i, n) {
  xl <- seq(x[i - 1], x[i], (x[i] - x[i - 1]) / n)
  yl <- exp(log(y[i]) + (xl - x[i]) * (log(y[i]) - log(y[i - 1])) / (x[i] - x[i - 1]))
  return(data.frame(x = xl, y = yl))
}

# log_line <- function(x, y, n = 5) {
#   l <- lapply(2:length(x), fa, n)
#   return(do.call(rbind, l))
# }

# gni <- list()
# for (i in 1:3)
# {
# gni[[i]] <- map(2:9, fa, x = da$x[(9*(i - 1) + 1):(9*(i - 1) + 9)], y = da$ll[(9*(i - 1) + 1):(9*(i - 1) + 9)], n = 5)
# gni[[i]] <- Reduce(rbind, gni[[i]])
# }
# gni <- Reduce(rbind, gni)


# MC_MSE <- filter(MC_MSE, bandwidth != 5)    # for the binary_x case

MC_graph <- function(MC, j, exo = c(FALSE, TRUE))  # graph of different statistics according to the value of the bandwidths
{

  # MC should be a list of all the different stats we want to graph: bias, variance, MSE, MAD, median bias,...
  b <- length(bandrange) 
  statistics <- c("Bias", "Variance", "MSE")
  # If looking at different bias measures
  # statistics <- c("Bias", "Median bias", "Median absolute deviation")
  da <- list()
  
  for (i in 1:3){
    # for ( i in 3){      # when making the zooms of the MSE ones
    data_to_graph <- MC[[i]] [ , -(1:5)]
    data_to_graph <- dplyr::filter(MC[[i]], sample == unique_designs[j,]$sample, 
                                   varu == unique_designs[j ,]$varu, 
                                   fs_fun == unique_designs[j ,]$fs_fun,
                                   varv == unique_designs[j,]$varv,
                                   covuv == unique_designs[j,]$covuv)
    data_to_graph <- data_to_graph %>% mutate(stat = rep(statistics[i], b)) %>% mutate(x = c(1:b))
    
    da[[i]] <- data_to_graph
  }
  da <- reduce(da, rbind)
  da$stat <- factor(da$stat, levels = c("Bias", "Variance", "MSE"))
  # If looking at different bias measures
  # da$stat <- factor(da$stat, levels = c("Bias", "Median bias", "Median absolute deviation"))
  if(exo == FALSE)
  {ex <- NULL}
  else{ex <- ", with exogenous variables"}
  
  name_graph <- paste(unique_designs[j, ]$fs_fun,
                      ", n=", unique_designs[j, ]$sample,
                      ", varu=", unique_designs[j, ]$varu,
                      ", varv=", unique_designs[j, ]$varv,
                      ", covuv=", unique_designs[j, ]$covuv, 
                      ex,
                      ".jpg", sep = "" )
  
  
  # Some shenanigans to get straight lines in log scale plots
  npoints <- 20
  gnilc <- list()
  for (i in 1:3)
  {
    gnilc[[i]] <- map(2:9, fa, x = da$x[(9*(i - 1) + 1):(9*(i - 1) + 9)], y = da$lc[(9*(i - 1) + 1):(9*(i - 1) + 9)], n = npoints)
    gnilc[[i]] <- Reduce(rbind, gnilc[[i]])
  }
  gnilc <- Reduce(rbind, gnilc)  
  names(gnilc) <- c("x", "ylc")
  
  gnill <- list()
  for (i in 1:3)
  {
    gnill[[i]] <- map(2:9, fa, x = da$x[(9*(i - 1) + 1):(9*(i - 1) + 9)], y = da$ll[(9*(i - 1) + 1):(9*(i - 1) + 9)], n = npoints)
    gnill[[i]] <- Reduce(rbind, gnill[[i]])
  }
  gnill <- Reduce(rbind, gnill) 
  names(gnill) <- c("x", "yll")
  
  gni <- cbind(gnilc, gnill)
  gni <- gni[!duplicated(as.list(gni))]
  stat <- rep(c("Bias", "Variance", "MSE"), each = (length(bandrange) - 1)*(npoints + 1))
  datu <- cbind(gni, stat)
  
  # Ggplot time
  
  ggplot(data = da, aes(x = x)) +
    geom_line (aes( y = WMDF, color = "WMDF"), size = 1 ) +
    geom_line (aes( y = ols, color =  "OLS"), size = 1 ) +
    geom_line (aes( y = tsls, color =  "2SLS" ), size = 1 ) +
    geom_line (aes( y = np_tsls, color =  "NP2SLS" ), size = 1) +
    # geom_line (aes( y = np_endo, color =  "NP2SLS" ), size = 1) +
    
    # geom_line (aes( y = LIML, color =  "LIML" ), size = 1) +
    
    geom_line (aes( y = CSA, color =  "CSA" ), size = 1) +
    geom_line (aes( y = DN, color =  "DN" ), size = 1) +
    
    # omit for the linear design case
    # geom_point (aes( y = lc, color = "mybeta_lc") , size = 2, shape = 15)  +
    # geom_line(data = datu, aes( x = x, y = ylc, color = "mybeta_lc"), size = 0.8)  +
    
    # geom_point(aes( y = ll, color = "mybeta_ll") , size = 2)  +
    # geom_line(aes( y = ll, color = "mybeta_ll") , size = 0.8)  +
    geom_point (aes( y = ll, color = "mybeta_ll") , size = 2, shape = 16)  +
    geom_line(data = datu, aes( x = x, y = yll, color = "mybeta_ll"), size = 0.8)  +
    
    
    # geom_point (aes( y = ll_endo, color = "mybeta_ll") , size = 2)  +
    # geom_line(aes( y = ll_endo, color = "mybeta_ll") , size = 0.8)  +
    
    
    xlab("Bandwidth") +
    ylab("") +
    theme( axis.line = element_line(colour = "darkblue", 
                                    size = 1, linetype = "solid")) +
    # scale_x_discrete(breaks = as.character(1:b),
    #                  labels = as.character(bandrange) ) +
    scale_x_discrete(limits = as.character(1:b), breaks = as.character(1:b) ,
                     labels =  c(".001", ".005", ".01", ".05", ".1",
                                 ".5", "1", 
                                 #"2", 
                                 #"5",
                                 "10",
                                 #"20",
                                 "100")
    ) +
    scale_color_manual( breaks = c( "2SLS", 
                                    "CSA",
                                    "DN",
                                    # "IV_K", 
                                    # "LIML",
                                    # "mybeta_lc",     # omit in the linear design case
                                    "mybeta_ll",
                                    "NP2SLS",
                                    "OLS",
                                    # "P2SLS", 
                                    "WMDF"),       # Put them in alphabetical order,it is easier for colours.
                        labels = list(TeX('$\\hat{\\beta}_{2SLS}$'),
                                      TeX('$\\hat{\\beta}_{CSA}$'),
                                      TeX('$\\hat{\\beta}_{DN}$'),
                                      # TeX('$\\hat{\\beta}_{LIML}$'),
                                      # TeX('$\\hat{\\beta}_{lc}$'),       # omit in the linear design case
                                      TeX('$\\hat{\\beta}_{ll}$'),
                                      TeX('$\\hat{\\beta}_{NP}$'),
                                      TeX('$\\hat{\\beta}_{OLS}$'), 
                                      # "P2SLS", 
                                      TeX('$\\hat{\\beta}_{WMDF}$') ),
                        values = c("#E69F00",
                                   "saddlebrown",   # CSA
                                   "deeppink4",   # DN
                                   # "darkgreen",    # LIML
                                   # "blue",                 # for beta_lc    # omit in the linear design case
                                   "orangered3",            # for beta_ll
                                   #"purple", 
                                   # "#00BFC4",   # other color for betahat_NP
                                   "yellow2",
                                   "olivedrab3", # "#00BA38", # "lawngreen",  
                                   # "#00BFC4", 
                                   "#00BFC4" # "#3300CC" # "#C77CFF"    # for WMDF
                        )) +      
    
    
    theme(axis.text.x = element_text(angle = 45)) +   # to make x-axis values bent
    
    theme(legend.title=element_blank(),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 6)) +
    
    facet_wrap(~stat , nrow = 1, scales = "free")   +
    coord_trans(y = "log10") +
    guides(colour = guide_legend(override.aes = list(shape = c(NA,     # modifying legend by hand
                                                               NA,     # CSA
                                                               NA,     # DN
                                                               # NA,     # LIML
                                                               # 15,   # lc     # omit in the linear design case
                                                               16,     # ll
                                                               NA,
                                                               NA,
                                                               NA)   ))) +
    labs(colour = "Estimators") +
    
    theme(axis.title.x = element_text(family = "serif"),   # Changes fonts into times new roman
          axis.title.y = element_text(family = "serif"),
          legend.text = element_text(family = "serif"),
          legend.title = element_text(family = "serif"),
          strip.text.x = element_text(family = "serif")) +
    
    
    theme(aspect.ratio=1)     # for the non-zoomed ones
  
  ggsave(filename = name_graph, device = "jpg", height = 6, units = "cm" )
  
  # on.exit(dev.off())
  graphics.off()
  return("done !")
  
}

# Hypothesis tests & power -----------

t_test <- function(betahat, se, para, level = 0.05)
{
  if(level <0 || level >1)
    error("the level should be between 0 and 1 !!")
  tvalue <- (betahat - para)/se 
  decision <- abs(tvalue) > qnorm(1 - level/2)
  decision <- ifelse(decision == TRUE, "reject", "not reject") 
  return(decision)
}
# Example
# t_test(betahat = est, se = se, para = 1, level = 0.05)


H0_test <- function(parameter, estimate, std_error, level = 0.05, df)
{
  name <- paste("H0(", parameter, ")", sep = "")
  df <- df %>% mutate(name = t_test(estimate, se = std_error, para = parameter, level = level))
  df[, ncol(df)] <- ifelse(df[, ncol(df)] == "reject", 1, 0) 
  names(df)[ncol(df)] <- name
  return(df[names(df)==name])
}

# listu <- grid %>% map(H0_test, estimate = listu$TV, std_error = listu$TV_se, level = 0.05, df = listu)

test_to_dat <- function(parameter, estimate, std_error, level = 0.05, df, nsim)
{
  dat <- parameter %>% map(H0_test, estimate = estimate, std_error = std_error, level = level, df = df)
  dat <- purrr::reduce(dat, cbind)
  
    dat <- cbind(df [, which(grepl("sample|varu|varv|covuv|fs_fun|bandwidth", names(df)) == TRUE)], dat)
    tests <- dat [, which(grepl("sample|varu|varv|covuv|fs_fun|bandwidth|H0", names(dat)) == TRUE)] %>%
             group_by( sample, varu, varv, covuv, fs_fun, bandwidth ) %>% 
             summarise_all( freq, S = nsim )  
 # nsim <- nrow(df) / nrow(tests)
  #tests[, which(grepl("H0", names(tests)) == TRUE)] <- tests[, which(grepl("H0", names(tests)) == TRUE)] / nsim
  return(tests)
}

# Power graphs (general estimators) ------------------------
# This function computes the rejection frequencies for each design, each bandwidth value and some Hypothesis H0
rej <- function(data, H0, est, se, level = 0.05, name = NULL)
{
  nsim <- max(data$simu)
 # Use left_join or right_join when adding the new power columns
  tests <- data %>% mutate( power = ifelse(t_test(betahat = est, se = se, para = H0, level = level) == "reject", 1, 0)) 
  
  rejections <- tests %>%dplyr::select(sample, varu, varv, covuv, fs_fun, bandwidth, ends_with("power") )%>%
    group_by(sample, varu, varv, covuv, fs_fun, bandwidth)%>%
    summarize(power = sum(power)/nsim)
  
if(is.null(name) == FALSE)
{
  names(rejections)[names(rejections) == 'power'] <- name
}
  
  rejections <- rejections %>% mutate(H0 = H0)
  return(rejections)
}

# Example
# powa_ols <- rej(data = listu, H0 = 1, est = listu$ols, se = listu$ols_se, level = 0.05)
# powa_tsls <- rej(data = listu, H0 = 1, est = listu$tsls, se = listu$tsls_se, level = 0.05)


power_graph <- function(j, data, H0_range, level = 0.05, by = c("sample", "bandwidth"))
{
  
 # b <- length(bandrange) 
  # statistics <- c("Bias", "Variance", "MSE")
#  da <- list()
  
 # for (i in 1:b){
    # for ( i in 3){      # when making the zooms of the MSE ones
  #   data_to_graph <- MC[[i]] [ , -(1:5)]
  #   data_to_graph <- dplyr::filter(MC[[i]], sample == unique_designs[j,]$sample, 
  #                                  varu == unique_designs[j ,]$varu, 
  #                                  fs_fun == unique_designs[j ,]$fs_fun,
  #                                  varv == unique_designs[j,]$varv,
  #                                  covuv == unique_designs[j,]$covuv)
  #   data_to_graph <- data_to_graph %>% mutate(stat = rep(statistics[i], b)) %>% mutate(x = c(1:b))
  #   
  #   da[[i]] <- data_to_graph
  # }
  # da <- reduce(da, rbind)
  # da$stat <- factor(da$stat, levels = c("Bias", "Variance", "MSE"))
  
# That selects a design and a bandwidth
unique_designs <- unique(data[, c(2:6)])    # As the same designs show up multiple times due to the different values of the bandwidth, I report the unique designs
dat <- list()
if( by == "sample")
{
# Power curves for each sample, one bandwidth
for (i in 1:length(N))
{
  datum <- dplyr::filter( data, sample == N[i],
                          varu == unique_designs[j ,]$varu,
                          fs_fun == unique_designs[j ,]$fs_fun,
                          varv == unique_designs[j,]$varv,
                          covuv == unique_designs[j,]$covuv,
                          bandwidth == bandrange[1] )  # We loop over sample sizes, not bandwidth in this chunk of code
  
  # dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols, se = datum$ols_se, level = level, name = "ols")
  # dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls, se = datum$tsls_se, level = level, name = "tsls")
  # dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_endo, se = sqrt(datum$np_endo_var), level = level, name = "np")
  # dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_endo, se = sqrt(datum$ll_endo_var), level = level, name = "ll")
  # dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_endo, se = sqrt(datum$lc_endo_var), level = level, name = "lc")
 
  # If we look at the exogenous variables simulation
  dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols_x, se = sqrt(datum$ols_x_var), level = level, name = "ols")
  dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls_x, se = sqrt(datum$tsls_x_var), level = level, name = "tsls")
  dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_x, se = sqrt(datum$np_x_var), level = level, name = "np")
  dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_x, se = sqrt(datum$ll_x_var), level = level, name = "ll")
  dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_x, se = sqrt(datum$lc_x_var), level = level, name = "lc")
  
   
  dat[[i]] <- dat_ols %>% left_join(dat_tsls, by = NULL) %>%
    left_join(dat_np, by = NULL) %>%
    left_join(dat_ll, by = NULL) %>%
    left_join(dat_lc, by = NULL) 
}
  name_graph <- paste("power graph, ", unique_designs[j, ]$fs_fun,
                      # "n=", unique_designs[j, ]$sample,
                      ", varu=", unique_designs[j, ]$varu,
                      ", varv=", unique_designs[j, ]$varv,
                      ", covuv=", unique_designs[j, ]$covuv, 
                      ", bandwidth=", bandrange[1],
                      ", samples",
                      ".jpg", sep = "")
}
else if (by == "bandwidth")
{
# Power curves for each bandwidth, one sample size
for (i in 1:length(bandrange))
{
    datum <- dplyr::filter( data, sample == unique_designs[j,]$sample,
                                     varu == unique_designs[j ,]$varu,
                                     fs_fun == unique_designs[j ,]$fs_fun,
                                     varv == unique_designs[j,]$varv,
                                     covuv == unique_designs[j,]$covuv,
                                     bandwidth == bandrange[i] )  # The bandwidth is fixed a well as the power graphs considered don't depend on the bandwidth

  # dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols, se = datum$ols_se, level = level, name = "ols")
  # dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls, se = datum$tsls_se, level = level, name = "tsls")
  # dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_endo, se = sqrt(datum$np_endo_var), level = level, name = "np")
  # dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_endo, se = sqrt(datum$ll_endo_var), level = level, name = "ll")
  # dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_endo, se = sqrt(datum$lc_endo_var), level = level, name = "lc")
# If we look at the exogenous variables simulation
  dat_ols <- H0_range %>% map_dfr(rej, data = datum, est = datum$ols_x, se = sqrt(datum$ols_x_var), level = level, name = "ols")
  dat_tsls <- H0_range %>% map_dfr(rej, data = datum, est = datum$tsls_x, se = sqrt(datum$tsls_x_var), level = level, name = "tsls")
  dat_np <- H0_range %>% map_dfr(rej, data = datum, est = datum$np_x, se = sqrt(datum$np_x_var), level = level, name = "np")
  dat_ll <- H0_range %>% map_dfr(rej, data = datum, est = datum$ll_x, se = sqrt(datum$ll_x_var), level = level, name = "ll")
  dat_lc <- H0_range %>% map_dfr(rej, data = datum, est = datum$lc_x, se = sqrt(datum$lc_x_var), level = level, name = "lc")
  
  
  
dat[[i]] <- dat_ols %>% left_join(dat_tsls, by = NULL) %>%
                  left_join(dat_np, by = NULL) %>%
                  left_join(dat_ll, by = NULL) %>%
                  left_join(dat_lc, by = NULL)
}
  name_graph <- paste("power graph, ", unique_designs[j, ]$fs_fun,
                      ", n=", unique_designs[j, ]$sample,
                      ", varu=", unique_designs[j, ]$varu,
                      ", varv=", unique_designs[j, ]$varv,
                      ", covuv=", unique_designs[j, ]$covuv,
                      ", bandwidths",
                      ".jpg", sep = "")
}
dat <- reduce(dat, rbind)

  #unique_designs <- unique(dat[, 1:5])    # As the same designs show up multiple times due to the different values of the bandwidth, I report the unique designs
  
#  data_to_graph <- dplyr::filter(dat, sample == unique_designs[j,]$sample, varu == unique_designs[j ,]$varu, fs_fun == unique_designs[j ,]$fs_fun,
                                 # varv == unique_designs[j,]$varv,
                                 # covuv == unique_designs[j,]$covuv)

  
p <-  ggplot(data = dat, aes(x = H0)) +
    geom_line (aes( y = ols, color = "ols"), size = 1 ) +
    geom_line (aes( y = tsls, color = "tsls"), size = 1 ) +
    geom_line (aes( y = np, color = "np_power"), size = 1 ) +
    # geom_line (aes( y = ll, color = "ll"), size = 1 ) +
    # geom_line (aes( y = lc, color = "lc"), size = 1 ) +
    
    
    # xlab("H0") +
    xlab(TeX('$H_0$')) +
    ylab("Probability of rejection") +
    theme( axis.line = element_line(colour = "darkblue", 
                                    size = 1, linetype = "solid")) +
    # scale_x_discrete(breaks = as.character(1:b),
    #                  labels = as.character(bandrange) ) +
    # scale_x_discrete(limits = as.character(1:b), breaks = as.character(1:b) ,
    #                  labels =  c(".001", ".005", ".01", ".05", ".1",
    #                              ".5", "1", 
    #                              #"2", 
    #                              #"5",
    #                              "10",
    #                              #"20",
    #                              "100")
    # ) +
  scale_color_manual( breaks = c(  #"lc",
                                   # "ll",
                                   "np_power",
                                   "ols",
                                   "tsls" ),       # Put them in alphabetical order,it is easier for colours.
                      labels = list( #TeX('$\\hat{\\beta}_{lc}$') ,
                                     #TeX('$\\hat{\\beta}_{ll}$') ,
                                     TeX('$\\hat{\\beta}_{NP}$'),
                                     TeX('$\\hat{\\beta}_{OLS}$') ,
                                     TeX('$\\hat{\\beta}_{2SLS}$')
                                     ),
                      
                      values = c(
                        # "#663300", 
                        #"#3300CC", 
                        # "yellow",
                        # "#00BFC4",   # Some kind of green-blue     # for beta_lc
                        # "blue",          # for beta_lc
                        # "orangered3",       # for beta_ll         
                        # "orange",       # For for betahat_NP      
                        "yellow2",       # For for betahat_NP      
                        
                       # "purple", 
                        # "#00BFC4",   
                        # "cyan",       # for OLS
                         # "yellow2",    # 2SLS
                        "olivedrab3",   # OLS
                       "#E69F00"    # 2SLS
                       # "#00BA38", # "lawngreen",  
                        # # "#00BFC4", 
                        # "#00BFC4" # "#3300CC" 
                       # "#C77CFF"    # for WMDF
                      )) +      
     labs(color = "Estimators") +
     scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1)) +
     theme(axis.text.x = element_text(angle = 45))

if(by == "sample")  
  
{    p <- p + facet_wrap(~ sample, scales = "free", labeller = as_labeller(c(`100` = "n = 100", `1000` = "n = 1000") ) ) } 
  else if (by == "bandwidth")    
{    
    p <- p + facet_wrap(~ bandwidth, scales = "free", labeller = as_labeller(c(`0.001` = "h = 0.001", `0.005` = "h = 0.005", `0.01` = "h = 0.01", `0.05` = "h = 0.05", `0.1` = "h = 0.1", `0.5` = "h = 0.5", `1` = "h = 1", `10` = "h = 10", `100` = "h = 100") ) )    
}     
p <- p + theme(axis.title.x = element_text(family = "serif"),   # Changes fonts into times new roman
            axis.title.y = element_text(family = "serif"),
            legend.text = element_text(family = "serif"),
            legend.title = element_text(family = "serif"),
            strip.text.x = element_text(family = "serif")) +    
      theme(aspect.ratio=1)     # for the non-zoomed ones
if(by == "sample")  
  
{    ggsave(filename = name_graph, device = "jpg", height = 6, units = "cm" ) } 
else if (by == "bandwidth")    
{    
  ggsave(filename = name_graph, device = "jpg")
}      
  # ggsave(filename = name_graph, device = "jpg", height = 6, units = "cm" )
  # ggsave(filename = name_graph, device = "jpg")
  
  # on.exit(dev.off())
  graphics.off()
  return("done !")
    
}




# Example
# n <- 100
# VCOV = matrix(c( 1 , 0.5 ,
#                  0.5 , 1 ), 2, 2 )
# mu <- c( 0 , 0 )  # mean of u and v for the joint distribution
# RNG <- MASS::mvrnorm( n , mu , Sigma = VCOV )   # or any other distribution. A zero covariance implies independence in this case though, a very nice feature
# u <- RNG [ , 1 ]
# v <- RNG [ , 2 ]
# 
# z <- rnorm(n)
# x <- z + v
# y <- x + u
# # seems to work
# band <- c(0.01, 0.1, 0.5, 1, 1.5, 2, 3, 5, 10, 100)
# band %>% map_dbl(loo_MSE, x = x, y = y, z = z, regtype = "ll")
# loo_choice(x = x, y = y, z = z, h = band, regtype = "ll")

# Small exercise about leave one out estimators using the formula provided in Wasserman and Raffaele Saggio----------------
 
# n <- 100
# VCOV = matrix(c( 1 , 0.5 ,
#                  0.5 , 1 ), 2, 2 )
# mu <- c( 0 , 0 )  # mean of u and v for the joint distribution
# RNG <- MASS::mvrnorm( n , mu , Sigma = VCOV )   # or any other distribution. A zero covariance implies independence in this case though, a very nice feature
# u <- RNG [ , 1 ]
# v <- RNG [ , 2 ]
# 
# z <- rnorm(n)
# x <- z + v
# y <- x + u
# # 
# # 
# h <- 100
# # # xhat <- lcll( dependent = z, explanatory = x, bandwidth = h, type = "ll", loo = FALSE )
# #   
# xhat <- np::npreg(bws = as.vector(h), txdat = z, tydat = x, ckertype = "gaussian", regtype = "ll")$mean
# H <- kernel_hat_matrix(dependent = x, explanatory = z, bandwidth = h, type = "ll", loo = FALSE )
# #   
# beta_ols <- lm(y ~ -1 +x)$coeffcients
# # betahat <-  myest(expl = y, endo = x, instru = z, h = h, type = "ll")
# # betahat <- solve(t(x)%*%t(H)%*%x)%*%t(x)%*%t(H)%*%y
# betahat <- solve(t(xhat)%*%x)%*%t(xhat)%*%y
# PZ <- z%*%solve(t(z)%*%z)%*%t(z)
# beta_iv <- solve(t(x)%*%PZ%*%x)%*%t(x)%*%PZ%*%y
# # HAT <- x%*%solve(t(x)%*%t(H)%*%x)%*%t(x)%*%t(H)
# HAT <- x%*%solve(t(xhat)%*%x)%*%t(xhat)
# # 
# LHS <- matrix(0, n, 1)
# betal <- matrix(0, n, 1)
# beta_iv_l <- matrix(0, n, 1)
# beta_ols_l <- matrix(0, n, 1)
# for (i in 1:n)
# {
# beta_ols_l <- lm(y[-i] ~ -1 + x[-i])$coefficients
# Pz <- z[-i]%*%solve(t(z[-i])%*%z[-i])%*%t(z[-i])
# beta_iv_l[i] <- solve(t(x[-i])%*%Pz%*%x[-i])%*%t(x[-i])%*%Pz%*%y[-i]
# # betal[i] <- myest(expl = y[-i], endo = x[-i], instru = z[-i], h = h, type = "ll")
# # formula from Raffaele Saggio's paper
# Sk <- t(xhat)%*%x
# betal[i] <- solve(Sk - xhat[i]*x[i])%*%t(xhat[-i])%*%y[-i]    # that is the right formula !!
# LHS[i] <- y[i] - x[i]*betal[i]
# }
# # RHS <- (y - x*betahat)/(1 - diag(HAT) )
# # 
# # RHS - LHS
# PX <- x%*%solve(t(x)%*%x)%*%t(x)
# beta_ols_loo <- (x*beta_ols - y*diag(PX) ) / ( (1 - diag(PX))*x )
# beta_ols_loo - beta_ols_l # the difference is 0
# # the difference is pretty much 0
# beta_iv_loo <- (x*beta_iv - y*diag(x%*%solve(t(x)%*%PZ%*%x)%*%t(x)%*%PZ) ) / ( (1 - diag(x%*%solve(t(x)%*%PZ%*%x)%*%t(x)%*%PZ))*x )
# beta_iv_loo - beta_iv_l
# # the difference is too big...
# betaloo <- (x*betahat - y*diag(HAT) ) / ( (1 - diag(HAT))*x )
# betaloo - betal
