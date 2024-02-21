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