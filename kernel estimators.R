# Kernel estimators: local constant and local linear
# Ridging is applied when numerical issues are likely. The ridging follows the npreg functions from the np package of Hayfield and Racine: https://cran.r-project.org/web/packages/np/vignettes/np.pdf
# The function returns the predictions yhat = fhat(x)

source("auxiliary functions.R")  # Need some of the functions from the auxiliary functions script

lcll <- function ( dependent, explanatory, bandwidth, type = c("ll","lc"), kertype = c("gaussian", "epanechnikov"))
{
  type <- match.arg(type) 
  kertype <- match.arg(kertype)
  
  n <- nrow(as.matrix(explanatory))
  
  firsterm <- rep.col(as.matrix(explanatory), n)
  secterm <- t(firsterm)
  U <- (firsterm-secterm)/bandwidth
  
  if (kertype == "epanechnikov")
  {
    K <- epanechnikov(U) # matrix [n,n]
  }
  else{
    K <- gaussian(U)
  }
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
  
  fit <- L%*%as.matrix(dependent)
  return(fit)
}


# Example:
# n <- 100
# x <- rnorm(n)
# u <- rnorm(n)
# y <- sin(x) + u 
# lcll(dependent = y, explanatory = x, bandwidth = 1, type = "ll", kertype = "epanechnikov")
