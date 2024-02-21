# this function implements the weighted minimum distance estimator proposed by Antoine and Lavergne (2014)
# Works with one instrument only for now

wmd <- function( x , y , z , exo = NULL, intercept = c("TRUE", "FALSE"))
{
  if(missing(intercept)){
    intercept <- TRUE
  }
  n <- nrow(as.matrix(x))
  
  if(intercept == TRUE){
    e <- matrix(1, nrow = n, ncol = 1 )
    Y1star <- cbind( e, x )
  }
  else {
    Y1star <- x
  }
  
  if(is.null(exo) == FALSE){
    Y1star <- cbind( Y1star, exo )
  }
  
  Ystar <- cbind( y , Y1star )
  
  if (is.null(exo) == FALSE)    # If there are exogenous variables 
  {
    Z <- cbind(z, exo)
    Z<- Z/rep.row((sqrt(diag(var(Z)))), n )    # Doing some re-scaling for the weights
    Ktilde <- matrix(1,nrow = nrow(Z),ncol = n)
    for(i in 1:ncol(Z))
    {
      K <- dnorm(rep.col(Z[,i],nrow(Z))-rep.row(t(Z[,i]) , n))
      K <- K - diag(diag(K))
      Ktilde <- Ktilde*K
    }
  } else {  
    z <- z/sd(z)
    Ktilde <- dnorm(rep.col(z,n)-rep.row(z,n))
    Ktilde <- Ktilde - diag(diag(Ktilde))        # the diagonal elements have to be 0
  }
  
  eig <- eigen(solve(t(Ystar)%*%Ystar)%*%(t(Ystar)%*%Ktilde%*%Ystar))
  lambdatilde <- min(eig$values)
  
  lambda_WMDF <- (lambdatilde - (1 - lambdatilde)/n)/(1 - (1 - lambdatilde)/n)
  
  WMD <- solve(t(Y1star)%*%(Ktilde-lambdatilde*diag(n))%*%Y1star)%*%(t(Y1star)%*%(Ktilde-lambdatilde*diag(n))%*%y)
  WMDF <- solve(t(Y1star)%*%(Ktilde-lambda_WMDF*diag(n))%*%Y1star)%*%(t(Y1star)%*%(Ktilde-lambda_WMDF*diag(n))%*%y)
  
  res <- list(WMD, WMDF)
  names(res) <- c("WMD", "WMDF")
  return(res)
}  

# Example
# n <- 100
# z1 <- rnorm(n)
# z2 <- rnorm( 100 , 2 )
# z <- z1
# u <- rnorm( n )
# v <- rnorm(n) + u
# w <- rnorm(n)
# x <- z1 + v
# y <-  2*x + w + u
# # 
# wmd( x = x , y = y , z = z1, exo = w, intercept = FALSE )