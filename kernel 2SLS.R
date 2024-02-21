# Two stage least squares estimator, where the first stage takes the fitted values of a Nadaraya-Watson or local linear estimator, and uses them in the second stage.
# The kernel regression in the first stage uses the np package of Hayfield and Racine

source("auxiliary functions.R")

kernel_2SLS <- function(Y, X, W = NULL, Z, h = "cv", type = c("lc", "ll"), intercept = c(TRUE, FALSE), kernel_type = c("gaussian", "epanechnikov"))
{
  type <- match.arg(type)
  kernel_type <- match.arg(kernel_type)
  library(np)
  
  if (h == "cv")
  {
    stage_1_bw <- npregbw( X ~ Z, regtype = type, ckertype = kernel_type)
    stage_1 <- npreg(bws = stage_1_bw)
    stage_1_bw <- stage_1_bw$bw
  }
  else {
    stage_1_bw <- h
    stage_1 <- npreg( X ~ Z, bws = stage_1_bw, regtype = type, ckertype = kernel_type) 
  }
  xhat <- as.matrix( stage_1$mean )
  
  X <- cbind(X, W)
  xhat <- cbind(xhat, W)
  if (intercept == TRUE)
  {
    int <- 1
    X <- cbind(int, X)
    xhat <- cbind(1, xhat)
    colnames(xhat)[1:2] <- c("intercept", "X")
  }
  colnames(xhat)[1] <- "X"
  bhat <- solve( t(xhat)%*%X )%*%( t(xhat)%*%Y ) 
  sig2hat <- mean((Y - X%*%c(bhat))^2)
  var_homo <- sig2hat*solve(( t(xhat)%*%xhat ))
  res <- Y - X%*%c(bhat)
  var_hetero <- solve(( t(xhat)%*%xhat ))%*% (t(xhat*rep.col(res^2, ncol(xhat)) )%*%xhat) %*%solve(( t(xhat)%*%xhat ))
  
  var_homo <- rbind(diag(var_homo))
  var_hetero <- rbind(diag(var_hetero))

  colnames(var_homo) <- rownames(bhat)
  colnames(var_hetero) <- rownames(bhat)
  
  results <- data.frame(Estimates = bhat,
                        'Std errors (homo)'  = t(sqrt(var_homo)),
                        'Std errors (hetero)'  = t(sqrt(var_hetero)),
                        check.names = FALSE   # check.names makes sure the spaces are not filled in the column names
                        )
  rownames(results) <- rownames(bhat)
  
  print(results)
  
  if (type == "lc")
  {type <- "Local constant"}
  else {type <- "Local linear"}
  
  cat("\n", fill = TRUE)  
  cat("Type of kernel regression in the first stage:")
  cat("\n")
  cat(type)
  cat("\n", fill = TRUE)  
  cat("Bandwidth:")
  cat("\n")
  cat(stage_1_bw)
  cat("\n", fill = TRUE)
  cat("Kernel function:")
  cat("\n")
  cat(kernel_type)
  
  return( list(beta_hat = bhat[, 1], varhat_homoskedasticity = var_homo[1,], varhat_heteroskedasticity = var_hetero[1,], Bandwidth = stage_1_bw) )
}

# Example:
n <- 100
VCOV <- matrix(c(1,0.5,0.5,1), 2, 2)
errors <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = VCOV  )
u <- errors[, 1]
v <- errors[, 2]
z <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x1 <- z + v
y <- 1 + x1 + x2 + x3 + u  # The true coefficients are all equal to 1

model_optimal_h <- kernel_2SLS(Y = y, X = x1, W = cbind(x2, x3), Z = z, type = "ll", intercept = TRUE, kernel_type = "epanechnikov")
print(model_optimal_h)

model <- kernel_2SLS(Y = y, X = x1, W = cbind(x2, x3), Z = z, h = 1, type = "ll", intercept = TRUE, kernel_type = "epanechnikov")
print(model)
