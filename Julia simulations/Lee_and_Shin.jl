# This script is the translation of the R code provided in Lee & Shin (2019)

# Lee and Shin (2018) codes -------------------------------

function gen_pi_con(K, rho_z, R_sq_f)
  pi_con = sqrt( R_sq_f / ( ( K + K * ( K - 1 ) * rho_z ) * ( 1 - R_sq_f ) ) )*ones(1, K)
  return pi_con
end

function gen_pi_dec(K, rho_z, R_sq_f)
  sum1 = sum2 = 0
  for ind_k = 1:K
    sum1 = sum1 + ( 1 - ind_k / (K+1) )^8
    Set_for_j = collect(1: K)
    Set_for_j = setdiff(Set_for_j, ind_k)    #[-ind_k]
    for ind_j in Set_for_j
    sum2 = sum2 + ( 1 - ind_k / ( K + 1 ) )^4 * ( 1 - ind_j / (K+1) )^4 * rho_z
    end
  end

  c_K_dec = sqrt( ( R_sq_f / ( 1-R_sq_f ) ) * ( 1 / (sum1 + sum2) ) )

  pi_dec = c_K_dec .* ( 1 .- collect(1:K) ./ ( K + 1 ) ).^4

  return pi_dec
end


function gen_pi_irr(K, rho_z, R_sq_f)
  half_K = K / 2
  sum1 = sum2 = 0
  for ind_k = (half_K+1):K
    sum1 = sum1 + ( 1 - (ind_k-half_K) / (half_K+1) )^8
    Set_for_j = collect((half_K+1):K)
    Set_for_j = setdiff(Set_for_j,ind_k)   # setdiff is a R function for elements in the first term, not in the seond term. Change that in Matlab
    for (ind_j in Set_for_j) {
      sum2 = sum2 + ( 1 - (ind_k - half_K) / ( half_K + 1 ) )^4 * ( 1 - (ind_j - half_K) / (half_K+1) )^4 * rho_z
    end
  end

  c_K_irr = sqrt( ( R_sq_f / ( 1-R_sq_f ) ) * ( 1 / (sum1 + sum2) ) )

  K_vec = collect(1:K)

  pi_irr = c_K_irr * ( 1 - (K_vec - half_K ) / ( half_K + 1 ) )^4 * ( K_vec > half_K )

  return pi_irr
end


# dgp_LS -------------------------------------
# PURPOSE:
#     Generate data for simulation in Lee and Shin
#
# USAGE:
#     dgp_LS(N,K,rho_eu,rho_z,pi_0,bt_0)
#
# MODEL:
#
#     Structural Eq:
#           y = x * beta + e
#     1st Stage Eq:
#           x = z' * pi  + u
#
# INPUT:
#     N: a scalar; the sample sizee
#     K: a scalar; the number of instruments
#     rho: a scalar; Cov(e, u)
#                            | 1        rho_eu |
#     Cov(e,u) = Sigma_eu =  | rho_eu   1      |
#      2 x 2
#
#     rho_z: a scalar; parameter for the correlation of z
#                          |   1       rho_z  ___   rho_z  |
#     Var(Z) = Sigma_Z  =  |   rho_z     1    ___   rho_z  |
#    ( K x K )             |   ___     ___    ___   ___    |
#                          |   rho_z   rho_z  rho_z  1     |
#     pi_0: K x 1 vector; pi
#     bt_0: a scalar; beta
#     high: a scalar in [0,100] ; higher percentile
#     low: a scalar in [0,100] ; lower percentile
#
# OUPUT:
#     y: N x 1; a dependenvt variable
#     x_end: N x 1; an endogenous variable
#     z: N x K; instruments
#
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/11/17      Y_ shin             Original Code
dgp_LS = function(N, K, rho_eu, rho_z, pi_0, bt_0)
end



# Projection matrix -----------------------------------------------------------
# USAGE:
#     proj_m(z)
#
# PURPOSE:
#     Calculate the projection matrix of z
#
# INPUT:
#     z: N x K, a matrix
#
# OUPUT:
#     proj_z: N x n; the projection matrix of z
#            proj_z = z(z'z)^{-1}z'
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/11/17      Y_ shin             Original Code
function proj_m(z)
  #library('gpuR')
  #z = gpuMatrix(z)
  #inv_z = (z * as_matrix( inv( t(z) * z ) ) * t(z))
  #return(as_matrix(inv_z))
  z * ( inv( z'*z ) ) * z'
end

# function proj2_m(z_in)
#   library('gpuR')
#   z = vclMatrix(z_in)
#   as_matrix((z * ( inv( t(z) * z ) ) * t(z)))
# end


# function trace(X)
#   sum(diag(X))
# end
#
# function matmul(X)
#   X * X
# end

# ols function -----------------------------------------------------------
# USAGE:
#     ls(y,x)
#
# PURPOSE:
#     Calculate the least squares estimator bt_hat for
#         y = x' bt + eps
#
# INPUT:
#     y: N x 1; a dependent variable
#     x: N x K; regressors
#
# OUPUT:
#     bt_hat: K x 1, the estimate for bt
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/11/17      Y_ shin             Original Code
# ols = function(y,x){
#   bt_hat = inv( x' * x ) * x' * y
#   return bt_hat
# }



# MSE -----------------------------------------------------------
# USAGE:
#     mb(ap_hat, ap_0)
#
# PURPOSE:
#     Calculate the mean squared error of the estimator
#
# INPUT:
#     ap_hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#     ap_0: (d x 1) vector; the true value of ap
#
# OUPUT:
#     mse: (d x 1) vector; mse of each parameter values
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/11/17      Y_ shin             Original Code
function mse( ap_hat, ap_0 )
  # Change it into a matrix just in case it is a numeric vector
  # ap_hat = as_matrix(ap_hat)

  # Define constants
  R = nrow(ap_hat)   # The number of replications
  d = ncol(ap_hat)   # The number of parameters

  # Construct a (Rxd) matrix of true parameter values
  # ones = ones(R,1)               # (Rx1) vector of 1's
  # mat_ap_0 = t(ap_0) %x% ones   # (Rxd) matrix of true values
  ap_0 = repeat(ap_0, R , 1)
  # Calculate the squared error
  sq_err = ( ap_hat - ap_0 ) .^ 2

  # Calculate the mean squared error
  mse = mean(sq_err)
  return mse
end


# Median Absolute Deviation (MAD)-----------------------------------------------------------
# USAGE:
#     mad(ap_hat)
#
# PURPOSE:
#     Calculate the median abolute deviation
#
#
# INPUT:
#     ap_hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#
#
# OUPUT:
#     mad: (d x 1) vector; median of asolute deviation distribution
#
# REVISION LOG:
#     Date        Programmer          Description
#   2017-02-11    Y_ shin             Original Code
#   2018-11-22
function mad( ap_hat )

  # Change it into a matrix just in case it is a numeric vector
  ap_hat = as_matrix(ap_hat)

  # Define constants
  R = nrow(ap_hat)   # The number of replications
  d = ncol(ap_hat)   # The number of parameters
  med_ap_hat = apply(ap_hat,2,median)

  # Construct a (Rxd) matrix of true parameter values
  ones = rep(1,R)                     # (Rx1) vector of 1's
  mat_ap_0 = t(med_ap_hat) %x% ones   # (Rxd) matrix of true values

  # Calculate the absolute values of the bias
  abs_bias = abs(ap_hat-mat_ap_0)

  # Calculate te median of the abolute biases
  median_absolute_dev = apply(abs_bias,2,median)
  return(as_numeric( median_absolute_dev ))
end


# inter-percentile range -----------------------------------------------------------
# USAGE:
#     range(ap_hat, high, low)
#
# PURPOSE:
#     Calculate inter-percentile range of the estimator
#
# INPUT:
#     ap_hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#     high: a scalar in [0,100] ; higher percentile
#     low: a scalar in [0,100] ; lower percentile
#
# OUPUT:
#     inter_per: (d x 1) vector; (high-low) percentile range
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/11/17      Y_ shin             Original Code

function range(ap_hat, high, low)
  # Change it into a matrix just in case it is a numeric vector
  ap_hat = as_matrix(ap_hat)


  # Calculate the percentiles at high and low
  high_per = apply(ap_hat, 2, quantile, probs=high/100)
  low_per = apply(ap_hat, 2, quantile, probs=low/100)

  # Calculate and return the inter percentile
  inter_per = high_per - low_per
  return inter_per
end


# Median Bias -----------------------------------------------------------
# USAGE:
#     mb(ap_hat, ap_0)
#
# PURPOSE:
#     Calculate the median of the bias distribution
#
# INPUT:
#     ap_hat: (R x d) matrix; the estimates in each replication
#             where R is the number of replication
#                   d is the size of estimator
#     ap_0: (d x 1) vector; the true value of ap
#
# OUPUT:
#     mb: (d x 1) vector; median of bias distribution
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/11/17      Y_ shin             Original Code
function mb( ap_hat, ap_0 )
  # Change it into a matrix just in case it is a numeric vector
  ap_hat = as_matrix(ap_hat)

  # Define constants
  R = nrow(ap_hat)   # The number of replications
  d = ncol(ap_hat)   # The number of parameters

  # Construct a (Rxd) matrix of true parameter values
  ones = rep(1,R)               # (Rx1) vector of 1's
  mat_ap_0 = t(ap_0) %x% ones   # (Rxd) matrix of true values

  # Calculat the bias
  bias = ap_hat-ap_0

  # Calculate the median bias and return the value
  mb = apply(bias, 2, median)
  return mb
end

# Bias --------------------------------------
bias = function bias( ap_hat, ap_0 )
  # Change it into a matrix just in case it is a numeric vector
  ap_hat = as_matrix(ap_hat)

  # Define constants
  R = nrow(ap_hat)   # The number of replications
  d = ncol(ap_hat)   # The number of parameters

  # Construct a (Rxd) matrix of true parameter values
  ones = rep(1,R)               # (Rx1) vector of 1's
  mat_ap_0 = t(ap_0) %x% ones   # (Rxd) matrix of true values

  # Calculat the bias
  bias = ap_hat-ap_0

  # Calculate the median bias and return the value
  bias = apply(bias, 2, mean)
  return bias
end

# 2SLS -----------------------------------------------------------
# USAGE:
#     tsls_est(y,x_end,x_exo,z)
#
# PURPOSE:
#     Calculate the two-stage least squares estimator bt_hat for
#         y = x' bt + eps
#         where x = ( x_end, x_exo )'
#
# INPUT:
#     y: N x 1, a dependent variable
#     x_end: N x d1, endogenous regressors
#     x_exo: N x d2, exogenous regressors
#     z: N x K; instrudments excluding x_exo
#
# OUPUT:
#     bt_hat: (d1+d2) x 1, the estimate for bt
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/11/17      Y_ shin             Original Code
function tsls_est(y, x_end, x_exo, z)

  # Combine z and x_exo if x_exo exists
  if x_exo == nothing                 #( !is_null(x_exo) )
      z_big = hcat( ones(size(x_end, 1), 1), z )
      x = x_end
  else
      z_big = hcat( ones(length(x_end), 1), x_exo, z )
      x = hcat(x_end, x_exo)
  end

  # Construct a projection matrix with z_big
  p = proj_m( z_big )
  # Calculate 2SLS estimator for b_hat
  bt_hat = inv( x' * p * x ) * x' * p * y
  return bt_hat
end


function pre_est(method_pre, y, x_end, x_exo, z_pre, P_all, ld)   # Only "one_step" and "two_step" are available now. Not "lasso"
  # cat('---------------------- \n')
  # cat('Preliminary Estimation \n')
  # cat('---------------------- \n')
  # cat('Method = ', method_pre, '\n')

  if x_exo == nothing
      x = x_end
      d_2 = 0
      iv_pre = z_pre
  else
      x = hcat(x_end, x_exo)
      d_2 = size(x_exo, 1)
      iv_pre = hcat(x_exo, z_pre)
  end

  d_x = size(x, 2)
  # dim of included exogenous regressors
  # d_2 =ncol(as_matrix(x_exo))
  N = size(x, 1)
  K_excl = size(z_pre, 2) - d_2
  if ld == nothing
  ld = ones(d_x, 1)./d_x
  end

  if method_pre == "one_step"
    P_pre = proj_m(iv_pre)
    f_pre = P_pre * x
    H_pre = ( x' * P_pre * x ) / N
    H_inv = inv(H_pre)
    u_pre = ( Diagonal(vec(ones( N, 1 ) ) ) - P_pre ) * x
    u_ld_pre = u_pre * H_inv * ld

    bt_pre = tsls_est( y, x_end, x_exo, z_pre )
    eps_pre = y - x * bt_pre
    sig2_eps =  ( eps_pre' * eps_pre ) / N
    sig_ld_eps = ( u_ld_pre' * eps_pre ) / N
    sig2_ld = (u_ld_pre' * u_ld_pre) / N
    sig_u_eps = (u_pre' * eps_pre) / N

    # cat('sig2_eps   = ', sig2_eps, '\n')
    # cat('sig2_ld    = ', sig2_ld, '\n')
    # cat('sig_ld_eps = ', sig_ld_eps, '\n')


     elseif (method_pre == "two_step")
    P_pre = proj_m(iv_pre)
    H_pre = ( x' * P_pre * x ) / N
    H_inv = inv(H_pre)
    u_pre = ( Diagonal(vec(ones( N, 1 ) ) ) - P_pre ) * x
    u_ld_pre = u_pre * H_inv * ld
    sig2_ld = ( u_ld_pre' * u_ld_pre ) / N

    Mallows_k = zeros(K_excl, 1)
    for  i_pre = 1 : K_excl
      P_k = P_all[ :, :, i_pre]
      u_pre_k = ( Diagonal(vec(ones( N, 1 ) ) ) - P_k ) * x
      u_ld_k = u_pre_k * H_inv * ld
      Mallows_k[ :,i_pre ] = (  u_ld_k' * u_ld_k ) / N + sig2_ld * ( 2  *i_pre / N )
    end
    k_m = argmin(Mallows_k) # Choice of k after applying the Mallows criterion
    k_m = k_m[2]
    # cat('k_m          =', k_m, '\n')

    # Re-estimate, eps, u, H with k_m
    z_m = z[ : , ( 1 : (d_2 + k_m ) ) ]
    P_m = P_all[: , : , k_m]
    f_pre = P_m * x
    H_m = ( x' * P_m * x ) / N
    H_inv = inv(H_m)
    u_m = u_pre = (Diagonal(vec(ones( N, 1 ) ) ) - P_m ) * x
    u_ld_m = u_m * H_inv * ld

    bt_m = tsls_est( y, x_end, x_exo, z_m[ : , (d_2 + 1): end] )    # -c( 1 : d_2 ) ] )
    eps_m = y - x * bt_m
    sig2_eps =  ( eps_m' * eps_m ) / N
    sig_ld_eps =  ( u_ld_m' * eps_m ) / N
    sig2_ld = u_ld_m' * u_ld_m / N
    sig_u_eps = ( u_m' * eps_m ) / N

    # cat('sig2_eps   = ', sig2_eps, '\n')
    # cat('sig_ld    = ', sig2_ld, '\n')
    # cat('sig_ld_eps = ', sig_ld_eps, '\n')

  end
   # else  if (method_pre == 'lasso'){
   #  library( 'glmnet' )
   #  if ( ncol( as_matrix( x_end ) == 1 ) ){
   #    cv_m = cv_glmnet( x = z, y = x_end )
   #    ld_1se = cv_m$lambda_1se
   #    m_first_stage = glmnet( x = z, y = x_end, family = 'gaussian', intercept = T, penalty_factor = c( rep( 0, d_x ), rep( 1, K_excl ) ), lambda = ld_1se )
   #    pi_1st = coef( m_first_stage )[ - 2]
   #    selected = which( pi_1st!= 0 )
   #    cat('1st Stage Lasso Selection =', selected, '\n')
   #    cat('index for x_exo           =', ( 1 : d_2 ), '\n')
   #    cat('index for z_excl           =', ( ( d_2 + 1 ) : ( ncol(z) ) ), '\n' )
   #    cat('No of selected z          =', length( selected ), '\n' )
   #
   #    # Calculate eps, u, H
   #    z_sel = z[ , selected ]
   #    P_pre = proj_m(z_sel)
   #    f_pre = P_pre * x
   #    H_pre = ( x' * P_pre * x ) / N
   #    H_inv = inv(H_pre)
   #    u_pre = ( Diagonal(vec(ones( N, 1 ) ) ) - P_pre ) * x
   #    u_ld_pre = u_pre * H_inv * ld
   #
   #    bt_pre = tsls_est( y, x_end, x_exo, z = z_sel[ , - c( 1 : d_2 ) ] )
   #    eps_pre = y - x * bt_pre
   #    sig2_eps = drop( ( eps_pre' * eps_pre ) / N )
   #    sig_ld_eps = drop( ( u_ld_pre' * eps_pre ) / N )
   #    sig2_ld = (u_ld_pre' * u_ld_pre ) / N
   #    sig_u_eps = ( u_pre' * eps_pre ) / N
   #
   #    # Evaluate the obj_ fun for different k
   #    cat('sig2_eps   = ', sig2_eps, '\n')
   #    cat('sig2_ld    = ', sig2_ld, '\n')
   #    cat('sig_ld_eps = ', sig_ld_eps, '\n')
   #  } else {
   #    cat('To be added \n')
   #    break
   #  }

 return (sig2_eps = sig2_eps,
              sig2_ld = sig2_ld,
              sig_ld_eps = sig_ld_eps,
              sig_u_eps = sig_u_eps,
              H_inv = H_inv,
              u_pre = u_pre,
              f_pre = f_pre  )
end

function gen_p_triangle(K)  # Generate Pascal's triangle, the lower triangular matrix with combinations of row per column (if row is k and column is n, it ic C(n, k))
  p_triangle =  zeros(Int64, K, K)
  for i = 1:K
    for j = 1:i
      p_triangle[i, j] = binomial(i, j)     # Combinations, à la choose(i,j) in R
    end
  end
  return p_triangle
end

function tsls_CSA_get_P(y, x_end, z_excl, sub_K; x_exo = nothing, R = Inf, ld = nothing)
  # sub_k is the size of the instruments used
  if  x_exo  == nothing
      d_2 = 0
      z = z_excl
      x = x_end
  else
      d_2 = size( x_exo, 2 )
      z = hcat(x_exo, z_excl)
      x = hcat(x_end, x_exo)
  end
  d_x = size(x, 2)

 # dim of included exogenous regressors
  # d_2 =ncol(as_matrix(x_exo))
  d_x = size(x, 2)
  N = size(z, 1)
  K = size(z, 2)

  # If lambda (ld) is not defined we set it to 1/dim(x), i.e. equal weights.
  if ld == nothing
    ld = ones( 1, d_x)./d_x
  end

  # Generate P(N,N,K) array, all projection matrices with averaging
  # Calculate u_hat_k
  K_excl = K - d_2
  aveP =  zeros(N, N, K_excl)             # replicate(K_excl, matrix(0,N,N), simplify = F)
  n_subset = gen_p_triangle(K_excl)[K_excl, :]                             # gen_p_triangle(K_excl)[K_excl, ]

 #   if use_par == TRUE
 #     library(doMC)
 #     library(foreach)
 #     registerDoMC(cores = n_core)
 #     aveP_par= (foreach(i_aveP = 1:K_excl) %dopar% loop_for_aveP(i_aveP, x_exo, z_excl, R, n_subset, K_excl))
 #     return(aveP_par)
 #   }
 # else{
    for i_aveP = 1:sub_K
      # cat('No_ of IVs =', i_aveP, '\n')
      # cat('No_ of Subsets =', n_subset[i_aveP], '\n')
      if n_subset[i_aveP] <= R
        sum_sub = zeros(N, N)
        subset_ind = collect(combinations(1:K_excl, i_aveP) )
        for  i_subset = 1:n_subset[i_aveP]
          z_sub = z_excl[: , subset_ind[i_subset]]
          if x_exo == nothing
              proj_sub = proj_m(z_sub)
          else
              proj_sub = proj_m(hcat(z_sub, x_exo))
          end
          # proj_sub = proj_m(hcat(z_sub, x_exo))
          sum_sub = sum_sub + proj_sub
        end
        aveP[:, :, i_aveP] = sum_sub ./ n_subset[i_aveP]

      else   # if the number of subset if bigger than R, then we conduct only R random selections
        # cat('Random Selection', '\n')
        sum_sub = zeros(N, N)
        for i_sub = 1:R
          sub_rnd = sort(Random.randperm(K_excl)[1:i_aveP])                             # sort(sample(1:K_excl, i_aveP, replace=F))
          z_sub = z_excl[: , sub_rnd]
          if x_exo == nothing
              proj_sub = proj_m(z_sub)
          else
              proj_sub = proj_m(hcat(z_sub, x_exo))
          end
          sum_sub = sum_sub + proj_sub
        end
        aveP[:, :, i_aveP] = sum_sub ./ R
      end
    end
    return aveP
end

 function tsls_CSA_aveP_in(y, x_end, z_excl, aveP, method_pre, z_pre; x_exo, ld = nothing)
  # if(missing(x_exo)) {x_exo = NULL}
  # Declare constants
  if x_exo == nothing
      z = z_excl
      x = x_end
      d_2 = 0
      iv_pre = z_pre
  else
      z = hcat(x_exo,z_excl)
      x = hcat(x_end, x_exo)
      d_2 = size(x_exo, 1)
      iv_pre = hcat(x_exo, z_pre)
  end
  N = size(z, 1)
  K = size(z, 2)
  K_excl = size(z_excl, 2)
  # dim of included exogenous regressors
  # d_2 =ncol(as_matrix(x_exo))
  d_x = size(x, 2)

  # If lambda (ld) is not defined and use the equal weight
  if ld == nothing
    ld = ones(d_x, 1)./d_x
  end

  # Construct projection matrices of each model and its weighted averages_
  P_all = zeros(N, N, K_excl)                        # array( NA, c(N, N, K_excl) )
  for i = 1 : K_excl
    P_m = proj_m( z[: , (1:(d_2 + i) ) ] )
    P_all[ : , : , i] = P_m
  end

  # cat('ld = ', ld, '\n')
  # cat('d_x = ', d_x, '\n')
  # cat('ld = ', length(ld), '\n')


  # Preliminary Regression to get sig2's
  prelim = pre_est(method_pre, y, x_end, x_exo, z_pre, P_all, ld)

  sig2_eps = prelim.sig2_eps[1]
  sig2_ld = prelim.sig2_ld
  sig_ld_eps = prelim.sig_ld_eps[1]
  sig_u_eps = prelim.sig_u_eps
  H_inv = prelim.H_inv
  H_ld = H_inv * ld
  u_pre = prelim.u_pre
  f_pre = prelim.f_pre
  P_f_pre = proj_m(f_pre)

  # Calculate S_hat_ld_k for k=1,___, K
  S_hat_ld_k = zeros(1, K_excl - 1)       #array(0,K_excl-1)
  bias_k = zeros(1, K_excl - 1)
  var_k = zeros(1, K_excl - 1)

  for  j = 1:(K_excl-1)
    P_K = aveP[ :, :, j]
    I_N = Diagonal(vec(ones( N, 1 ) ) )
    bias_k[j] = (sig_ld_eps^2 * j^2/N)[1]

    #--- Calculate S_hat function
    Sigma_u_tilde = (u_pre' * u_pre ) / N
    S_term_1 = (x' * (I_N - P_K)* (I_N - P_K)  * x) / N - Sigma_u_tilde * (-2*j/N + sum(diag(P_K*P_K))/N)
    S_term_2 = (x' * (I_N - P_K)  * x) / N + Sigma_u_tilde * ( (j/N) - 1 )
    mid_term = S_term_1 - S_term_2 * H_inv * S_term_2
    var_k[j] = (sig2_eps * ( H_ld' *  mid_term * H_ld ))[1]
    S_hat_ld_k[j] = bias_k[j] + var_k[j]
  end
  opt_k = argmin(S_hat_ld_k)[2]
  # cat('opt_k =', opt_k, '\n')

  bt_hat_CSA = inv( x' * aveP[ :, :, opt_k] * x ) * ( x' * aveP[ :, :, opt_k] *y )

  # cat('bt_hat_CSA = ',bt_hat_CSA,'\n')
  return (    bt_hat = bt_hat_CSA,
              opt_k=opt_k,
              obj_v=S_hat_ld_k,
              aveP=aveP)
end

# CSA estimator-----------------------------------------------------------
# USAGE:
#     tsls_CSA(y, x_end, z_excl, method_pre, z_pre; x_exo = nothing, R = Inf, ld = nothing )
#
# PURPOSE:
#     Calculate the two-stage least squares estimator using subsets of instruments for
#         y = x' bt + eps
#         where x = ( x_end, x_exo )'
#
# INPUT:
#     y: N x 1, a dependent variable
#     x_end: N x d1, endogenous regressors
#     x_exo: N x d2, exogenous regressors
#     z_excl: N x K, instruments excluding x_exo
#
# OUPUT:
#     bt_hat: (d1+d2) x 1, the estimate for bt
#     opt_k: the optimal subset size
#     obj_v: the values of the objective functions
#     aveP: the average projection matrices for each subset size
# REVISION LOG:
#     Date        Programmer          Description
#    02/11/17      Y_ shin             Original Code
#    08/07/19      Thomas Vigié        Code translated in Julia    ]


function tsls_CSA(y, x_end, z_excl, method_pre, z_pre; x_exo = nothing, R = Inf, ld = nothing)
   if x_exo == nothing
       d_2 = 0
       z = z_excl
       x = x_end
       iv_pre = z_pre
   else
       d_2 = size( x_exo, 2 )
       z = hcat(x_exo,z_excl)
       x = hcat(x_end, x_exo)
       iv_pre = hcat(x_exo, z_pre)
   end
  # Declare constants

  N = size(z, 1)
  K = size(z, 2)
  K_excl = size(z_excl, 2)
 # dim of included exogenous regressors
  # d_2 =ncol(as_matrix(x_exo))
  d_x = size(x, 2)

  # If lambda (ld) is not defined and use the equal weight
  if ld == nothing
    ld = ones(d_x, 1)./d_x
  end

  # Construct projection matrices of each model and its weighted averages_
  P_all = zeros(N, N, K_excl)# array( NA, c(N, N, K_excl) )
  for  i =  1:K_excl
    P_m = proj_m(z[ :, (1: (d_2 + i ) ) ] )
    P_all[ :, :, i] = P_m
  end

  # Generate all projection matrices with averaging: list of K matrices whose dimension is (N by N)
  aveP = tsls_CSA_get_P(y, x_end, z_excl, K_excl, x_exo = x_exo, R = Inf, ld = ld)
  # The CSA-2SLS estimator with aveP
  m_csa = tsls_CSA_aveP_in(y, x_end, z_excl, aveP, method_pre, z_pre, x_exo = x_exo, ld = ld,)

  return (bt_hat = m_csa.bt_hat,
              opt_k = m_csa.opt_k,
              obj_v = m_csa.obj_v,
              aveP = aveP )
end

function tsls_DN(y, x_end, z_excl, method_pre, z_pre; x_exo = nothing, ld = nothing)
 # Required functions
 if x_exo == nothing
     d_2 = 0
     z = z_excl
     x = x_end
     iv_pre = z_pre
 else
     d_2 = size( x_exo, 2 )
     z = hcat(x_exo,z_excl)
     x = hcat(x_end, x_exo)
     iv_pre = hcat(x_exo, z_pre)
 end
 # Declare constants

 N = size(z, 1)
 K = size(z, 2)
 K_excl = size(z_excl, 2)
 d_x = size(x, 2)

 # If lambda (ld) is not defined and use the equal weight
 if ld == nothing
   ld = ones(d_x, 1)./d_x
 end

 # Collect all projection matrices
 P_all = zeros(N, N, K_excl)
 # array( NA, c(N, N, K_excl) )
 for i = 1 : K_excl
   P_m = proj_m( z[:, (1:(d_2+i) ) ] )
   P_all[: ,: , i] = P_m
 end


 prelim = pre_est(method_pre, y, x_end, x_exo, z_pre, P_all, ld)
 sig2_eps = prelim.sig2_eps[1]
 sig2_ld = prelim.sig2_ld
 sig_ld_eps = prelim.sig_ld_eps
 H_inv = prelim.H_inv

 # Evaluate the obj_ fun for different k
 S_DN_k = zeros(1, K_excl)    #rep(NA, K_excl)
 for i_DN = 1:K_excl
   P_k = P_all[ :, :, i_DN]
   u_pre_k = ( Diagonal(vec(ones( N, 1 ) ) ) - P_k ) * x
   u_ld_k = u_pre_k * H_inv * ld
   S_DN_k[i_DN] = (sig_ld_eps.^2 * (i_DN^2 / N) + sig2_eps^2 * ( (u_ld_k' * u_ld_k) / N + sig2_ld .*(i_DN/N) ) )[1]
 end

 k_opt = argmin(S_DN_k)[2]
 z_opt = z[ :, (1 : (d_2+k_opt) ) ]
 P_DN = proj_m( z_opt )
 if d_2 == 0
   bt_DN = tsls_est(y, x_end, x_exo, z_opt )
 else
   bt_DN = tsls_est(y, x_end, x_exo, z_opt[ :, d_2 : end ] )
 end
 # cat('k_opt = ', k_opt,'\n')
 # cat('bt_hat_DN = \n', bt_DN,'\n')
 return( opt_k = k_opt,
         bt_hat = bt_DN,
         P_DN = P_DN  )
end
# tsls_DN -----------------------------------------------------------
# USAGE:
#     tsls_DN(y, x_end, z_excl, method_pre, z_pre; x_exo = nothing, ld = nothing )
#
# PURPOSE:
#     Calculate the two-stage least squares estimator bt_hat for
#         y = x' bt + eps
#         where x = ( x_end, x_exo )' with selection of the instruments according to Donald and Newey (2001)
#
# INPUT:
#     y: N x 1, a dependent variable
#     x_end: N x d1, endogenous regressors
#     z_excl: N x K, the instruments
#     x_exo: N x d2, exogenous regressors  (default to nothing)
#     ld: 1 x d1 vector of weights for the endogenous variables, if more than one (default to nothing, i.e. equal weight on all the endogenous variables)
#
# OUPUT:
#     bt_hat: d x 1, the estimate for bt
#     opt_k: a scalar, the optimal K_star_hat, estimated value, arg_max S_hat(K)
#     P_DN: the projection matrix of the selected instruments

# REVISION LOG:
#     Date        Programmer          Description
#  02/13/2017     Y_ Shin             Original Code
#  02/16/2017     Y_ Shin             For R_hat_ld_K, use updated sig_ld instead of sig_hat_ld_pre (for Mallows)
#  08/27/2019     Thomas Vigié        Function coded in Julia


# Example:
using CSV
using IterTools
using Combinatorics
include("auxiliary functions.jl")
data = CSV.File("test data.csv")
psiz = my_poly(data.z, order = 5)
gna = tsls_CSA(data.y, data.x, psiz, "one_step", data.z; x_exo = nothing, R = Inf, ld = nothing)
tsls_CSA(data.y, data.x, psiz, "two_step", data.z; x_exo = nothing, R = Inf, ld = nothing)
gna = tsls_DN(data.y, data.x, psiz, "one_step", data.z, x_exo = nothing, ld = nothing)



# -------------------------------------------------------------------

loop_for_aveP = function(i_aveP, x_exo, z_excl, R, n_subset,  K_excl){

  N = size(z_excl, 1)

  cat('No_ of IVs =', i_aveP, '\n')
  cat('No_ of Subsets =', n_subset[i_aveP], '\n')
  if (n_subset[i_aveP] <= R){
    sum_sub = matrix(0, N, N)
    subset_ind = combn(K_excl, i_aveP)
    for ( i_subset in c(1:n_subset[i_aveP]) ){
      z_sub = z_excl[,subset_ind[,i_subset]]
      proj_sub = proj_m(cbind(x_exo, z_sub))
      sum_sub = sum_sub + proj_sub
    }
    aveP_k = sum_sub / n_subset[i_aveP]
  } else {  # if the number of subset if bigger than R, then we conduct only R radom selections
    cat('Random Selection', '\n')
    sum_sub = matrix(0, N, N)
    for (i_sub in (1:R)){
      sub_rnd = sort(sample(1:K_excl, i_aveP, replace=F))
      z_sub = z_excl[,sub_rnd]
      proj_sub = proj_m(cbind(x_exo, z_sub))
      sum_sub = sum_sub + proj_sub
    }
    aveP_k = sum_sub / R
  }
  return aveP_k
end




# CV objective function -----------------------------------------------------------
# USAGE:
#     R_hat_cv( K, x, z, H_tilde )
#
# PURPOSE:
#     Calculate the objective function value of CV
#
# INPUT:
#     K: the number of IVs where the ftn will be evaluated
#     x: N x 1, endogenous variable
#     z: N x K; instrudments excluding e_exo
#
# OUPUT:
#     the evaluated funcation value
#
# REVISION LOG:
#     Date        Programmer          Description
#    2/12/17      Y_ shin             Original Code
#    [Consider a multidensional delta and add 'ld' laster]
#
function R_hat_cv( K, x, z, H_tilde)

  # Declare some constants
  # z = as_matrix(z)
  N = size(z, 1)

  z_K = z[,(1:K)]
  P_K = proj_m(z_K)
  u_hat_K = ( Diagonal(vec(ones( N, 1 ) ) ) - P_K ) * x
  u_hat_ld_K = u_hat_K * inv(H_tilde)

  R_hat_cv_value=mean( u_hat_ld_K^2/((1-diag(P_K))^2) )

  return R_hat_cv_value
end



# Mallows objective function -----------------------------------------------------------
# USAGE:
#     R_hat_m( K, x, z, H_tilde, sig_hat_ld )
#
# PURPOSE:
#     Calculate the objective function value of the Mallows criterion
#
# INPUT:
#     K: the number of IVs where the ftn will be evaluated
#     x: N x 1, endogenous variable
#     z: N x K; instrudments excluding e_exo
#
# OUPUT:
#     the evaluated funcation value
#
# REVISION LOG:
#     Date        Programmer          Description
#   2/12/2017     Y_ Shin             Original Code
#    [Consider a multidensional delta and add 'ld' laster]
#
 function R_hat_mallows( K, x, z, H_tilde, sig_hat_ld, ld )
  # Declare constants
  # z = as_matrix(z)
  N = size(z, 1)

  z_K = z[,(1:K)]
  P_K = proj_m(z_K)
  u_hat_K = (Diagonal(vec(ones( N, 1 ) ) )-P_K) * x
  u_hat_ld_K = u_hat_K * inv(H_tilde) * ld

  R_hat_m_value = ( ( t(u_hat_ld_K) * u_hat_ld_K ) / N ) + ( ( 2 * sig_hat_ld * K )/ N )
  return R_hat_m_value
end


# 1st stage GOF criterion function -----------------------------------------------------------
# USAGE:
#     first_stage_GOF( y, x, z, criteria )
#     # first_stage_GOF( y, x_end, x_exo, z, criteria )
#
# PURPOSE:
#     Calculate the preliminary number of instruments K_opt
#
# INPUT:
#     y: N x 1, a dependent variable
#     x: N x 1, endogenous regressors
#     # x_end: N x d1, endogenous regressors
#     # x_exo: N x d2, exogenous regressors
#     z: N x K; instrudments excluding e_exo
#     H_tilde: d x d, matrix, for fixed K_
#     criteria: characters, 'CV' (cross-validation) or 'Mallows'
#
# OUPUT:
#     opt_K: a scalar, the preliminary number of instruments
#     obj_value: a scalar, the objective function values_ The 'opt_K'
#           is the minimizer of this function valuef
#
# REVISION LOG:
#     Date        Programmer          Description
#   2017-02-12    Y_ shin             Original Code
#   2017-10-17                        Knitro Part error corrected (get_weight)
#
#    [Consider a multidensional delta and add 'ld' laster]
#
 function first_stage_GOF( y, x, z, H_tilde, sig_hat_ld, criteria, ld, d_2 )
  # z = as_matrix(z)  # z matrix should be [x_exo, z_excl]
  N = size(z, 1)
  K = size(z, 2)

  obj_value = rep(NA,K)

  # When the criteria is Cross-validations
  if ( criteria == 'CV' ){
    for ( i in (1:K) ){
      obj_value[i] = R_hat_cv( K=i, x=x, z=z, H_tilde=H_tilde )
    }
  } else if ( criteria == 'Mallows' ) {
    for (i in (1:K)){
      obj_value[i] = R_hat_mallows( K=i, x=x, z=z, H_tilde=H_tilde, sig_hat_ld=sig_hat_ld, ld=ld )
    }
  }
  opt_K = which_min(obj_value[-c(1:d_2)]) + d_2   # Find the minimum value outside (1:d_2)

  return( list( opt_K=opt_K, obj_value=obj_value ) )
end


# This functions returns estimates for residuals required for estimating a MSE estimate ---------------------------------
 function get_sigma_DN( y, x, z, opt_K, H_tilde, ld )
  N = size(y, 1)
  Z_1st = z[:, (1:opt_K)]
  P_1st = proj_m(Z_1st)
  bt_1st = tsls_est(y=y, x_end=x, x_exo=NULL, z=Z_1st)
  eps_tilde = y - x * bt_1st
  u_tilde = ( Diagonal(vec(ones(N, 1))) - P_1st ) * x
  u_ld = u_tilde* inv(H_tilde) * ld
  sig_eps2 = ( eps_tilde'* eps_tilde ) / N
  sig_ld2 = ( u_ld'* u_ld ) / N
  sig_ld_eps2 = ( u_ld'*eps_tilde ) / N

  return( list( eps2=sig_eps2, ld2=sig_ld2, ld_eps2=sig_ld_eps2, bt_1st=bt_1st ) )
end





#Functions required for KO -----------------------------------------------------------------------------
get_weights = function(y, x_end, x_exo, z_excl, H_tilde, sig_hat_ld, P_all, Gm, initial_z='Mallows', ld){
  library("quadprog")

  # Declare constants
  z = hcat(x_exo, z_excl)
  K = size(z, 2)
  N = size(z, 1)
  d_2 = size(x_exo, 2)

  x = hcat(x_end, x_exo)
  opt_w = rep(NA, K)

  # To get the preliminary estimate bt_tilde
  #r1 = first_stage_GOF(y=y, x=x, z=z, H_tilde=H_tilde, sig_hat_ld=sig_hat_ld, criteria='Mallows')
  #r1 = first_stage_GOF( y=y, x=x, z=z, H_tilde=H_tilde, sig_hat_ld=sig_hat_ld, criteria='Mallows', ld=ld, d_2=d_2 )

  #K_tilde = r1$opt_K
  #if (GF=='all') {K_tilde=K}
  #cat('KO estimator: K_tilde =', K_tilde, '\n')

  # Get preliminary K_opt using the first stage Goodness of fit criteria
  if (initial_z == 'all') {
    K_tilde=K
  } else {
    r1 = first_stage_GOF( y=y, x=x, z=z, H_tilde=H_tilde, sig_hat_ld=sig_hat_ld, criteria=initial_z, ld=ld, d_2=d_2 )
    K_tilde = r1$opt_K
  }
  cat('KO estimator: K_tilde = ',K_tilde,'\n')


  # Update H_tilde using K_tilde
  H_tilde_update = ( t(x) * P_all[ , , K_tilde] * x ) / N

  # Calculate sig_hat_eps, sig_hat_ld, sig_hat_ld_eps, H_hat using 'get_sigma_DN'
  #r2 = get_sigma_DN(y=y ,x=x ,z=z, opt_K=K_tilde, H_tilde=H_tilde_update)
  r2 = get_sigma_DN(y=y, x=x, z=z, opt_K=K_tilde, H_tilde=H_tilde_update, ld)

  sig_hat_eps = as_numeric(r2$eps2)
  sig_hat_ld = as_numeric(r2$ld2)
  sig_hat_ld_eps = as_numeric(r2$ld_eps2)


  # Calculate u_hat_ld and U_hat
  u_stack = matrix(0, N, K)
  for i = 1:K
    u_stack[ , i] = ( P_all[ , , K] - P_all[ , , i] ) * x * inv(H_tilde_update) * ld
  end
  U_hat = u_stack'*u_stack              # (K x K) matrix for S(K)



  # Generate a (K x 1) vector K_vec s_t_ K_vec[i]=i
  K_vec = c(1:K)

  # Generate a (K x 1) vector of 1's
  one_K = rep(1,K)

  # Generate a (K x 1) vector of 0's
  zero_K = rep(0,K)


  #  This expression comes from Equation (2_6) in KO_
  #  c = as_numeric(2 * sig_hat_eps * sig_hat_ld * t(K_vec)  ) /N
  #  H = 2 * ( sig_hat_ld_eps^2 * (K_vec * t(K_vec)  ) + (sig_hat_eps * U_hat  ) - (sig_hat_eps * sig_hat_ld * Gm )  ) / N

  # Equation (2_5) is rewritten: H is a matrix for the quadratic term; c is a matrix for the linear term
  c = as_numeric(-(B_ld_N) * t(K_vec) ) /N
  H = ( sig_hat_ld_eps^2 * (K_vec * t(K_vec)  ) +  (sig_hat_ld_eps^2 * Gm )  + (sig_hat_eps * U_hat ) ) / N

  Dmat = 2*H
  sc = norm(Dmat, "2")
  dvec = -c
  # Rescale the objective function to reinv the overflow issue
  Dmat = Dmat / sc
  dvec = dvec / sc
  Amat = vcat(one_K, diag(one_K), -diag(one_K))
  Amat = Amat'
  bvec = c(1, zero_K, -one_K)

  r_qp = inv_QP(Dmat, dvec, Amat, bvec, meq=1)

  obj_val=  r_qp$value
  opt_w = round(r_qp$solution,4)


  return(list(opt_w=opt_w, obj_val=obj_val))
end





#obj_fnt=function(x, c, H) {c * x + 0_5* t(x) * H * x }


get_weights_knitro = function(y, x_end, x_exo, z, H_tilde, sig_hat_ld, P_all, Gm){
  library("KnitroR")

  # Declare constants
  z = as_matrix(z)
  K = ncol(z)
  N = nrow(z)

  x = cbind(x_end, x_exo)
  opt_w = rep(NA, K)

  # To get the preliminary estimate bt_tilde
  r1 = first_stage_GOF(y=y, x=x, z=z, H_tilde=H_tilde, sig_hat_ld=sig_hat_ld, criteria='Mallows')
  r1 = first_stage_GOF(y=y, x=x, z=z, H_tilde=1, sig_hat_ld=sig_hat_ld*H_tilde^2, criteria='Mallows')
  K_tilde = r1$opt_K
  # R_hat_all = r1$obj_value    # (K x 1) vector for R_hat (K), required

  # Update H_tilde using K_tilde
  H_tilde_update = ( t(x) * P_all[ , , K_tilde] * x ) / N

  # Calculate sig_hat_eps, sig_hat_ld, sig_hat_ld_eps, H_hat using 'get_sigma_DN'
  r2 = get_sigma_DN(y=y ,x=x ,z=z, opt_K=K_tilde, H_tilde=H_tilde_update)
  sig_hat_eps = as_numeric(r2$eps2)
  sig_hat_ld = as_numeric(r2$ld2)
  sig_hat_ld_eps = as_numeric(r2$ld_eps2)


  # Calculate u_hat_ld and U_hat
  u_stack = matrix(0, N, K)
  for (i in (1:K)){
    u_stack[ , i]=( P_all[ , , K] - P_all[ , , i] ) * x * inv(H_tilde_update)
  }
  U_hat = t(u_stack) * u_stack              # (K x K) matrix for S(K)



  # Generate a (K x 1) vector K_vec s_t_ K_vec[i]=i
  K_vec = c(1:K)




  # QP problem using Knitro

  # These parameters and functions are dropping (H^{-1})^2 as in the Ox codes of KO_
  # They should be equivalent to the one that uses the original formular_
  # If you want to use these formular then change the objective function in the knitro call from obj_fnt to obj_fnt_KO

    # c_KO = as_numeric(2 * sig_hat_eps * sig_hat_u * t(K_vec)  )
  # H_KO = 2 * ( sig_hat_u_eps^2 * (K_vec * t(K_vec)  ) + (sig_hat_eps * U_hat_KO  ) - (sig_hat_eps * sig_hat_u * Gm ) )
  # obj_fnt_KO=function(x) {t(c_KO) * x + 0_5* t(x) * H_KO * x }
  #---------------------------------------------

  c_fnt = function(x) {
    K_ones = rep(1,K)
    return(K_ones * x)
  }
  #  This expression comes from Equation (2_6) in KO_
  #  c = as_numeric(2 * sig_hat_eps * sig_hat_ld * t(K_vec)  ) /N
  #  H = 2 * ( sig_hat_ld_eps^2 * (K_vec * t(K_vec)  ) + (sig_hat_eps * U_hat  ) - (sig_hat_eps * sig_hat_ld * Gm )  ) / N

  # Equation (2_5) is rewritten: H is a matrix for the quadratic term; c is a matrix for the linear term
  c = as_numeric(   -2*t(K_vec)*(sig_hat_eps*sig_hat_ld + 4*sig_hat_ld_eps^2)  + 2 * sig_hat_eps * sig_hat_ld * t(K_vec)  ) /N
  H = 2 * ( sig_hat_ld_eps^2 * (K_vec * t(K_vec)  ) + (sig_hat_eps * U_hat  ) + (sig_hat_ld_eps^2 * Gm )  ) / N


  obj_fnt=function(x) {t(c) * x + 0_5* t(x) * H * x }


  eval_grad_f = function(x) {
    grad_f = H * x + c
    return( grad_f )
  }

  # Jacobian callback
  eval_jac_g = function( x ) {
    return( rep(1,K) )
  }

  cL = c(1_0)
  cU = rep(1_0)

  xL = rep(0,K)
  xU = rep(1,K)

  #
  # Knitro algorithms: 1 = Interior/Direct
  #                    2 = Interior/CG
  #                    3 = Active Set
  #                    4 = SQR
  # ms_enable: multiple starting points, T or F
  # URL: https://www_artelys_com/tools/knitro_doc/2_userGuide/algorithms_html
  #

  #r7=knitro(nvar=K, ncon=K+1, nnzJ=-1, x0=w0, objective = obj_fnt, constraints=c_fnt ,cL=cL, cU=cU, options=knitro_opt)
  r7=knitro(nvar=K, ncon=1, gradient=eval_grad_f,  jacobian=eval_jac_g, objective = obj_fnt, constraints=c_fnt ,cL=cL,  cU=cU, xL=xL, xU=xU)
  sol7=r7$x
  obj_val=obj_fnt(sol7)
  opt_w = round(sol7,4)


  return(list(opt_w=opt_w, obj_val=obj_val))
}





function tsls_KO(y, x_end, x_exo, z_excl, ld=NULL, method_pre, z_pre) {
  library('quadprog')

  # Declare constants
  z = hcat(x_exo,z_excl)
  x = hcat(x_end, x_exo)
  iv_pre = hcat(x_exo, z_pre)
  N = size(z, 1)
  K = size(z, 2)
  K_excl = size(z_excl, 2)
  d_2 = size(x_exo, 2)  # dim of included exogenous regressros
  if (is_null(d_2)) d_2=0
  d_x = size(x, 2)

  # Generate a (K x 1) vector K_vec s_t_ K_vec[i]=i
  K_vec = c(1:K_excl)

  # Generate a (K x 1) vector of 1's
  one_K = rep(1,K_excl)

  # Generate a (K x 1) vector of 0's
  zero_K = rep(0,K_excl)

  # Generate Gamma matrix such that GM[i,j] = min{i,j}
  Gm = matrix(0, K_excl, K_excl)
  for  i = 1:K_excl
    for  j = 1:K_excl
      Gm[i, j] = minimum(c(i,j))
  end
end

  # Construct projection matrices of each model and its weighted averages_
  P_all = array( NA, c(N, N, K_excl) )
  for  i = 1 : K_excl
    P_m = proj_m(z[,(1:(d_2+i))])
    P_all[ , , i] = P_m
end

  # If lambda (ld) is not defined and use the equal weight
  if (is_null(ld))
    ld = rep(1/d_x, d_x)
 end

  # Preliminary Regression to get sig2's
  prelim = pre_est(method_pre, y, x_end, x_exo, z_pre, z, P_all, ld, d_2)
  sig2_eps = prelim$sig2_eps
  sig2_ld = prelim$sig2_ld
  sig_ld_eps = prelim$sig_ld_eps
  sig_u_eps = prelim$sig_u_eps
  H_inv = prelim$H_inv
  u_pre = prelim$u_pre
  f_pre = prelim$f_pre


  # Evaluate the obj_ fun
  # Calculate U_hat
  u_stack = matrix(0, N, K_excl)
  for (i = 1:K_excl
    u_stack[ , i]=( P_all[ , , K_excl] - P_all[ , , i] ) * x * H_inv * ld
  end
  U_hat = u_stack'* u_stack

  # Calculate Sig_u
  Sig_u = (u_pre'* u_pre) / N

  # Calculate B_ld_N
  B_N_1 = sig2_eps * Sig_u
  B_N_2 = d_x * (sig_u_eps * t(sig_u_eps))
  B_N_3 = matrix(0, d_x, d_x)
  for (i_B = (1:N)
    f_i = f_pre[i_B, : ]
    term_1 = f_i * sig_u_eps'* H_inv * sig_u_eps * f_i'
    term_2 = f_i *sig_u_eps' * H_inv * f_i * sig_u_eps'
    term_3 = sig_u_eps * f_i' * H_inv * sig_u_eps * f_i'
    B_N_3 = B_N_3 + term_1 + term_2 + term_3
  }
  B_N_3 = B_N_3 / N
  B_N = 2*(B_N_1 + B_N_2 + B_N_3)
  B_ld_N = drop( ld' * H_inv * B_N * H_inv * ld )

  # Equation (2_5) is rewritten: H is a matrix for the quadratic term; c is a matrix for the linear term
  c = as_numeric(-(B_ld_N) * t(K_vec) ) /N
  H = ( sig_ld_eps^2 * (K_vec * K_vec'  ) +  (sig_ld_eps^2 * Gm )  + (sig2_eps * U_hat ) ) / N

  Dmat = 2*H
  sc = norm(Dmat, "2")
  dvec = -c
  # Rescale the objective function to reinv the overflow issue
  Dmat = Dmat / sc
  dvec = dvec / sc
  Amat = rbind(one_K, diag(one_K), -diag(one_K))
  Amat = t(Amat)
  bvec = c(1, zero_K, -one_K)

  r_qp=inv_QP(Dmat, dvec, Amat, bvec, meq=1)

  obj_val=r_qp$value
  opt_w = round(r_qp$solution,4)

  cat('Optimal Weights = ', opt_w,'\n')

  # Calculate the averaged projection matrix with the optimal weight_
  P_W = matrix(0, N, N)
  for  i = 1 : K_excl
    P_W = P_W + opt_w[i] * P_all[ , ,i]
  end

  # Calculate the 2SLS estimator with P_W and return the result
  bt_hat = inv(x'*P_W*x)*(x'*P_W*y)
  cat('bt_hat (KO) = ', bt_hat, '\n')

  return(list(bt_hat = bt_hat, opt_w = opt_w, obj_val = obj_val, P_KO = P_W))
}


cluster_se = function(id, bt, x, y, P){
  eps_hat = y - x * bt
  id_index = list()
  P_g = list()
  eps_g_hat = list()

  G = length(table(id))
  id_set = names(table(id))
  d = size(x, 1)
  mid_sum_G = matrix(0,d,d)

  for i in id_set
    id_index[[i]] = which(id==i)
    P_g[[i]] = P[id_index[[i]], , drop=FALSE]
    eps_g_hat[[i]] = eps_hat[id_index[[i]]]
    mid_sum_G = mid_sum_G + t(x) * t(P_g[[i]]) * eps_g_hat[[i]] * t(eps_g_hat[[i]]) * P_g[[i]] * x
  end

  outer_part = inv(x' * P * x)
  Sig_hat = N*outer_part * (mid_sum_G) * outer_part
  vcov = Sig_hat/N
  se = sqrt(diag(vcov))
  cat('hetero-cluster robust se =', se,'\n')
  return(list(se=se, vcov=vcov))
end


# Count the number of products that have inelastic demand
function inelastic(ap, price, share)
  sum(ap*(price)*(1-share) > -1)
end

function coverage(ap, se, ap_0)
  upper = ap + 1_96*se
  lower = ap - 1_96*se
  cover = (upper>ap_0) & (lower<ap_0)
  #re=cbind(ap,se,upper,lower,cover)
  #print(head(re))
  mean(cover)
end
