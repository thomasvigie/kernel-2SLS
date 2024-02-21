# This script lists classical estimators used to estimate linear models, and estimators using kernels ---------------------------

# OLS function. Include a vector of ones in x if one wants to include an intercept
function ols(x, y; intercept = "false")
    n = size(y, 1)
    x = reshape(x, size(x, 1), size(x, 2))   # If x is a vector, it turns it into a matrix so inv can be used
    if intercept == "true"
      x = hcat(ones(n, 1), x)
    else x = x
    end
     beta =inv(x'*x)* x'*y
	 res = y - x*beta
R2 = 1 - sum(res.^2 )/sum( (y .- mean(y)).^2 )

    K = size(x, 2)
	 sig2hat = (sum(res.^2 ))/(n - K)
     var_homo = diag(sig2hat*inv( x'*x ) )
#     var_homo = reshape(var_homo, length(var_homo), 1)
     std_errors = sqrt.(var_homo)

#     var_hetero = inv( x'*x )*(x'*rep.col(res.^2, ncol(x)) *x)*inv( x'*x )
t_stats = beta ./ std_errors
p_values = (1 .- cdf.(Normal(0, 1), abs.(t_stats)) ).*2
result = DataFrame( hcat(beta, std_errors, t_stats, p_values), :auto )
rename!(result, [:"Estimates", :"Standard errors", :"t-statistics", :"p-values"])
#insertcols!(result, 1, :variable => ["β0","β1","β2","β3"])	# Naming the rows
return result#, R2#, tstats = t_stats, p_values = p_values) #, Decision = decisions , var_hetero = diag(var_hetero))
end

# Example
#=
using Statistics
using LinearAlgebra
using DataFrames

n = 1000
u = reshape( randn(n), n, 1 )
x1 = reshape( randn(n), n, 1 )
x2 = reshape( randn(n), n, 1 )
y = 1 .+ 0.1 .*x1 .- 3 .*x2 .+ u
x = hcat(x1, x2)
ols(x, y; intercept = "true")
=#

# The Two Stage Least Squares (TSLS) function. Include a vector of ones in x and/or z if one wants to include an intercept the first and/or second stage
function tsls(x, y, z; intercept = "true")
    x = reshape(x, size(x, 1), size(x, 2))   # If x is a vector, it turns it into a matrix so inv can be used
    z = reshape(z, size(z, 1), size(z, 2))   # If z is a vector, it turns it into a matrix so inv can be used
    n = size(y, 1)
    z = hcat(ones(size(z, 1), 1), z)
    pz = z*inv(z'*z)*z'

    if intercept == "true"
      x = hcat(ones(n, 1), x)
    else x = x
    end

    beta = inv(x'*pz*x)* x'*pz*y
	res = y - x*beta
    K = size(x, 2)

    sig2hat = (sum(res.^2 ))/(n - K)
    R2 = 1 - sum(res.^2 )/sum( (y .- mean(y)).^2 )

    Q = z'*z
    Gam = x'*z
	var_homo = diag(sig2hat*inv(( Gam*inv(Q)*Gam' )))
    std_errors = sqrt.(var_homo)
    t_stats = beta ./ std_errors
    p_values = (1 .- cdf.(Normal(0, 1), abs.(t_stats)) ).*2
    result = DataFrame( hcat(beta, std_errors, t_stats, p_values), :auto )
    rename!(result, [:"Estimates", :"Standard errors", :"t-statistics", :"p-values"])
    #insertcols!(result, 1, :variable => ["β0","β1","β2","β3"])	# Naming the rows
    return result#, R2#, tstats = t_stats, p_values = p_values) #, Decision = decisions

    #Omega = (z*rep.col(res.^2, ncol(xhat)) )*z'    # For the heteroskedastic variance
#	var_hetero = inv(( Gamma*inv(Q)*Gamma' ))*Gamma*inv(Q)*Omega*inv(Q)*Gamma*inv(( Gamma*inv(Q)*Gamma' ))
end

# The k-class estimator of Theil(1958). It includes the OLS (k = 0) and TSLS (k = 1) as special cases
function kclass(x, y, z; k = 1)
    n = size(y, 1)
    Pz = z*inv(z'*z)*z'
    Mz = Matrix{Float64}(I, n, n) - Pz
    betahat = inv(x'*(Matrix{Float64}(I, n, n) - k*Mz)*x)*x'*(Matrix{Float64}(I, n, n) - k*Mz)*y
    return betahat
end

# The Limited Information Maximum Likelihood (LIML) estimator, also a special case of the k-class estimator.
function liml(x, y, z)
    n = size(y, 1)
    Pz = z*inv(z'*z)*z'
    Mz = Matrix{Float64}(I, n, n) - Pz
    w = hcat(y, x)
    function obj(k)
        expr = det(w'*w - k.*w'*Mz*w)
    end
    kliml = nlsolve(obj, vec([0.5])).zero
    betahat = inv(x'*(Matrix{Float64}(I, n, n) - kliml.*Mz)*x)*x'*(Matrix{Float64}(I, n, n) - kliml.*Mz)*y
    return betahat
end

# The WMDF estimator of Antoine & Lavergne (2014)
function wmdf(x, y, z; intercept = "false")
    n = size(y, 1)
     if intercept == "true"
       Y1star = hcat(ones(n, 1), x)
     else Y1star = x
     end
      Ystar = hcat( y , Y1star )
      z = z./repeat(std(z, dims = 1), n, 1)  # some rescaling

    if size(z, 2)>1
      Ktilde = ones(n, n)
      for i = 1:size(z, 2)
        K = dnorm( repeat(z[:, i], 1, n) -  repeat(z[:, i]', n, 1) )
        K = K - Diagonal(K)
        Ktilde = Ktilde.*K
      end
     else
        Ktilde = dnorm( repeat(z, 1, n) -  repeat(z', n, 1) )
        Ktilde = Ktilde - Diagonal(Ktilde)
    end

      eig = eigvals(inv(Ystar'*Ystar)*(Ystar'*Ktilde*Ystar))
      lambdatilde = minimum(eig)

      lambda_WMDF = (lambdatilde - (1 - lambdatilde)/n)/(1 - (1 - lambdatilde)/n)
      In = Diagonal(vec(ones(n, 1)))
      WMD = inv(Y1star'*(Ktilde-lambdatilde.*In)*Y1star)*(Y1star'*(Ktilde-lambdatilde*In)*y)
      WMDF = inv(Y1star'*(Ktilde-lambda_WMDF.*In)*Y1star)*(Y1star'*(Ktilde-lambda_WMDF*In)*y)
      return (WMD = WMD, WMDF = WMDF)
end

# Example
#=
 include("auxiliary functions.jl")
 n = 1000
 mu = [0.0, 0]
 sigma = [1.0, 0.5, 0.5, 1.0]
 sigma = reshape(sigma, 2, 2)
 # Generate some data where the coeffecieint of interest, beta_0, is equal to 1.
 data = dgp(n, 1, mu, sigma; design = "linear")

 OLS = ols(data.x, data.y)
 TSLS = tsls(data.x, data.y, data.z)
 LIML = liml(data.x, data.y, data.z)
 KCLASS_0 = kclass(data.x, data.y, data.z, k = 0)    # returns the same as OLS
 KCLASS_1 = kclass(data.x, data.y, data.z, k = 1)    # returns the same as TSLS
 WMDF = wmdf(data.x, data.y, data.z, intercept = "false")
 NP = np_tsls(data.x, data.y, data.z, option = localconstant)
=#


# Estimators computing or using kernels -------------------------------------

# This function computes the Nadaraya-Watson (type = "lc") and local linear (type = "ll") kernel regression estimators
# A Gaussian kernel is used
# No ridging to handle small bandwidths, instead a numerical adjustment to avoid dividing by 0. It is not ideal.
function lcll( dep, expl; bw = 0.1, type = "ll" )
    if type != "lc" && type != "ll"    # Check how to make that statement work
		println("Mate, the type of kernel regression should be either 'lc' (local constant) or 'll' (local linear) ")
	end
  n = size(expl, 1)
  firsterm = repeat(expl, 1, n)
  secterm = firsterm'
  U = (firsterm-secterm)./bw
  K = gaussian(U)   # Replace that by any other kernel function if desired
  if type == "ll"

    Sn1K = sum(K.*(firsterm-secterm), dims = 2)
    Sn1K = repeat(Sn1K, 1, n)

    Sn2K = sum(K.*(firsterm-secterm).^2, dims = 2)
    Sn2K = repeat(Sn2K, 1, n)

# Using the notation from Cheng (1997) that implements ridging towards the local cosntant estimator
# Racine and Hayfiled have epsilon = 1/n
eta = 1/n
 w_K = K.*(Sn2K - (firsterm-secterm).*Sn1K)

top = w_K
bottom = sum(w_K, dims = 2)
Ksum = sum(K, dims = 2)

#      println("Careful mate ! I applied some ridging")

    ind = findall( x -> x == 0, bottom)
     ind_1 = [i[1] for i in ind]   # extract all the row numbers with zeros
#if ind_1 != 0  # returns a message even when no ridging seems to have been implemented
#	println(string("Careful mate ! I applied some ridging on rows ", ind_1))
#end

# Ridging around the Nadaraya Watson estimator like Racine and Hayfiled in R
# Ridging using eta = 1/n
      bottom[ind]= bottom[ind] + eta.*Ksum[ind]
      top[ind_1,:] = top[ind_1, :] + eta.*K[ind_1, :]


L =  top ./repeat(bottom, 1, n)

else
	top = K
	bottom = sum(K, dims = 2)
	L =  top ./repeat(bottom, 1, n)
end
  fit = L*dep
  return(fit)
end

# Example
#=
using Random      # for sampling randomly in a vector
using Statistics
using DataFrames
using Distributions
using CSV
using SpecialFunctions

beta_0 = 1
n = 100
mu = [0, 0]
sigma = [1.0, 0.5, 0.5, 1.0]
sigma = reshape(sigma, 2, 2)
RNG = MvNormal( mu, sigma )
RNG = rand(RNG, n)'
u = reshape( RNG[:, 1], n, 1 )
v = reshape( RNG[:, 2], n, 1 )
z = reshape( randn(n), n, 1 )
x =  z + v
y = x*beta_0 + u
data = hcat(y, x, z, u, v)
data = DataFrame(data, [:y, :x, :z, :u, :v])
#CSV.write("F:\\Kernel 2SLS\\Julia simulations\\try data.csv", data)
=#


# Extract the "hat" matrix when using the Nadaraya-Watson or local linear kernel estimators.
function kernel_hat_matrix( dep, expl; bw = 0.1, type = "lc" )
      n = size(expl, 1)
      firsterm = repeat(expl, 1, n)
      secterm = firsterm'
      U = (firsterm-secterm)./bw
if type == "ll"

      K = gaussian(U)
      # if type == "ll"

      Sn1K = sum(K.*(firsterm-secterm), dims = 2)
      Sn1K = repeat(Sn1K, 1, n)

      Sn2K = sum(K.*(firsterm-secterm).^2, dims = 2)
      Sn2K = repeat(Sn2K, 1, n)

  # Using the notation from Cheng (1997) that implements ridging using a student density
  # S = stud(U, 3)   # Replace that by any other kernel function if desired
  # Sn1L = sum(S.*(firsterm-secterm), dims = 2)
  # Sn1L = repeat(Sn1L, 1, n)
  #
  # Sn2L = sum(S.*(firsterm-secterm).^2, dims = 2)
  # Sn2L = repeat(Sn2L, 1, n)

  # eta = sqrt(n)
  # Racine and Hayfiled have epsilon = 1/n
  eta = 1/n
  #eta = 0

   w_K = K.*(Sn2K - (firsterm-secterm).*Sn1K)
  # w_L = S.*(Sn2L - (firsterm-secterm).*Sn1L)

  top = w_K
  bottom = sum(w_K, dims = 2)
  Ksum = sum(K, dims = 2)

  #top_ridge = w_K + eta.*w_L
  #bottom_ridge = sum(w_K + eta.*w_L, dims = 2)

  #      println("Careful mate ! I applied some ridging")

      ind = findall( x -> x == 0, bottom)
       ind_1 = [i[1] for i in ind]   # extract all the row numbers with zeros
  # Ridging as in Cheng  et al (1997) using a Student density
  #      bottom[ind]= bottom_ridge[ind]
  #      top[ind_1,:] = top_ridge[ind_1, :]

  # Ridging around the Nadaraya Watson estimator like Racine and Hayfiled in R
        bottom[ind]= bottom[ind] + eta.*Ksum[ind]
        top[ind_1,:] = top[ind_1, :] + eta.*K[ind_1, :]

  L =  top ./repeat(bottom, 1, n)
else
  top = K
  bottom = sum(K, dims = 2)
  L =  top ./repeat(bottom, 1, n)
end
  return(L)
end

# Example using the try data set
#L = kernel_hat_matrix( data.x, data.z; bw = 0.1, type = "lc" )
#diag(L)
#1 .- diag(L)


# Computes the leave-one-out cross validation criterion
function loocv( dep, expl, bw, type = "ll" )
    n = size(expl, 1)
    L = kernel_hat_matrix(dep, expl; bw = bw, type = type)
#    ind = findall( x -> x == 1, diag(L))   # Sort this out
#    for i in ind    # if a diagonal element is equal to 1, it will screw up cross validation...
#    L[i, i] = 1 - 0.001
#    end

    # LOO = zeros(n, n)
    # for i in 1:n
    #  LOO[i, :] = L[i, :]/sum(L[i, 1:end .!= i]) # Fix that issue here
    #  LOO[i, i] = 0
    # end
 # LOO[LinearAlgebra.diag(LOO)] = 0
  fit = L*dep
  # objective = sum((fit - dep).^2)
  objective = (1/n)*sum( ( (dep - fit)./(1 .- diag(L) ) ).^2)
  return(objective)
end

# Example witht the try data set
#loocv(data.x, data.z, 0.001, "lc" )
#loocv(data.x, data.z, 1, "lc" )
#loocv(data.x, data.z, 9539536, "lc" )


function opti_bandwidth(dep, expl, type = "ll" )
# Might be used as starting values for the minimization
hlb = 0.1
hub = 1000000
opti_h = Optim.minimizer(Optim.optimize(bw->loocv(dep, expl, bw, type), hlb, hub))
# opti_h = optimize(loocv, [1])    # starting value at 1
return(opti_h)
end

# Example using the try data set
#op_h = opti_bandwidth(data.x, data.z, "lc" )
#loocv(data.x, data.z, op_h, "lc" )
#red = lcll( data.x, data.z; bw = op_h, type = "lc" )  # don't mind the "lc", right now only the local linear is coded
#red

# Implements the AIC criterion function for the bandwidth from Hurvic, Simonoff & Tsai (1999)
 function hurvic_AIC(dep, expl; bw = 0.1, type = "ll")
  H = kernel_hat_matrix(dep, expl, bw = bw, type = type)
  n = size(dep, 1)
  sig2 = ( dep'*( Matrix{Float64}(I, n, n) - H )'*( Matrix{Float64}(I, n, n) - H )'*dep ) / n
  crit = log(sig2) + (1 +  tr(H)/ n)/( 1 - (tr(H) + 2) / n)
  return crit
end

# Function that computes the optimal bandwidth by minimizing to the iv_AIC function. Still in the works, use at your own risk.
function min_AIC(dep, expl; lower = 0.01, upper = 50.0, type = "ll") #, loo = FALSE )
   # bboptimize(iv_AIC; SearchRange = [lower, upper], NumDimensions = 1, expl = expl, endo = endo, instru = instru, type = type)   # Pkg.add("BlackBoxOptim")  DE type algorithm
   n = size(expl, 1)
   # H = kernel_hat_matrix(endo, instru, bw = h, type = type) #, loo)
   # # In the IV case, the effective hat matrix is (X'H'X)^(-1) X'H'
   # H = endo*inv(endo'*H'*endo)*endo'*H'
   # sig2 = ( expl'*( Matrix{Float64}(I, n, n) - H )'*( Matrix{Float64}(I, n, n) - H )*expl ) / n
   # # crit = log(sig2) + (1 +  tr(H)/ n)/( 1 - (tr(H) + 2) / n)

   function obj(h)  # I define the criterion function directly inside the min_AIC_iv function, so that it depends on the bandwidth only. Makes the minimization code easier
       # n = size(expl, 1)
       H = kernel_hat_matrix(dep, expl, bw = h, type = type)
       sig2 = ( dep'*( Matrix{Float64}(I, n, n) - H )'*( Matrix{Float64}(I, n, n) - H )'*dep ) / n
       crit = log(sig2) + (1 +  tr(H)/ n)/( 1 - (tr(H) + 2) / n)
       return crit
   end
   # With a classical optimizer
   sol = optimize(obj, lower, upper, Brent() )
   minimum = Optim.minimum(sol)
   minimizer = Optim.minimizer(sol)
   # With a Differential Evolution algorithm
   # sol = bboptimize(obj; SearchRange = (0.01, 50.0), NumDimensions = 1  )
   # minimum = best_fitness(sol)
   # minimizer = best_candidate(sol)
  # convergence = sol$convergence
  return(opt_obj = minimum, opt_h = minimizer)
end

# An extension of the AIC function of Hurvic, Simonoff & Tsai (1999) for the IV case
function iv_AIC(expl, endo, instru; bw = 0.1, type = "ll")
  n = size(expl, 1)
  H = kernel_hat_matrix(endo, instru, bw = bw, type = type)
  # In the IV case, the effective hat matrix is (X'H'X)^(-1) X'H'
  H = endo*inv(endo'*H'*endo)*endo'*H'
  sig2 = ( expl'*( Matrix{Float64}(I, n, n) - H )'*( Matrix{Float64}(I, n, n) - H )*expl ) / n
  crit = log(sig2) + (1 +  tr(H)/ n)/( 1 - (tr(H) + 2) / n)
  return(crit)
end

# Example:
#=
include("auxiliary function.jl")
 n = 1000
 mu = [0.0, 0]
 sigma = [1.0, 0.5, 0.5, 1.0]
 sigma = reshape(sigma, 2, 2)
 data = dgp(n, 1, mu, sigma)
 AIC_h = hurvic_AIC(data.y, data.x, bw = 0.1)
 AIC_iv_h = iv_AIC(data.y, data.x, data.z; bw = 0.1, type = "ll")
=#

# Function that computes the optimal bandwidth by minimizing to the iv_AIC function. Still in the works, use at your own risk.
function min_AIC_iv( expl, endo, instru; lower = 0.01, upper = 50.0, type = "ll") #, loo = FALSE )
   # bboptimize(iv_AIC; SearchRange = [lower, upper], NumDimensions = 1, expl = expl, endo = endo, instru = instru, type = type)   # Pkg.add("BlackBoxOptim")  DE type algorithm
   n = size(expl, 1)
   # H = kernel_hat_matrix(endo, instru, bw = h, type = type) #, loo)
   # # In the IV case, the effective hat matrix is (X'H'X)^(-1) X'H'
   # H = endo*inv(endo'*H'*endo)*endo'*H'
   # sig2 = ( expl'*( Matrix{Float64}(I, n, n) - H )'*( Matrix{Float64}(I, n, n) - H )*expl ) / n
   # # crit = log(sig2) + (1 +  tr(H)/ n)/( 1 - (tr(H) + 2) / n)

   function obj(h)  # I define the criterion function directly inside the min_AIC_iv function, so that it depends on the bandwidth only. Makes the minimization code easier
       # n = size(expl, 1)
        H = kernel_hat_matrix(endo, instru, bw = h, type = type) #, loo)
       # # In the IV case, the effective hat matrix is (X'H'X)^(-1) X'H'
        H = endo*inv(endo'*H'*endo)*endo'*H'
        sig2 = ( expl'*( Matrix{Float64}(I, n, n) - H )'*( Matrix{Float64}(I, n, n) - H )*expl ) / n
       crit = log(sig2) + (1 +  tr(H)/ n)/( 1 - (tr(H) + 2) / n)
       return crit
   end
   # With a classical optimizer
   sol = optimize(obj, lower, upper, Brent() )
   minimum = Optim.minimum(sol)
   minimizer = Optim.minimizer(sol)
   # With a Differential Evolution algorithm
   # sol = bboptimize(obj; SearchRange = (0.01, 50.0), NumDimensions = 1  )
   # minimum = best_fitness(sol)
   # minimizer = best_candidate(sol)
  # convergence = sol$convergence
  return(opt_obj = minimum, opt_h = minimizer)
end

# This function proposes a leave-one-out criterion for betahat
# Formula taken from Raffaele Saggio (econometrica)
function loo_MSE(x, y, z, h, regtype = "ll")
  xhat = lcll(x, z; bw = h, type = regtype)  # Ridging has been implemented
  # xhat <- np::npreg(bws = as.vector(h), txdat = z, tydat = x, ckertype = "gaussian", regtype = regtype)$mean
  HAT = x*inv(xhat'*x)*xhat'
  betahat =  myest(y, x, z, h; regtype = regtype).Estimates   # add the variance etc like in myest_np
  betas_loo = (x.*betahat .- y.*diag(HAT) ) ./ ( (1 .- diag(HAT)).*x )
  MSE = mean((betas_loo .- betahat ).^2)
  return MSE
end

function loo_opti(x, y, z; regtype = "ll")
  opti_h = optimize(h->loo_MSE(x, y, z, h, regtype), 0.001, 1000)
  # Check if the algorithm converged and throw an error otherwise
  # converged(opti_h) || error("Failed to converge in $(iterations(opti_h)) iterations")
  opt_h = opti_h.minimizer
  value = opti_h.minimum
  betahat =  myest(y, x, z, opt_h, regtype = regtype)
#  return DataFrame(betahat = betahat, optimal_bandwidth = opt_h, minimum_objective = value)
return betahat
end


# This function computes optimal bandwidh for the kernel-IV estimator
# I need to include the iv_AIC criterion (from R) as an option as well
# Change the optimizer
function myest_bw(x, y, z; regtype = "ll")
	opti_h = optimize(bw->loocv(x, z, bw, regtype), 0.001, 1000)
    # Check if the algorithm converged and throw an error otherwise
#    converged(opti_h) || error("Failed to converge in $(iterations(opti_h)) iterations")
    bandwidth = opti_h.minimizer
	value = opti_h.minimum
    return DataFrame(bandwidth = bandwidth, value = value)
  end

# This function computes the kernel-IV estimator (called myest in R) for a general bandwidth h
# add all the options (variance, possibility of an intercept and exogenous variables etc)
function myest(x, y, z, h; regtype = "ll")
     xhat = lcll( x, z; bw = h, type = regtype )
    # xhat = lcll(x, z, h, regtype, kertype = "gaussian" )  # Ridging has been implemented
     beta = inv(xhat'*x)* xhat'*y

	 res = y - x*beta
     R2 = 1 - sum(res.^2 )/sum( (y .- mean(y)).^2 )

	 sig2hat = mean( (y - x*beta) .^2 )
     var_homo = sig2hat*inv(( xhat'*xhat ))
#     var_hetero = inv(( xhat'*xhat ))* (xhat'*rep.col(res.^2, ncol(xhat))*xhat) *inv(( xhat'*xhat ))
# Include var_hetero later
std_errors = sqrt.(var_homo)

#     var_hetero = inv( x'*x )*(x'*rep.col(res.^2, ncol(x)) *x)*inv( x'*x )
t_stats = beta ./ std_errors
p_values = (1 .- cdf.(Normal(0, 1), abs.(t_stats)) ).*2
result = DataFrame( hcat(beta, std_errors, h, t_stats, p_values), :auto )
rename!(result, [:"Estimates", :"Standard errors", :"Bandwidth", :"t-statistics", :"p-values"])
#insertcols!(result, 1, :variable => ["β0","β1","β2","β3"])	# Naming the rows
return result#, R2#, tstats = t_stats, p_values = p_values) #, Decision = decisions , var_hetero = diag(var_hetero))
end

#=
function gna(x, y; option = "oui")
a = x + y
b = x*y
	 return DataFrame(gnu=a, gni=b)
end
u = gna(1, 2, option = "non")
=#

# Example
#=
include("auxiliary function.jl")
 n = 1000
 mu = [0.0, 0]
 sigma = [1.0, 0.5, 0.5, 1.0]
 sigma = reshape(sigma, 2, 2)
 data = dgp(n, 1, mu, sigma)

 OLS = ols(data.x, data.y)
 TSLS = tsls(data.x, data.y, data.z)
 LIML = liml(data.x, data.y, data.z)
 KCLASS_0 = kclass(data.x, data.y, data.z, k = 0)
 KCLASS_1 = kclass(data.x, data.y, data.z, k = 1)
 WMDF = wmdf(data.x, data.y, data.z, intercept = "true")
 NP = np_tsls(data.x, data.y, data.z, option = localconstant)
 TV = kernel_tv(data.x, data.y, data.z; option = "ll", bw = 1)
 # Note that if the bandwidth is very high, the estimator coincides with teh TSLS one
 TV_500 = kernel_tv(data.x, data.y, data.z; option = "ll", bw = 500)
=#


# Example:
# mu = c(0,0)
# VCOV = matrix(c( 1 , 0.5 ,
#                  0.5 , 1 ), nrow = 2, ncol = 2)
# errors = MASS::mvrnorm(100, mu, Sigma = VCOV )
# z = rnorm(100)
# x = z + errors[, 2]
# y = 2*x + errors[, 1]
#
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1)
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1, loo = TRUE)
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1, type = "ll")
# iv_AIC(expl = y, endo = x, instru = z, h = 0.1, type = "ll", loo = TRUE)
