# This script lists some auxiliary function used in the simulation study --------------------

#=
 Pkg.add("NLsolve")
 Pkg.add("Statistics")
 Pkg.add("Random")
 Pkg.add("Plots")
 Pkg.add("KernelEstimator")
 Pkg.add("SpecialFunctions")
 Pkg.add("CSV")
 Pkg.add("DataFrames")
 Pkg.add("LinearAlgebra")
 Pkg.add("Distributions")
 Pkg.add("Optim")
 Pkg.add("HypothesisTests")

=#

 using Random      # for sampling randomly in a vector
 using Statistics
 # using NLsolve
 using DataFrames
 using Distributions
 # using BlackBoxOptim
 using LinearAlgebra
 #using HypothesisTests
 # using Optim
 using CSV
 # using IterTools
 # using Combinatorics
 using SpecialFunctions
 using KernelEstimator
 # using Plots
 # using OptimPackNextGen
 using Optim
 using Plots

# Gram Schmidt orthogonalization
function GS(mat; norm = "yes")
newmat = zeros(size(mat))
newmat[:, 1] = mat[:, 1]
proj = zeros(size(mat, 1), 1)
   for i = 2:size(mat, 2)
       v = mat[:, i]
       u = newmat[:, i - 1]
       proj = proj + ((v'*u)/(u'*u)) * u
       newmat[:, i] = mat[:, i] - proj
       if norm == "yes"
       newmat[:, i] = newmat[:, i] /sqrt(dot(newmat[:, i], newmat[:, i]) )
       end
   end
   if norm == "yes"
 newmat[:, 1] = newmat[:, 1] /sqrt(dot(newmat[:, 1], newmat[:, 1]) )
   end
    return newmat
end

# Gaussian kernel function
function gaussian(u)
  ker = (1/sqrt(2*pi)).*exp.( -u.^2 ./2)
  return (ker)
end

# Gamma function (for student density)
#=
function gamma_f(n)
    gam = factorial(n - 1)
    return(gam)
end
=#

# Normal density
function dnorm(u)
  dens = (1/sqrt(2*pi)).*exp.( -u.^2 ./2)
  return (dens)
end

# Student density
function stud(x, nu)
    nu = 3
    dens = ( ( gamma((nu + 1)/2) )/(sqrt(pi*nu)*gamma(nu/2)) ).*(1 .+ x.^2 ./nu).^(-(nu + 1)/2)
    return (dens)
end

# A Data Generating Process function, with linear, quadratic and cubic designs. The DGP of Antoine & Lavergne (2019) is also an option under design = "A&L".
function dgp(n, beta_0, mu, sigma; design = "linear")
    RNG = MvNormal( mu, sigma)
    RNG = rand(RNG, n)'
    u = reshape( RNG[:, 1], n, 1 )
    v = reshape( RNG[:, 2], n, 1 )
    z = reshape( randn(n), n, 1 )
    if design == "linear"
    x =  z + v
    elseif design == "log"
	x = log.(abs.(z)) + v
    elseif design == "A&L"
    x = (1/n^(1/4))*(3 .*z-z.^3) + (1/sqrt(n))*(z.^2 - ones(n, 1)) + v
elseif design == "cube"
	x = z - z.^2/3 + z.^3/6 + v
elseif design == "square"
    x = z.^2 + v
elseif design == "normal density"
	x = (1/sqrt(2*pi))*exp.(-0.5*(z.^2)) + v
elseif design == "binary_x"
    x = z + v .- 0.5
	x .> 0
    end
    y = x*beta_0 + u
    data = DataFrame( hcat(y, x, z, u, v), :auto )
    rename!(data, [:y, :x, :z, :u, :v])
    return data
end

# Function that gives names according to some character and the numeric variable k
function name_k(k, name)
   name = string(name, "_", k)
   return name
end

# Example
#=
 n = 100
 mu = [0.0, 0]
 sigma = [1.0, 0.5, 0.5, 1.0]
 sigma = reshape(sigma, 2, 2)
 data = dgp(n, mu, sigma)
 names(data)
# Select particular variables in the dataframe
 data.x
 data.y
=#
