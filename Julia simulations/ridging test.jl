# This script tests the kernel esitmators "faits maison"
using CSV
using DataFrames
using KernelEstimator

include("auxiliary functions.jl")
include("Estimators.jl")

data = CSV.read("F:\\Kernel 2SLS\\Julia simulations\\try data.csv", DataFrame)   # Laptop computer
#data = CSV.read("G:\\Kernel 2SLS\\Julia simulations\\try data.csv", DataFrame)   # Desktop computer
# OLS
ols(data.x, data.y; intercept = "true")
# 2SLS
tsls(data.x, data.y, data.z; intercept = "true")
# WMDF
wmdf(data.x, data.y, data.z; intercept = "false")
# Kernel regression
pred = lcll( data.x, data.z; bw = 0.001)
print(pred[15, :])
print(pred[21, :])
print(pred[25, :])
print(pred[31, :])

np_pred = npr(data.z, data.x, xeval=data.z, reg=locallinear, kernel=gaussiankernel, h = 0.001)

show(stdout, "text/plain", pred)
show(stdout, "text/plain", np_pred)

# Look at the leave-one-out CV criterion function
loocv(data.x, data.z, 0.001, "lc" )
loocv(data.x, data.z, 1, "lc" )
loocv(data.x, data.z, 9539536, "lc" )

# Look at the optimal bandwidth using leave-one-out CV
op_h = opti_bandwidth(data.x, data.z, "lc" )
loocv(data.x, data.z, op_h, "lc" )
pred = lcll( data.x, data.z; bw = op_h, type = "lc" )  # don't mind the "lc", right now only the local linear is coded
pred

# Look at the information criterion from Hurvic Simonoff and Tsai (1997)
iv_AIC(data.y, data.x, data.z; bw = 0.1, type = "lc")
min_AIC_iv(data.y, data.x, data.z; lower = 0.01, upper = 50.0, type = "lc") #, loo = FALSE )

# Look at the optimal bandwidh with a leave-one-out MSE criterion function
kernel_TV = loo_opti(data.x, data.y, data.z, regtype = "ll")
bw_TV = kernel_TV.Bandwidth[1]

# Estimate with the CV criterion function
myest_bw(data.x, data.y, data.z; regtype = "lc")


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
