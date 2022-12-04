# S&P 500 data, the Covariance matrix is not positive define
#https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda

using CodecXz
using RData
using EfficientFrontier
using LinearAlgebra
using Serialization

# download https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda to /tmp
sp5h = load("/tmp/muS.sp500.rda")

E = values(sp5h["muS.sp500"]["mu"])
V = values(sp5h["muS.sp500"]["covar"])

display(norm(V - V'))

V = (V+V')/2    #make sure symetry
maximum(V - V') # 0.0
display(norm(V - V'))
display(rank(V))    #263


#=  #v0.1.0
EfficientFrontier.noShortsale(E, V)
t0 = time()
aCL = EfficientFrontier.ECL()
t1 = time()
display(aCL)
display(t1 - t0)    #0.25 seconds
Z = EfficientFrontier.CornerP()
=#

#v0.2.0
t0 = time()
P = Problem(E, V; equilibrate=false)
ts = @elapsed aCL = EfficientFrontier.ECL(P)
aEF = EfficientFrontier.eFrontier(aCL, P)
t1 = time()
display(ts)   #0.20 seconds
display(t1 - t0)    #0.25 seconds
display(aEF.Z)