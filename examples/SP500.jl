# S&P 500 data, the Covariance matrix is not positive define
#https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda

using CodecXz
using RData
using EfficientFrontier
using LinearAlgebra
using Serialization

println("\n--- Status-Segment Method vs Markowitz's CLA  ---\n")

# download https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda to /tmp
sp5h = load("/tmp/muS.sp500.rda")

E = values(sp5h["muS.sp500"]["mu"])
V = values(sp5h["muS.sp500"]["covar"])

#display(norm(V - V'))

V = (V+V')/2    #make sure symetry
#maximum(V - V') # 0.0
#display(norm(V - V'))
println("rank(V):  ", rank(V))  #263


if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    include("./Markowitz.jl")
    using .Markowitz
end

m = markowitz(E, V)
unit_sum(m) # total weight = 100%
ts = @elapsed f = frontier(m)
println("Markowitz CLA:  ", ts, "  seconds")     #0.14 seconds
display(f.weights)

t0 = time()
P = Problem(E, V; equilibrate=false)
ts = @elapsed aCL = EfficientFrontier.ECL(P)
aEF = eFrontier(aCL, P)
t1 = time()
println("Status-Segment Method:  ", ts, "  seconds")    #0.039446974  seconds
#display(t1 - t0)    #0.083 seconds
display(aEF.Z)
