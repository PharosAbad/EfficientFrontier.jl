# S&P 500 data, the Covariance matrix is not positive define
#https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda

using CodecXz
using RData
using EfficientFrontier
using LinearAlgebra
using Serialization

println("\n--- connecting Critical Line Segments vs Markowitz's CLA  ---\n")

# download https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda to /tmp
sp5h = load("/tmp/muS.sp500.rda")

E = values(sp5h["muS.sp500"]["mu"])
V = values(sp5h["muS.sp500"]["covar"])

#display(norm(V - V'))

V = (V+V')/2    #make sure symetry
#maximum(V - V') # 0.0
#display(norm(V - V'))
println("rank(V):  ", rank(V))  #263

#sort the data improve the speed in this example
ip = sortperm(vec(E), rev=true)
E = vec(E[ip])
V = V[ip, ip]


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
println("connecting Critical Line Segments:  ", ts, "  seconds")    #0.057  seconds if sortperm
#display(t1 - t0)    #0.10 seconds  if sortperm
display(aEF.Z)
