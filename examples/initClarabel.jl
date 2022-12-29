#When there is an upper bound, we may use the `Clarabel.jl` (an interior point numerical solver) to find a CL
# v1.0  the default Simplex method is the best. Faster speed and less resource than Clarabel

# S&P 500 data, the Covariance matrix is not positive define
#https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda

using CodecXz
using RData
using EfficientFrontier
using LinearAlgebra
using Serialization


println("\n pls download https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda to /tmp \n")
sp5h = load("/tmp/muS.sp500.rda")

E = values(sp5h["muS.sp500"]["mu"])
V = values(sp5h["muS.sp500"]["covar"])
V = (V+V')/2


println("\n--- connecting Critical Line Segments vs Markowitz's CLA  ---\n")

if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    include("./Markowitz.jl")
    using .Markowitz
end

u = 3/32*ones(length(E))

m = markowitz(E, V; upper=u)
unit_sum(m) # total weight = 100%
ts = @elapsed f = frontier(m)
println("Markowitz CLA:  ", ts, "  seconds")     #0.14 seconds



Pu = Problem(E, V, u)

#DO NOT do this, unless you have a quantum computer
#ts = @elapsed aCLu = EfficientFrontier.ECL(Pu; init=cbCL!)
# with upper bound, if N is large, e.g., N>30  (the solution space may reach 3^N, and at least 2^N for combinatorial search)

if length(filter((x) -> x == :uClarabel, names(Main, imported=true))) == 0
    include("./uClarabel.jl")
    using .uClarabel
end
println("\n--- connecting Critical Line Segments: init by `Clarabel.jl`   ---\n")

ts = @elapsed aCLu = EfficientFrontier.ECL(Pu; init=ClarabelCL!)   #using numerical solver
aEFu = eFrontier(aCLu, Pu)
println("connecting Critical Line Segments:  ", ts, "  seconds")    #0.303391205  seconds


println("\n--- connecting Critical Line Segments: init by default `SimplexCL!`  ---\n")
ts = @elapsed aCL = EfficientFrontier.ECL(Pu)   #v1.0   using Simplex solver
aEF = eFrontier(aCL, Pu)
println("connecting Critical Line Segments:  ", ts, "  seconds")    #0.061639485  seconds


#=
using LinearAlgebra
aCL = aCLu
tolNorm = 2^-26

nL = lastindex(aCL)
W = trues(nL)
K = 0
for k in eachindex(aCL)
    if norm(aCL[k].beta) < tolNorm
        W[k] = false
    end
    global K = max(K, aCL[k].K)
end

display( findall(.!W) )
#CL 1 is singular, the starting point of frontier, CL 9 is also a singular, hence a kink

println("CLA  ", size(f.weights,1), "\ncCLS ", size(aEFu.Z,1), "\nmax K: ", K)
=#
