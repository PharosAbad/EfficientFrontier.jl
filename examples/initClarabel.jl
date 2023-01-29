#When there is an upper bound, we may use the `Clarabel.jl` (an interior point numerical solver) to find a CL
# v1.0  the default Simplex method is the best. Faster speed and less resource than Clarabel

# S&P 500 data, the Covariance matrix is not positive define
#https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda

using EfficientFrontier, LinearAlgebra, CodecXz, Serialization, Downloads, TranscodingStreams
println("\n--- Status-Segment Method vs Markowitz's CLA  ---\n")

#=
#------------ rda
using RData

# download https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda to /tmp
# sp5h = load("/tmp/muS.sp500.rda")
sp500 = Downloads.download("https://gitlab.math.ethz.ch/maechler/CLA/-/raw/master/data/muS.sp500.rda")
sp5h = load(sp500)  #RData

E = values(sp5h["muS.sp500"]["mu"])
V = values(sp5h["muS.sp500"]["covar"])

#display(norm(V - V'))

V = (V+V')/2    #make sure symetry
#maximum(V - V') # 0.0
#display(norm(V - V'))
#------------ rda
=#


#---------- jls.xz
#xzFile = Downloads.download("https://github.com/PharosAbad/PharosAbad.github.io/raw/master/files/sp500.jls.xz")
xzFile = "/tmp/sp500.jls.xz"
if !isfile(xzFile)
    Downloads.download("https://github.com/PharosAbad/PharosAbad.github.io/raw/master/files/sp500.jls.xz", xzFile)
end

io = open(xzFile)
io = TranscodingStream(XzDecompressor(), io)
E = deserialize(io)
V = deserialize(io)
close(io)
#---------- jls.xz



if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    include("./Markowitz.jl")
    using .Markowitz
end

#u = 3/32*ones(length(E))
u = fill(3 / 32, length(E))

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
println("\n--- Status-Segment Method: init by `Clarabel.jl`   ---\n")

settings = uClarabel.SettingsCl(Pu; verbose=false)
ts = @elapsed aCLu = EfficientFrontier.ECL(Pu; init=ClarabelCL!, settings=settings)   #using numerical solver
aEFu = eFrontier(aCLu, Pu)
println("Status-Segment Method:  ", ts, "  seconds")    #0.303391205  seconds


println("\n--- Status-Segment Method: init by default `SimplexCL!`  ---\n")
ts = @elapsed aCL = EfficientFrontier.ECL(Pu)   #v1.0   using Simplex solver
aEF = eFrontier(aCL, Pu)
println("Status-Segment Method:  ", ts, "  seconds")    #0.061639485  seconds
aCLu[7].S == aCL[7].S

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
