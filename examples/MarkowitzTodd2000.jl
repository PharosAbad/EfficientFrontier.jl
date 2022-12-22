# lower and upper bounds, 2 inequality constraints

using EfficientFrontier
#using Pkg; Pkg.add("EfficientFrontier")
#using Pkg; Pkg.add("https://github.com/PharosAbad/EfficientFrontier.jl.git")

# from  CLA-Data.csv https://github.com/mdengler/cla
V = [0.40755159 0.03175842 0.05183923 0.05663904 0.0330226 0.00827775 0.02165938 0.01332419 0.0343476 0.02249903
   0.03175842 0.9063047 0.03136385 0.02687256 0.01917172 0.00934384 0.02495043 0.00761036 0.02874874 0.01336866
   0.05183923 0.03136385 0.19490901 0.04408485 0.03006772 0.01322738 0.03525971 0.0115493 0.0427563 0.02057303
   0.05663904 0.02687256 0.04408485 0.19528471 0.02777345 0.00526665 0.01375808 0.00780878 0.02914176 0.01640377
   0.0330226 0.01917172 0.03006772 0.02777345 0.34059105 0.00777055 0.02067844 0.00736409 0.02542657 0.01284075
   0.00827775 0.00934384 0.01322738 0.00526665 0.00777055 0.15983874 0.02105575 0.00518686 0.01723737 0.00723779
   0.02165938 0.02495043 0.03525971 0.01375808 0.02067844 0.02105575 0.68056711 0.01377882 0.04627027 0.01926088
   0.01332419 0.00761036 0.0115493 0.00780878 0.00736409 0.00518686 0.01377882 0.95526918 0.0106553 0.00760955
   0.0343476 0.02874874 0.0427563 0.02914176 0.02542657 0.01723737 0.04627027 0.0106553 0.31681584 0.01854318
   0.02249903 0.01336866 0.02057303 0.01640377 0.01284075 0.00723779 0.01926088 0.00760955 0.01854318 0.11079287]
E = [1.175, 1.19, 0.396, 1.12, 0.346, 0.679, 0.089, 0.73, 0.481, 1.08]

#EfficientFrontier
E0 = vec(E)
V0 = copy(V)
N = length(E)

A = ones(1, N)
b = ones(1)
d = zeros(N)
d[1] = 0.1
d[5] = 0.1
u = 0.3 * ones(N)
g = [-0.2; 0.5]
G = zeros(length(g), N)
G[1, 1:3] .= -1.0
G[1, 4] = -0.5
G[2, 4] = 0.5
G[2, 5:7] .= 1.0

println("--- Markowitz and Todd (2000), chapter 13, pp.337 --- ")
#=
if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
   include("./Markowitz.jl")
   using .Markowitz
end
m = markowitz(E, V, lower=d, upper=u)
unit_sum(m)
add_constraint(m, G[1,:], '<', g[1])
add_constraint(m, G[2,:], '<', g[2])
ts = @elapsed f = frontier(m)
println("Markowitz CLA:  ", ts, "  seconds")    #0.0007 seconds
display(f.weights)
=#


P = Problem(E0, V0, u, d, G, g)
ts = @elapsed aCL = EfficientFrontier.ECL(P)
aEF = eFrontier(aCL, P)
println("connecting Critical Line Segments:  ", ts, "  seconds")   #0.0006 seconds
display(aEF.Z)

# the CLA in Markowitz and Todd (2000) assumes either one IN or one OUT.
# Our code handles two or more IN and/or OUT.  The first CL in this example has: 
# Event[Event(UP, IN, 2, 2.494839663636366), Event(DN, IN, 10, 2.494839663636377)]. 
# Because the first Corner Portfolio is on the boundary, there is no IN assets
#aCL[1].I1[1].L - aCL[1].I1[2].L    #-1.1102230246251565e-14

#=
#BigFloat
Pb = Problem(convert(Vector{BigFloat}, E0), V0, u, d, G, g)
ts = @elapsed aCLb = EfficientFrontier.ECL(Pb; numSettings = Settings{BigFloat}(tolL = BigFloat(2)^-51))
#aCLc = EfficientFrontier.ECL(Pb; init=EfficientFrontier.ClarabelCL!, numSettings = Settings{BigFloat}(tolL = BigFloat(2)^-51))
aEFb = eFrontier(aCLb, Pb)
println("BigFloat:  ", ts, "  seconds")   #0.020 seconds
#display(aEFb.Z)
#maximum(abs.(aEFb.Z-aEF.Z))   #3.344546861683284e-15

#=
aCLb[1].I1
2-element Vector{Event{BigFloat}}:
 Event{BigFloat}(UP, IN, 2, 2.494839663636366245529450347001271614280350227856738443120461658250300005883108)
 Event{BigFloat}(DN, IN, 10, 2.494839663636366495420530397802951875034152818858742312229266642281434186751967)
=#
#aCLb[1].I1[1].L - aCLb[1].I1[2].L  #-2.498910800508016802607538025910020038691088049840311341808673741032753977364765e-16
#for BigFloat, default tolL = BigFloat(2)^-76 = 1.32348898008484427979425390731194056570529937744140625e-23
=#

#BigFloat from the source data: convert to Int (rationals) then to BigFloat to improve precision
V0 = convert(Matrix{Int64},round.(V*10^8))
E0 = convert(Vector{Int64},round.(E*10^3))
d0 = convert(Vector{Int64},round.(d*10))
u0 = convert(Vector{Int64},round.(u*10))
g0 = convert(Vector{Int64},round.(g*10))
G0 = convert(Matrix{Int64},round.(G*10))
Vb = V0 /BigFloat(10^8)
Eb = E0 /BigFloat(10^3)
db = d0 /BigFloat(10)
ub = u0 /BigFloat(10)
gb = g0 /BigFloat(10)
Gb = G0 /BigFloat(10)
#println("\n xxx  ", maximum(abs.(Eb-E)))
println("\n E[4]=1.12,  but BigFloat(\"1.12\")-E[4] is: ", BigFloat("1.12")-E[4])
println("now convert the raw data to BigFloat to improve precision (raw data to BigFloat, not raw data to Float64 then to BigFloat)")


#S = Status[UP, IN, DN, UP, DN, DN, DN, DN, DN, IN, OE, OE]
S = aCL[1].S
Pb = Problem(Eb, Vb, ub, db, Gb, gb)
aCLb = Vector{sCL{typeof(Pb).parameters[1]}}(undef, 0)
computeCL!(aCLb, S, Pb, Settings(Pb))
#computeCL!(aCLb, S, Pb, Settings{BigFloat}(tolL = BigFloat(2)^-51)) #if the data is first truncated by Float64, adjust the TolL


println("\n Assets 2 and 10 go IN at the same time, NOT the case that ONE at each time")
display(aCLb[1].I1)
println("using BigFloat, the difference on L: ", aCLb[1].I1[1].L -aCLb[1].I1[2].L)
#aCLb[1].I1[1].L -aCLb[1].I1[2].L        #-1.450876317255866697064907112950467127947488061225295272683982182988323422931288e-75



println("
Because the first Corner Portfolio is on the boundary, there is no IN assets (K=0).  there must be K>=M for each CL.
   Hence,  K>=2 assets go IN if M>=2 (M=1 in this example)

the CLA in Markowitz and Todd (2000) assumes either one IN or one OUT (if and only if one asset changes state).
Our code handles two or more IN and/or OUT (allow any number of changes concurrently).
")
