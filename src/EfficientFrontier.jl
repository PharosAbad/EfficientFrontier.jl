#__precompile__()

"Entire Efficient Frontier by Status-Segment Method "
module EfficientFrontier

export Status, IN, DN, UP, OE, EO, Event, sCL, Problem, Settings, sEF, SettingsQP, SettingsLP
export computeCL!, ECL!, cbCL!, eFrontier, ePortfolio

using LinearAlgebra, Combinatorics

using LightenQP: LightenQP, solveOOQP
import LightenQP: OOQP, fPortfolio #, Settings as SettingsQP
export OOQP, solveOOQP, fPortfolio #, Solution   #from LightenQP

include("./types.jl")

include("./criticalLine.jl")

include("./portfolio.jl")

include("./Simplex.jl")
#using .Simplex: Simplex, SimplexLP
using .Simplex
#import .Simplex: Settings as SettingsLP
export SimplexLP #

include("./EVdata.jl")


include("./ASQP.jl")
using .ASQP
export solveASQP, asQP, asCL!

#=
Base.precompile(Tuple{Core.Typeof(ECL), Problem{Float64}})
Base.precompile(Tuple{Core.Typeof(eFrontier), sCL{Float64}, Problem{Float64}})
Base.precompile(Tuple{Core.Typeof(ePortfolio), sEF, Float64})
Base.precompile(Tuple{Core.Typeof(ECL!), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{Core.Typeof(computeCL!), Vector{sCL{Float64}}, Vector{Status}, Problem{Float64}, Settings{Float64}})
Base.precompile(Tuple{Core.Typeof(SimplexCL!), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{Core.Typeof(fPortfolio), Problem{Float64}, Float64})
=#


#force precompile
E, V = EVdata(:Abad, false)
P = Problem(E, V)
aCL = ECL(P)
aEF = eFrontier(aCL, P)
nothing

end
