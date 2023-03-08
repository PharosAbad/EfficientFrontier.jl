__precompile__()

"Entire Efficient Frontier by Status-Segment Method "
module EfficientFrontier

export Status, IN, DN, UP, OE, EO, Event, sCL, Problem, Settings, sEF, SettingsQP, SettingsLP
export computeCL!, ECL!, cbCL!, eFrontier, ePortfolio, mu2L, L2mu

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
export solveASQP, asQP, asCL!, getSx


# #=
Base.precompile(Tuple{Core.Typeof(ECL), Problem{Float64}})
Base.precompile(Tuple{Core.Typeof(eFrontier), sCL{Float64}, Problem{Float64}})
Base.precompile(Tuple{Core.Typeof(ePortfolio), sEF, Float64})
Base.precompile(Tuple{Core.Typeof(ECL!), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{Core.Typeof(computeCL!), Vector{sCL{Float64}}, Vector{Status}, Problem{Float64}, Settings{Float64}})
Base.precompile(Tuple{Core.Typeof(SimplexCL!), Vector{sCL{Float64}}, Problem{Float64}})
#Base.precompile(Tuple{Core.Typeof(fPortfolio), Problem{Float64}, Float64})
Base.precompile(Tuple{Core.Typeof(SimplexLP), Problem{Float64}})
#Base.precompile(Tuple{Core.Typeof(Simplex.cDantzigLP), Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Status}})
# =#




#force precompile
#E, V = EVdata(:Abad, false)
#= E, V = EVdata(:Ungil, false)
P = Problem(E, V)
aCL = ECL(P)
aEF = eFrontier(aCL, P)
nothing =#



#=
using SnoopPrecompile
@precompile_setup begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    E, V = EVdata(:Abad, false)
    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        P = Problem(E, V)
        aCL = ECL(P)
        aEF = eFrontier(aCL, P)
    end
end
=#

end
