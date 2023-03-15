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
#precompile(Tuple{Type{Problem{T} where T<:AbstractFloat}, Array{Float64, 1}, Array{Float64, 2}})
#precompile(Tuple{Type{Problem{T} where T<:AbstractFloat}, Array{Float64, 1}, Array{Float64, 2}, Array{Float64, 1}})
Base.precompile(Tuple{typeof(ECL), Problem{Float64}})
Base.precompile(Tuple{typeof(SimplexCL!), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{typeof(SimplexLP), Problem{Float64}})
Base.precompile(Tuple{typeof(computeCL!), Vector{sCL{Float64}}, Vector{Status}, Problem{Float64}, Settings{Float64}})
Base.precompile(Tuple{typeof(ECL!), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{typeof(Simplex.cDantzigLP), Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Status}})
Base.precompile(Tuple{typeof(eFrontier), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{typeof(ePortfolio), sEF, Float64})
#Base.precompile(Tuple{typeof(fPortfolio), Problem{Float64}, Float64})
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
    u = fill(Inf, length(E))
    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        P = Problem(E, V)
        Pc = Problem(E, V, u)
        aCLc = sCL(P)
        SimplexCL!(aCLc, P; nS=Settings(P), settings=SettingsQP(P), settingsLP=SettingsLP(P))
        aCL = ECL(P)
        aEF = eFrontier(aCL, P)
    end
end
 =#

end
