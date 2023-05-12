#__precompile__()

"Entire Efficient Frontier by Status-Segment Method "
module EfficientFrontier

using LinearAlgebra, Combinatorics
#import StatusSwitchingQP as SSQ

using StatusSwitchingQP: Status, IN, DN, UP, OE, EO, Event, getRowsGJr, solveQP
using StatusSwitchingQP: StatusSwitchingQP as SS
import StatusSwitchingQP: SimplexLP, QP

export Status, IN, DN, UP, OE, EO, Event, sCL, Problem, Settings, sEF, SettingsQP, SettingsLP
export computeCL!, ECL!, cbCL!, eFrontier, ePortfolio, mu2L, L2mu
export SimplexLP, QP

include("./types.jl")

include("./SS.jl")

include("./criticalLine.jl")

include("./portfolio.jl")

#=
include("./Simplex.jl")
#using .Simplex: Simplex, SimplexLP
using .Simplex
#import .Simplex: Settings as SettingsLP
export SimplexLP #
=#

include("./EVdata.jl")

#=
include("./SSQP.jl")
using .SSQP
export solveQP, QP
=#


#=
#precompile(Tuple{Type{Problem{T} where T<:AbstractFloat}, Array{Float64, 1}, Array{Float64, 2}})
#precompile(Tuple{Type{Problem{T} where T<:AbstractFloat}, Array{Float64, 1}, Array{Float64, 2}, Array{Float64, 1}})
Base.precompile(Tuple{typeof(ECL), Problem{Float64}})
Base.precompile(Tuple{typeof(SimplexCL!), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{typeof(SimplexLP), Problem{Float64}})
Base.precompile(Tuple{typeof(computeCL!), Vector{sCL{Float64}}, Vector{Status}, Problem{Float64}, Settings{Float64}})
Base.precompile(Tuple{typeof(ECL!), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{typeof(Simplex.cDantzigLP), Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Status}})
Base.precompile(Tuple{typeof(eFrontier), Vector{sCL{Float64}}, Problem{Float64}})
Base.precompile(Tuple{typeof(ePortfolio), Float64, sEF})
#Base.precompile(Tuple{typeof(fPortfolio), Float64, Problem{Float64}})
=#

using PrecompileTools
#PrecompileTools.verbose[] = true

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster
    E, V = EVdata(:Abad, false)
    #u = fill(Inf, length(E))
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        P = Problem(E, V)
        #Pc = Problem(E, V, u)
        #aCLc = sCL(P)
        #SimplexCL!(aCLc, P; nS=Settings(P), settings=SettingsQP(P), settingsLP=SettingsLP(P))
        aCL = ECL(P)
        aEF = eFrontier(aCL, P)
        #t = computeCL!(aCLc, aCL[2].S, P, Settings(P))
    end
end




end
