"Entire Efficient Frontier by Status-Segment Method "
module EfficientFrontier

export Status, IN, DN, UP, OE, EO, Event, sCL, Problem, Settings, sEF
export computeCL!, ECL!, cbCL!, eFrontier, ePortfolio

using LinearAlgebra, Combinatorics

using LightenQP: LightenQP, solveOOQP
import LightenQP: OOQP, Settings as SettingsQP, fPortfolio
export OOQP, SettingsQP, solveOOQP, fPortfolio #, Solution   #from LightenQP

include("./types.jl")

include("./criticalLine.jl")

include("./portfolio.jl")

include("./Simplex.jl")
using .Simplex: Simplex, SimplexLP
import .Simplex: Settings as SettingsLP
export SettingsLP, SimplexLP

include("./EVdata.jl")

end
