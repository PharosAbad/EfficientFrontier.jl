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
using .Simplex: Simplex, SimplexLP
#import .Simplex: Settings as SettingsLP
export SimplexLP #

include("./EVdata.jl")

end
