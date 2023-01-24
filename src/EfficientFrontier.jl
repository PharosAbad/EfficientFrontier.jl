"Entire Efficient Frontier by Status-Segment Method "
module EfficientFrontier

#using LinearAlgebra, Clarabel, SparseArrays, Combinatorics
using LinearAlgebra, Combinatorics
using LightenQP: LightenQP, solveOOQP
import LightenQP: OOQP, Settings as SettingsQP


#export Status, Event, sCL, aCL, computeCL, ECL, CornerP
export Status, Event, sCL, IN, DN, UP, OE, EO
export Settings, Problem, sEF, eFrontier, ePortfolio
export computeCL!, ECL!, cbCL!
export OOQP, SettingsQP #, Solution, solveOOQP    #from LightenQP

include("./types.jl")

include("./criticalLine.jl")

include("./portfolio.jl")

include("./Simplex.jl")
using .Simplex: Simplex, SimplexLP
import .Simplex: Settings as SettingsLP
export SettingsLP

include("./EVdata.jl")

end
