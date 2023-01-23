"Entire Efficient Frontier by Status-Segment Method "
module EfficientFrontier

#using LinearAlgebra, Clarabel, SparseArrays, Combinatorics
using LinearAlgebra, Combinatorics
using LightenQP: optionsQP, solutionQP, solveOOQP
import LightenQP: OOQP

#export Status, Event, sCL, aCL, computeCL, ECL, CornerP
export Status, Event, sCL, IN, DN, UP, OE, EO
export Settings, Problem, sEF, eFrontier, ePortfolio
export computeCL!, ECL!, cbCL!
export OOQP, optionsQP, solutionQP, solveOOQP    #from LightenQP

include("./types.jl")

include("./criticalLine.jl")

include("./portfolio.jl")

include("./uSimplex.jl")
using .uSimplex

include("./EVdata.jl")

end
