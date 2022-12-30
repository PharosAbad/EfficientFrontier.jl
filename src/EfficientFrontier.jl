"Entire Efficient Frontier by Status-Segment Method "
module EfficientFrontier

#using LinearAlgebra, Clarabel, SparseArrays, Combinatorics
using LinearAlgebra, Combinatorics

#export Status, Event, sCL, aCL, computeCL, ECL, CornerP
export Status, Event, sCL, IN, DN, UP, OE, EO
export Settings, Problem, sEF, eFrontier, ePortfolio
export computeCL!, ECL!, cbCL!

include("./types.jl")

include("./criticalLine.jl")

include("./portfolio.jl")

#include("./uClarabel.jl")
#using .uClarabel

include("./uSimplex.jl")
using .uSimplex

end
