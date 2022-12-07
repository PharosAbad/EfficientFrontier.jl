"Full Efficient Frontier by connecting Critical Line Segments"
module EfficientFrontier

using LinearAlgebra, Clarabel, SparseArrays

#export Status, Event, sCL, aCL, computeCL, ECL, CornerP
export Status, Event, sCL, IN, DN, UP, OE, EO
export Settings, Problem, sEF, eFrontier, ePortfolio
export computeCL!, ECL!

include("./types.jl")

include("./criticalLine.jl")

include("./portfolio.jl")

end
