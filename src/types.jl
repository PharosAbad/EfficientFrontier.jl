
# Problem-Algorithm-Solver pattern
# see https://discourse.julialang.org/t/function-depending-on-the-global-variable-inside-module/64322/10


"""

        @enum Status

Status: assets go IN/OUT(DN or UP), or inequalities go binded (EO, as equality) or not (OE, original ineq), with fields:

            IN  #within the lower and upper bound
            DN  #down, lower bound
            UP  #upper bound
            OE  #original <= not binded
            EO  #edge, <= as =

"""
@enum Status begin
    IN
    DN  #down, lower bound
    UP  #upper bound
    OE  #original <= not binded
    EO  #edge, <= as =
end



"""

        struct Event{T<:AbstractFloat}

Events that assets go IN/OUT(DN or UP), or inequalities go binded (EO, as equality) or not (OE, original ineq), with fields:

            From::Status
            To::Status
            id::Int
            L::T

"""
struct Event{T<:AbstractFloat}
    From::Status
    To::Status
    id::Int
    L::T
end



"""

        struct sCL{T<:AbstractFloat}

Critical Line segment, with fields: default T = Float64

            S::Vector{Status}       # (N+J)x1
            K::Integer              #number of IN assets
            L1::T                   #higher lambda
            I1::Vector{Event{T}}    #go in/out events at L1
            L0::T                   #lower lambda
            I0::Vector{Event{T}}    #go in/out events at L0
            alpha::Vector{T}        # K x 1
            beta::Vector{T}         # K x 1

"""
struct sCL{T<:AbstractFloat}    #critical line segment
    S::Vector{Status}   # (N+J)x1
    K::Integer  #number of IN assets
    L1::T #higher lambda
    I1::Vector{Event{T}}   #go in/out events at L1
    L0::T #lower lambda
    I0::Vector{Event{T}}   #go in/out events at L0
    alpha::Vector{T}  # K x 1
    beta::Vector{T}   # K x 1
end

sCL(args...) = sCL{Float64}(args...)


"""

        Problem(E, V; u, d, G, g, A, b, equilibrate)
        Problem(E, V, u; equilibrate)
        Problem(E, V, u, d; equilibrate)
        Problem(E, V, u, d, G, g; equilibrate)
        Problem(E, V, u, d, G, g, A, b; equilibrate)

Setup a Portfolio Selection model: for mean vector E and variance matrix V

```math
	min     (1/2)z′Vz
	s.t.    z′E = μ
	        Az  = b     ∈ R^{M}
	        Gz  ≤ g     ∈ R^{J}
	        d ≤ z ≤ u   ∈ R^{N}
```
Default values: u = +∞, d = 0, G = [], g = [], A = ones(1,N), b = [1], and equilibrate = false

# Example
```
    using EfficientFrontier
    V = [1/100 1/80 1/100
        1/80 1/16 1/40
        1/100 1/40 1/25]    #variance matrix
    E = [109 / 100; 23 / 20; 119 / 100] #mean vector
    P = Problem(E, V)   #Default models, non-shortsale
    aCL = EfficientFrontier.ECL(P)  #compute all Critical Lines
    aEF = eFrontier(aCL, P)     #compute the Efficient Frontier
    mu = (aEF.mu[1]+aEF.mu[end])/2  #mu at the middle of Efficient Frontier
    z = ePortfolio(mu, aEF)     #weight of the efficient portfolio given mu

```

See [`Documentation for EfficientFrontier.jl`](https://github.com/PharosAbad/EfficientFrontier.jl/wiki)

See also [`EfficientFrontier.ECL`](@ref)

"""
struct Problem{T<:AbstractFloat}
    E::Vector{T}
    V::Matrix{T}
    u::Vector{T}
    d::Vector{T}
    G::Matrix{T}
    g::Vector{T}
    A::Matrix{T}
    b::Vector{T}
    N::Int32
    M::Int32
    J::Int32
    equilibrate::Bool    #true: enable equilibration (scale the variance)
    eE::T
    eV::T
end

function Problem(E, V;
    u = fill(Inf, length(E)),
    d = zeros(length(E)),
    G = ones(0, length(E)),
    g = ones(0),
    A = ones(1, length(E)),
    b = ones(1),
    equilibrate=false)

    FloatT = typeof(E).parameters[1]
    N::Int32 = length(E)
    (N, N) == size(V) || throw(DimensionMismatch("incompatible dimension: V"))
    sum(abs.(E)) == 0 && error("mean vector == 0")
    sum(abs.(V)) == 0 && error("variance matrix == 0")
    Eq = copy(vec(E))     #make sure vector and a new copy
    Vq = convert(Matrix{FloatT}, (V+V')/2)   #make sure symmetric
    #@assert det(Vq)>=-sqrt(eps(FloatT)) "variance matrix has negative determinant"
    @assert det(Vq)>=0 "variance matrix has negative determinant"
    eE::FloatT = one(FloatT)
    eV::FloatT = one(FloatT)
    if equilibrate
        eE = maximum(abs.(Eq))   #>0
        eV = maximum(abs.(Vq))   #>0
        Eq ./= eE
        Vq ./= eV
    end
    
    M::Int32 = length(b)
    J::Int32 = length(g)
    N == size(u, 1) || throw(DimensionMismatch("incompatible dimension: u"))
    N == size(d, 1) || throw(DimensionMismatch("incompatible dimension: d"))    
    (J, N) == size(G) || throw(DimensionMismatch("incompatible dimension: G"))
    (M, N) == size(A) || throw(DimensionMismatch("incompatible dimension: A"))

    #norm(A[1,:] .- b[1]) == 0 || error("First equality constraint must be 1'z = 1")
    sum(abs.(A[1, :] .- b[1])) == 0 || error("First equality constraint must be 1'z = 1")

    d[isinf.(d)] .= -1.0    #replace -Inf in lower bound by -1.0
    iu = u .< d
    u[iu] .= d[iu]   #make sure u >= d
    iu = u .< 0
    u[iu] .= 0  #make sure u >= 0


    
    @assert sum(d) < 1 "the sum of downside/lower bound is greater than 1"
    @assert sum(u) > 1 "the sum of upside/higher bound is less than 1"

    Problem{FloatT}(Eq, Vq,
        convert(Vector{FloatT}, copy(vec(u)) ),
        convert(Vector{FloatT}, copy(vec(d)) ),
        convert(Matrix{FloatT}, copy(G)),   #make a copy, just in case it is modified somewhere
        convert(Vector{FloatT}, copy(vec(g)) ),
        convert(Matrix{FloatT}, copy(A)),
        convert(Vector{FloatT}, copy(vec(b)) ), N, M, J, equilibrate, eE, eV)
end

Problem(E, V, u; equilibrate=false) = Problem(E, V; u = u, equilibrate=equilibrate)
Problem(E, V, u, d; equilibrate=false) = Problem(E, V; u = u, d = d, equilibrate=equilibrate)
Problem(E, V, u, d, G, g; equilibrate=false) = Problem(E, V; u = u, d = d, G = G, g = g, equilibrate=equilibrate)
Problem(E, V, u, d, G, g, A, b; equilibrate=false) = Problem(E, V; u = u, d = d, G = G, g = g, A = A, b = b, equilibrate=equilibrate)



"""

        Settings(P::Problem)        The default Settings to given Problem
        Settings(; kwargs...)       The default Settings is set by Float64 type
        Settings{T<:AbstractFloat}(; kwargs...)

kwargs are from the fields of Settings{T<:AbstractFloat} for Float64 and BigFloat

            tol::T         #general scalar
            tolNorm::T     #for norms
            tolS::T        #for s from Clarabel
            tolL::T        #for L
            tolG::T        #for Greeks (beta and gamma)
            muShft::T      #shift the max mu to (1-muShft)*mu    

"""
struct Settings{T<:AbstractFloat}
    tol::T         #general scalar
    tolNorm::T     #for norms
    tolS::T       #for s from Clarabel
    tolL::T        #for L
    tolG::T        #for Greeks (beta and gamma)
    muShft::T      #shift the max mu to (1-muShft)*mu    
end

Settings(; kwargs...) = Settings{Float64}(; kwargs...)
function Settings{Float64}(; tol = 2^-26, 
    tolNorm = 2^-26, 
    tolS = 2^-26, 
    tolL = 2^-26, 
    tolG = 2^-26, 
    muShft = 2^-18)
    Settings{Float64}(tol, tolNorm, tolS, tolL, tolG, muShft)
end

function Settings{BigFloat}(; tol = BigFloat(2)^-76 , 
    tolNorm = BigFloat(2)^-76, 
    #tolS = BigFloat(2)^-76,     #BigFloat(2)^-26
    tolS = BigFloat(2)^-26, #waiting for Clarabel.Settings to update for BigFloat
    tolL = BigFloat(2)^-76, 
    tolG = BigFloat(2)^-76, 
    #muShft = BigFloat(2)^-27)
    muShft = BigFloat(2)^-18)
    Settings{BigFloat}(tol, tolNorm, tolS, tolL, tolG, muShft)
end

function Settings(P::Problem)
    Settings{typeof(P).parameters[1]}()
end



"""
        struct sEF

Efficient Frontier, with fields:

            Ap::Matrix{Float64}     # v=a₂μ²+a₁μ+a₀, each row [a₂ a₁ a₀]
            mu::Vector{Float64}     #higher mean
            sgm::Vector{Float64}    #higher sigma
            Z::Matrix{Float64}      #weights, each corner portfolio in one row
            ic::Vector{Int64}       #id of related critical line

"""
struct sEF    #Efficient Frontier       Float64 is OK
    Ap::Matrix{Float64}  # v=a₂μ²+a₁μ+a₀
    mu::Vector{Float64} #higher mean
    sgm::Vector{Float64}   #higher sigma
    Z::Matrix{Float64}
    ic::Vector{Int64}
end
