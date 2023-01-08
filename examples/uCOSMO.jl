"using COSMO to do LP and QP"
module uCOSMO
# 1000 times slower than Clarabel
# max_iter=N*10000, otherwise Max_iter_reached

using LinearAlgebra, COSMO, SparseArrays, EfficientFrontier
export COSMoLP, COSMoQP
export COSMoCL!


function COSMoLP(E::Vector{T}, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, A::Matrix{T}, b::Vector{T}) where {T}
    N = length(E)
    #iu = u .!= Inf  #Float64(Inf) == BigFloat(Inf)
    iu = u .< Inf
    P = sparse(zeros(T, N, N))
    q = -E
    # Az = b
    c1 = COSMO.Constraint(A, -b, COSMO.ZeroSet)
    # Gz≤g
    c2 = COSMO.Constraint(-G, g, COSMO.Nonnegatives)
    # d≤z≤u
    #c3 = COSMO.Constraint(Matrix{T}(I, N, N), zeros(T, N), COSMO.Box(d, u))
    c3 = COSMO.Constraint([Matrix{T}(I, N, N); -Matrix{T}(I, N, N)[iu, :]], [-d; u[iu]], COSMO.Nonnegatives)
    settings = COSMO.Settings(verbose=false)
    model = COSMO.Model()
    assemble!(model, P, q, [c1; c2; c3], settings=settings)
    COSMO.optimize!(model)
end

function COSMoQP(E::Vector{T}, V::Matrix{T}, mu::T, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, A::Matrix{T}, b::Vector{T}) where {T}
    N = length(E)
    #iu = u .!= Inf
    iu = u .< Inf
    P = sparse(V)
    q = zeros(T, N) #Clarabel need P, q, A, b to be in type T    
    # z′μ=μ Az=b
    c1 = COSMO.Constraint([E'; A], -[mu;b], COSMO.ZeroSet)
    # Gz≤g
    c2 = COSMO.Constraint(-G, g, COSMO.Nonnegatives)
    # d≤z≤u
    #c3 = COSMO.Constraint(Matrix{T}(I, N, N), zeros(T, N), COSMO.Box(d, u))
    c3 = COSMO.Constraint([Matrix{T}(I, N, N); -Matrix{T}(I, N, N)[iu, :]], [-d; u[iu]], COSMO.Nonnegatives)
    settings = COSMO.Settings(verbose=false, max_iter=N*10000)
    model = COSMO.Model()
    assemble!(model, P, q, [c1; c2; c3], settings=settings)
    COSMO.optimize!(model)
end

"""

        COSMoCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where T

compute the Critical Line Segments by COSMO.jl (Conic operator splitting method), for the highest expected return. Save the CL to aCL if done


"""
function COSMoCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = PS
    (; tolS, muShft) = nS

    x = COSMoLP(E, u, d, G, g, A, b)
    #return x
    if x.status != :Solved
        error("Not able to find the max expected return of efficient frontier")
    end
    y = COSMoQP(E, V, x.obj_val * (muShft - 1), u, d, G, g, A, b)
    #return y
    if y.status != :Solved
        error("Not able to find the max expected return efficient portfolio")
    end

    Y = y.s
    iu = u .< Inf
    D = Y[(1:N).+(M+J+1)] .< tolS
    U = falses(N)
    U[iu] = Y[(1:sum(iu)).+(M+J+N+1)] .< tolS
    S = fill(IN, N + J)
    Sv = @view S[1:N]
    Sv[D] .= DN
    Sv[U] .= UP
    if J > 0
        Se = @view S[(1:J).+N]
        Se .= EO
        Og = Y[(1:J).+(M+1)] .> tolS
        Se[Og] .= OE
    end
    computeCL!(aCL, S, PS, nS)
    #return S
end
end