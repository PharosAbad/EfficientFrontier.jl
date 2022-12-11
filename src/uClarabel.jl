"using Clarabel to do LP and QP"
module uClarabel

using LinearAlgebra, Clarabel, SparseArrays     #, EfficientFrontier
export ClarabelQP, ClarabelLP

function ClarabelQP(E::Vector{T}, V::Matrix{T}, mu::T, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, Ae::Matrix{T}, be::Vector{T}) where {T}
    N = length(E)
    iu = u .!= Inf
    P = sparse(V)
    q = zeros(T, N) #Clarabel need P, q, A, b to be in type T
    A = sparse([E'; Ae; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    b = [mu; be; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(1 + length(be)), Clarabel.NonnegativeConeT(length(g) + N + sum(iu))]
    settings = Clarabel.Settings{T}()
    settings.verbose = false
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    Clarabel.solve!(solver)
end

function ClarabelLP(E::Vector{T}, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, Ae::Matrix{T}, be::Vector{T}) where {T}
    N = length(E)
    iu = u .!= Inf  #Float64(Inf) == BigFloat(Inf)
    P = sparse(zeros(T, N, N))
    q = -E
    A = sparse([Ae; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    b = [be; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(length(be)), Clarabel.NonnegativeConeT(length(g) + N + sum(iu))]
    settings = Clarabel.Settings{T}()
    settings.verbose = false
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    Clarabel.solve!(solver)
end


end