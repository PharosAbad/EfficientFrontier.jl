"using Clarabel to do LP and QP"
module uClarabel
# promising, but too BIG, if the dep Pardiso can be removed, it will be perfect

using LinearAlgebra, Clarabel, SparseArrays, EfficientFrontier
export ClarabelQP, ClarabelLP
export ClarabelCL!

function SettingsCl(PS::Problem{T}; kwargs...) where {T}
#function SettingsCl{T}(; kwargs...) where {T}
    settings = Clarabel.Settings{T}(; kwargs...)
    #settings.verbose = false
end

function ClarabelQP(PS::Problem{T}, mu::T; settings= SettingsCl(PS; verbose = false)) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = PS
    #N = length(E)
    iu = findall(u .< Inf)
    P = sparse(V)
    q = zeros(T, N) #Clarabel need P, q, A, b to be in type T
    Ac = sparse([E'; A; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    bc = [mu; b; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(1 + M), Clarabel.NonnegativeConeT(J + N + length(iu))]    
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, Ac, bc, cones, settings)
    Clarabel.solve!(solver)
end

function ClarabelQP0(E::Vector{T}, V::Matrix{T}, mu::T, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, Ae::Matrix{T}, be::Vector{T}) where {T}
    N = length(E)
    #iu = u .!= Inf
    iu = u .< Inf
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

# Global Minimum Variance Portfolio (GMVP)
function ClarabelQP(PS::Problem{T}; settings= SettingsCl(PS; verbose = false)) where {T}
    (; V, u, d, G, g, A, b, N, M, J) = PS
    iu = findall(u .< Inf)
    P = sparse(V)
    q = zeros(T, N) #Clarabel need P, q, A, b to be in type T
    Ac = sparse([A; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    bc = [b; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(M), Clarabel.NonnegativeConeT(J + N + length(iu))]
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, Ac, bc, cones, settings)
    Clarabel.solve!(solver)
end
function ClarabelQP0(V::Matrix{T}, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, Ae::Matrix{T}, be::Vector{T}) where {T}
    N = length(d)
    iu = findall(u .< Inf)
    P = sparse(V)
    q = zeros(T, N) #Clarabel need P, q, A, b to be in type T
    A = sparse([Ae; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    b = [be; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(length(be)), Clarabel.NonnegativeConeT(length(g) + N + length(iu))]
    settings = Clarabel.Settings{T}()
    settings.verbose = false
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    Clarabel.solve!(solver)
end


function ClarabelLP(PS::Problem{T}; settings= SettingsCl(PS; verbose = false)) where {T}
    (; E, u, d, G, g, A, b, N, M, J) = PS
    #N = length(E)
    iu = findall(u .< Inf)
    P = sparse(zeros(T, N, N))
    q = -E
    Ac = sparse([A; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    bc = [b; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(M), Clarabel.NonnegativeConeT(J + N + length(iu))]    
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, Ac, bc, cones, settings)
    Clarabel.solve!(solver)
end
function ClarabelLP0(E::Vector{T}, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, Ae::Matrix{T}, be::Vector{T}) where {T}
    N = length(E)
    #iu = u .!= Inf  #Float64(Inf) == BigFloat(Inf)
    iu = u .< Inf
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

"""

        ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl, kwargs...) where T

compute the Critical Line Segments by Clarabel.jl (Interior Point QP), for the highest expected return. Save the CL to aCL if done


"""
function ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl(PS; verbose = false), kwargs...) where {T}
#function ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), kwargs...) where {T}
    #(; E, V, u, d, G, g, A, b, N, M, J) = PS
    #(; E, V, u, d, G, g, A, b, N) = PS
    #(; tolS, muShft) = nS

    #x = ClarabelLP(E, u, d, G, g, A, b)
    x = ClarabelLP(PS; settings=settings)
    if Int(x.status) != 1   #SOLVED
        error("Not able to find the expected return of HMFP (Highest Mean Frontier Portfolio)")
    end
    #y = ClarabelQP(E, V, x.obj_val * (muShft - 1), u, d, G, g, A, b)
    y = ClarabelQP(PS, x.obj_val * (nS.muShft - 1); settings=settings)
    if Int(y.status) != 1   #SOLVED
        error("Not able to find a muShft to the HMFP (Highest Mean Frontier Portfolio)")
    end

    #= Y = y.s
    #iu = u .!= Inf
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
    end =#
    Y = y.s[PS.M+2:end]
    S = EfficientFrontier.getS(Y, PS, nS.tolS)
    computeCL!(aCL, S, PS, nS)
end

end