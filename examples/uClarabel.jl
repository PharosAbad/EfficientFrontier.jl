"using Clarabel to do LP and QP"
module uClarabel
# promising, but too BIG, if the dep Pardiso can be removed, it will be perfect

using LinearAlgebra, Clarabel, SparseArrays, EfficientFrontier
export ClarabelQP, ClarabelLP
export ClarabelCL!

function SettingsCl(PS::Problem{T}; kwargs...) where {T}
    Clarabel.Settings{T}(; kwargs...)
end

function ClarabelQP(PS::Problem{T}, mu::T; settings= SettingsCl(PS; verbose = false)) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = PS
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


#find the highest mean
function ClarabelLP(PS::Problem{T}; settings= SettingsCl(PS; verbose = false)) where {T}
    (; E, u, d, G, g, A, b, N, M, J) = PS
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


"""

        ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl, kwargs...) where T

compute the Critical Line Segments by Clarabel.jl (Interior Point QP), for the highest expected return. Save the CL to aCL if done


"""
function ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl(PS; verbose = false), kwargs...) where {T}
    x = ClarabelLP(PS; settings=settings)
    if Int(x.status) != 1   #SOLVED
        error("Not able to find the expected return of HMFP (Highest Mean Frontier Portfolio)")
    end
    y = ClarabelQP(PS, x.obj_val * (nS.muShft - 1); settings=settings)
    if Int(y.status) != 1   #SOLVED
        error("Not able to find a muShft to the HMFP (Highest Mean Frontier Portfolio)")
    end

    Y = y.s[PS.M+2:end]
    S = EfficientFrontier.getS(Y, PS, nS.tolS)
    computeCL!(aCL, S, PS, nS)
end

end