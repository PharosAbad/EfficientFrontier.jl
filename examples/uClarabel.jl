"using Clarabel (an interior point numerical solver) to do LP and QP"
module uClarabel
# promising, but too BIG, if the dep Pardiso can be removed, it will be perfect

using LinearAlgebra, Clarabel, SparseArrays, EfficientFrontier
export ClarabelQP, ClarabelLP
export ClarabelCL!

function SettingsCl(PS::Problem{T}; kwargs...) where {T}
    if isempty(kwargs)  #default setting
        return Clarabel.Settings{T}(tol_gap_abs=2^-52, tol_gap_rel=2^-47, verbose = false)
    end
    Clarabel.Settings{T}(; kwargs...)
end

function SettingsCl(; kwargs...)
    if isempty(kwargs)  #default setting
        return Clarabel.Settings(tol_gap_abs=2^-52, tol_gap_rel=2^-47, verbose = false)
    end    
    Clarabel.Settings(; kwargs...)
end

function ClarabelQP(PS::Problem{T}, mu::T; settings= SettingsCl()) where {T}
#function ClarabelQP(PS::Problem{T}, mu::T; settings= SettingsCl(PS; verbose = false)) where {T}
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


# LVEP (Lowest Variance Efficient Portfolio) == Global Minimum Variance Portfolio (GMVP)
function ClarabelQP(PS::Problem{T}; settings= SettingsCl()) where {T}
#function ClarabelQP(PS::Problem{T}; settings= SettingsCl(PS; verbose = false)) where {T}
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


#find the highest mean: -x.obj_val
function ClarabelLP(PS::Problem{T}; settings= SettingsCl()) where {T}
#function ClarabelLP(PS::Problem{T}; settings= SettingsCl(PS; verbose = false)) where {T}
    (; E, u, d, G, g, A, b, N, M, J) = PS
    iu = findall(u .< Inf)
    P = spzeros(T, N, N)
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
function ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl(), kwargs...) where {T}
#function ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl(PS; verbose = false), kwargs...) where {T}
    
    x = ClarabelLP(PS; settings=settings)
    if Int(x.status) != 1   #SOLVED
        error("Not able to find the expected return of HMFP (Highest Mean Frontier Portfolio)")
    end
    mu = -x.obj_val
    shft =  nS.muShft
    if  mu < -1 || mu > 1
        shft *= abs(mu)
    end
    mu -= shft
    #y = ClarabelQP(PS, x.obj_val * (nS.muShft - 1); settings=settings)
    y = ClarabelQP(PS, mu; settings=settings)
    if Int(y.status) != 1   #SOLVED
        error("Not able to find a muShft to the HMFP (Highest Mean Frontier Portfolio)")
    end

    Y = y.s[PS.M+2:end]
    S = EfficientFrontier.getS(Y, PS, nS.tolS) 
    

    #=
    # GMVP, fail to get correct S most of time. See https://github.com/oxfordcontrol/Clarabel.jl/issues/109 Corner Portfolio Blur 
    y = ClarabelQP(PS; settings=settings)   #may fail, leave it to computeCL!
    #= if Int(y.status) != 1   #SOLVED
        error("Clarabel: Not able to find the Lowest Variance Efficient Portfolio")
    end =#
    Y = y.s[PS.M+1:end]
    S = EfficientFrontier.getS(Y, PS, nS.tolS*128)  =#


    computeCL!(aCL, S, PS, nS)
end

end