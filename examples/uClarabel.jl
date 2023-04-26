"using Clarabel (an interior point numerical solver) to do LP and QP"
module uClarabel
# promising, but too BIG, if the dep Pardiso can be removed, it will be perfect

using LinearAlgebra, Clarabel, SparseArrays, EfficientFrontier
export ClarabelQP, ClarabelLP
export ClarabelCL!

function SettingsCl(PS::Problem{T}; kwargs...) where {T}
    if isempty(kwargs)  #default setting
        return Clarabel.Settings{T}(tol_gap_abs=2^-52, tol_gap_rel=2^-47, verbose=false)
    end
    Clarabel.Settings{T}(; kwargs...)
end

function SettingsCl(; kwargs...)
    if isempty(kwargs)  #default setting
        return Clarabel.Settings(tol_gap_abs=2^-52, tol_gap_rel=2^-47, verbose=false)
    end
    Clarabel.Settings(; kwargs...)
end

function ClarabelQP(mu::T, PS::Problem{T}; settings=SettingsCl()) where {T}
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

#=
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
=#

function ClarabelQP(PS::Problem{T}, L::T=0.0; settings=SettingsCl()) where {T}

    if isinf(L)
        min = L == Inf ? false : true
        x = ClarabelLP(PS; settings=settings, min=min)
        return ClarabelQP(x.obj_val, PS; settings=settings)
    end

    (; E, V, u, d, G, g, A, b, N, M, J) = PS
    iu = findall(u .< Inf)
    P = sparse(V)
    #q = zeros(T, N) #Clarabel need P, q, A, b to be in type T
    Ac = sparse([A; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    bc = [b; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(M), Clarabel.NonnegativeConeT(J + N + length(iu))]
    solver = Clarabel.Solver{T}()
    if L == 0
        qq = zeros(T, N)
    else
        qq = -L * E
    end

    #=
    if L == 0
        Clarabel.setup!(solver, P, zeros(T, N), Ac, bc, cones, settings)
        return Clarabel.solve!(solver)
    end
    if isfinite(L)
        Clarabel.setup!(solver, P, -L * E, Ac, bc, cones, settings)
        return Clarabel.solve!(solver)
    end
    =#
    Clarabel.setup!(solver, P, qq, Ac, bc, cones, settings)
    return Clarabel.solve!(solver)
end


#find the highest mean: -x.obj_val
function ClarabelLP(PS::Problem{T}; settings=SettingsCl(), min=true) where {T}
    #function ClarabelLP(PS::Problem{T}; settings= SettingsCl(PS; verbose = false)) where {T}
    (; E, u, d, G, g, A, b, N, M, J) = PS
    iu = findall(u .< Inf)
    P = spzeros(T, N, N)
    #q = -E
    #sgn = min == true ? 1 : -1
    q = min ? E : -E
    Ac = sparse([A; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    bc = [b; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(M), Clarabel.NonnegativeConeT(J + N + length(iu))]
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, Ac, bc, cones, settings)
    #Clarabel.setup!(solver, P, sgn * E, Ac, bc, cones, settings)
    x = Clarabel.solve!(solver)
    if !min
        x.obj_val = -x.obj_val
    end
    return x
end


function getS(Y, PS, tolS)
    #from slack variables
    (; u, N, J) = PS
    iu = findall(u .< Inf)
    D = Y[(1:N).+J] .< tolS
    U = falses(N)
    U[iu] = Y[(1:length(iu)).+(J+N)] .< tolS
    S = fill(IN, N + J)
    Sv = @view S[1:N]
    Sv[D] .= DN
    Sv[U] .= UP
    if J > 0
        Se = @view S[(1:J).+N]
        Se .= EO
        Og = Y[1:J] .> tolS
        Se[Og] .= OE
    end
    return S
end


"""

        ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl, kwargs...) where T

compute the Critical Line Segments by Clarabel.jl (Interior Point QP), for the highest expected return. Save the CL to aCL if done


"""
function ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl(), kwargs...) where {T}
    #function ClarabelCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsCl(PS; verbose = false), kwargs...) where {T}

    x = ClarabelLP(PS; settings=settings, min=false)
    if Int(x.status) != 1   #SOLVED
        error("Not able to find the expected return of HMFP (Highest Mean Frontier Portfolio)")
    end
    #mu = -x.obj_val
    mu = x.obj_val
    shft = nS.muShft
    if mu < -1 || mu > 1
        shft *= abs(mu)
    end
    mu -= shft
    y = ClarabelQP(mu, PS; settings=settings)
    if Int(y.status) != 1   #SOLVED
        error("Not able to find a muShft to the HMFP (Highest Mean Frontier Portfolio)")
    end

    Y = y.s[PS.M+2:end]
    #S = EfficientFrontier.getS(Y, PS, nS.tolS)
    S = getS(Y, PS, nS.tolS)


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