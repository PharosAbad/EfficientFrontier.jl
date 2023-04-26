"using OSQP (Operator Splitting Quadratic Program) numerical solver to do LP and QP"
module uOSQP

using LinearAlgebra, OSQP, SparseArrays, EfficientFrontier
export OpSpQP, OpSpLP, OpSpCL!

function SettingsOS(; kwargs...)
    if isempty(kwargs)  #default setting
        return SettingsOS(; eps_abs=2^-26, eps_rel=2^-26, eps_prim_inf=2^-26, eps_dual_inf=2^-26, verbose=false, polish=true, max_iter=771000)
    end
    kwargs
end


function OpSpQP(PS::Problem{T}, L::T=0.0; settings=SettingsOS()) where {T}
    if isinf(L)
        min = L == Inf ? false : true
        x = OpSpLP(PS; settings=settings, min=min)
        return OpSpQP(x.info.obj_val, PS; settings=settings)
    end

    (; E, V, u, d, G, g, A, b, N, J) = PS
    P = sparse(V)
    Ao = sparse([A; G; Matrix{T}(I, N, N)])
    uo = [b; g; u]
    lo = [b; fill(-Inf, J); d]
    model = OSQP.Model()
    #=
    if L == 0
        OSQP.setup!(model; P=P, q=zeros(T, N), A=Ao, l=lo, u=uo, settings...) # `settings...` to iterate pairs in `settings`
        return OSQP.solve!(model)
    end
    if isfinite(L)
        OSQP.setup!(model; P=P, q=-L * E, A=Ao, l=lo, u=uo, settings...) # `settings...` to iterate pairs in `settings`
        return OSQP.solve!(model)
    end
    sgn = L == Inf ? -1 : 1
    OSQP.setup!(model; P=P, q=sgn * E, A=Ao, l=lo, u=uo, settings...) # `settings...` to iterate pairs in `settings`
    =#
    if L == 0
        qq = zeros(T, N)
    else
        qq = -L * E
    end
    OSQP.setup!(model; P=P, q=qq, A=Ao, l=lo, u=uo, settings...) # `settings...` to iterate pairs in `settings`
    return OSQP.solve!(model)
end


#find the highest mean: -x.info.obj_val
function OpSpLP(PS::Problem{T}; settings=SettingsOS(), min=true) where {T}
    (; E, u, d, G, g, A, b, N, J) = PS
    P = spzeros(T, N, N)
    #q = -E
    q = min ? E : -E
    Ao = sparse([A; G; Matrix{T}(I, N, N)])
    uo = [b; g; u]
    lo = [b; fill(-Inf, J); d]
    model = OSQP.Model()
    OSQP.setup!(model; P=P, q=q, A=Ao, l=lo, u=uo, settings...)
    x = OSQP.solve!(model)
    if !min
        x.info.obj_val = -x.info.obj_val
    end
    return x
end

function OpSpQP(mu::T, PS::Problem{T}; settings=SettingsOS()) where {T}
#function OpSpQP(E::Vector{T}, V::Matrix{T}, mu::T, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, A::Matrix{T}, b::Vector{T}) where {T}
    (; E, V, u, d, G, g, A, b, N, J) = PS
    P = sparse(V)
    q = zeros(T, N)
    Ao = sparse([E'; A; G; Matrix{T}(I, N, N)])
    uo = [mu; b; g; u]
    lo = [mu; b; fill(-Inf, J); d]
    model = OSQP.Model()
    OSQP.setup!(model; P=P, q=q, A=Ao, l=lo, u=uo, settings...)
    OSQP.solve!(model)
end


"""

        OpSpCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsOS(), kwargs...) where T

compute the Critical Line Segments by OSQP (Operator Splitting Quadratic Program) numerical solver, at LMEP. Save the CL to aCL if done


"""
function OpSpCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsOS(), kwargs...) where {T}
    #function OpSpCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    (; u, d, G, g, N, J) = PS
    (; tolS, muShft) = nS

    # using OSQP for HMFP and muShft are a disaster, LP or QP
    #= x = OpSpLP(PS; settings=settings)

    if x.info.status_val != 1
        error("Not able to find the expected return of HMFP (Highest Mean Frontier Portfolio)")
    end =#

    #=
    #mu = -x.info.obj_val
    mu = x.info.obj_val
    shft =  muShft
    if  mu < -1 || mu > 1
        shft *= abs(mu)
    end
    mu -= shft
    #y = OpSpQP(E, V, mu, u, d, G, g, A, b)
    y = OpSpQP(mu, PS; settings=settings)
    if y.info.status_val != 1   #SOLVED
        error("Not able to find a muShft to the HMFP (Highest Mean Frontier Portfolio)")
    end =#

    y = OpSpQP(PS; settings=settings)   #relative less Accuracy, but still very good
    if y.info.status_val != 1   #SOLVED
        error("Not able to find the LMEP (Lowest Mean Efficient Portfolio)")
    end

    Y = y.x
    S = fill(IN, N + J)
    Sv = @view S[1:N]
    Sv[abs.(Y - d).<tolS] .= DN
    Sv[abs.(Y - u).<tolS] .= UP

    for m = 1:J
        S[N+m] = abs(g[m] - G[m, :]' * Y) < tolS ? EO : OE
    end
    #return mu, S
    computeCL!(aCL, S, PS, nS)
end

end