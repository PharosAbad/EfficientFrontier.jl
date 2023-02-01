"using OSQP (Operator Splitting Quadratic Program) numerical solver to do LP and QP"
module uOSQP

using LinearAlgebra, OSQP, SparseArrays, EfficientFrontier
export OpSpQP, OpSpLP, OpSpCL!

function SettingsOS(; kwargs...)
    #= sD = Dict{Symbol,Any}()
    if !isempty(kwargs)
        for (key, value) in kwargs
            sD[key] = value
        end
    end
    OSQP.Settings(sD) =#
    if isempty(kwargs)  #default setting
        return SettingsOS(; eps_abs=2^-26, eps_rel=2^-26, eps_prim_inf=2^-26, eps_dual_inf=2^-26, verbose=false, polish=true, max_iter=771000)
    end
    kwargs
end


#=
# LVEP (Lowest Variance Efficient Portfolio) == Global Minimum Variance Portfolio (GMVP)
function OpSpQP(PS::Problem{T}; settings=SettingsOS()) where {T}
    #function OpSpQP(PS::Problem{T}; settings=OSQP.Settings()) where {T}
    (; V, u, d, G, g, A, b, N, J) = PS
    P = sparse(V)
    q = zeros(T, N)
    Ao = sparse([A; G; Matrix{T}(I, N, N)])
    uo = [b; g; u]
    lo = [b; fill(-Inf, J); d]
    model = OSQP.Model()
    OSQP.setup!(model; P=P, q=q, A=Ao, l=lo, u=uo, settings...) # `settings...` to iterate pairs in `settings`
    OSQP.solve!(model)
end
=#

function OpSpQP(PS::Problem{T}; settings=SettingsOS(), L::T=0.0) where {T}
    #function OpSpQP(PS::Problem{T}; settings=OSQP.Settings()) where {T}
    (; E, V, u, d, G, g, A, b, N, J) = PS
    P = sparse(V)
    #q = zeros(T, N)
    Ao = sparse([A; G; Matrix{T}(I, N, N)])
    uo = [b; g; u]
    lo = [b; fill(-Inf, J); d]
    model = OSQP.Model()
    #OSQP.setup!(model; P=P, q=q, A=Ao, l=lo, u=uo, settings...) # `settings...` to iterate pairs in `settings`
    #OSQP.solve!(model)
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
    return OSQP.solve!(model)    
end


#find the highest mean: -x.info.obj_val
function OpSpLP(PS::Problem{T}; settings=SettingsOS()) where {T}
#function OpSpLP(E::Vector{T}, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, A::Matrix{T}, b::Vector{T}) where {T}
    (; E, u, d, G, g, A, b, N, J) = PS
    P = spzeros(T, N, N)
    q = -E
    Ao = sparse([A; G; Matrix{T}(I, N, N)])
    uo = [b; g; u]
    lo = [b; fill(-Inf, J); d]
    model = OSQP.Model()
    OSQP.setup!(model; P=P, q=q, A=Ao, l=lo, u=uo, settings...)
    OSQP.solve!(model)
end

function OpSpQP(PS::Problem{T}, mu::T; settings=SettingsOS()) where {T}
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


function OpSpCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsOS(), kwargs...) where {T}
    #function OpSpCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    (; u, d, G, g, N, J) = PS
    (; tolS, muShft) = nS

    # using OSQP for HMFP and muShft are a disaster, LP or QP
    #= x = OpSpLP(PS; settings=settings)

    if x.info.status_val != 1
        error("Not able to find the expected return of HMFP (Highest Mean Frontier Portfolio)")
    end =#

    #= mu = -x.info.obj_val
    shft =  muShft
    if  mu < -1 || mu > 1
        shft *= abs(mu)
    end
    mu -= shft
    #y = OpSpQP(E, V, mu, u, d, G, g, A, b)
    y = OpSpQP(PS, mu; settings=settings)
    if y.info.status_val != 1   #SOLVED
        error("Not able to find a muShft to the HMFP (Highest Mean Frontier Portfolio)")
    end =#

    y = OpSpQP(PS; settings=settings)   #relative less Accuracy, but still very good
    if y.info.status_val != 1   #SOLVED
        error("Not able to find the LVEP (Lowest Variance Efficient Portfolio)")
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