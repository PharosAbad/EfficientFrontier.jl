"Status Switching Method for QP"
module SSQP

using LinearAlgebra

using EfficientFrontier: EfficientFrontier, Problem, Status, Event, IN, DN, UP, OE, EO
import EfficientFrontier: SettingsLP

using EfficientFrontier.Simplex: cDantzigLP, maxImprvLP, Simplex

export solveQP, QP

"""


min   (1/2)z′Vz+q′z
s.t.   Az=b ∈ R^{M}
       Gz≤g ∈ R^{J}
       d≤z≤u ∈ R^{N}

some variable may be free, say -Inf < zi < +Inf
"""
struct QP{T<:AbstractFloat}    #standard QP, or structure of QP
    V::Matrix{T}
    A::Matrix{T}
    G::Matrix{T}
    q::Vector{T}
    b::Vector{T}
    g::Vector{T}
    d::Vector{T}
    u::Vector{T}
    N::Int
    M::Int
    J::Int
end

function QP(V::Matrix{T}; N=size(V, 1), #N=convert(Int32, size(V, 1)),
    q=zeros(N),
    u=fill(Inf, N),
    d=zeros(N),
    G=ones(0, N),
    g=ones(0),
    A=ones(1, N),
    b=ones(1)) where {T}

    #N::Int32 = length(q)
    #M::Int32 = length(b)
    #J::Int32 = length(g)
    M = length(b)
    J = length(g)

    (N, N) == size(V) || throw(DimensionMismatch("incompatible dimension: V"))
    V = convert(Matrix{T}, (V + V') / 2)   #make sure symmetric
    @assert det(V) >= 0 "variance matrix has negative determinant"
    (M, N) == size(A) || throw(DimensionMismatch("incompatible dimension: A"))
    (J, N) == size(G) || throw(DimensionMismatch("incompatible dimension: G"))
    N == size(q, 1) || throw(DimensionMismatch("incompatible dimension: q"))
    #M == size(b, 1) || throw(DimensionMismatch("incompatible dimension: b"))
    #J == size(g, 1) || throw(DimensionMismatch("incompatible dimension: g"))
    N == size(d, 1) || throw(DimensionMismatch("incompatible dimension: d"))
    N == size(u, 1) || throw(DimensionMismatch("incompatible dimension: u"))

    iu = u .< d
    if sum(iu) > 0
        @warn "swap the elements where u < d, to make sure u > d"
        t = u[iu]
        u[iu] .= d[iu]
        d[iu] .= t
    end
    #to do: J+ (num of finite d u) > 0  , when LP is introduce

    @assert !any(d .== u) "downside bound == upper bound detected"
    @assert J > 0 || any(isfinite.(d)) || any(isfinite.(u)) "no any inequalities or bounds"

    QP{T}(V, convert(Matrix{T}, copy(A)),
        convert(Matrix{T}, copy(G)),
        convert(Vector{T}, copy(q)),
        convert(Vector{T}, copy(vec(b))),
        convert(Vector{T}, copy(vec(g))),
        convert(Vector{T}, copy(d)),
        convert(Vector{T}, copy(u)), N, M, J)
end


function QP(P::Problem{T}, L::T=0.0) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = P
    q = -L * E
    return QP(V, A, G, q, b, g, d, u, N, M, J)
end

function QP(mu::T, P::Problem{T}) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = P
    q = zeros(T, N)
    M += 1
    Am = [E'; A]
    bm = [mu; b]
    return QP(V, Am, G, q, bm, g, d, u, N, M, J)
end


function SettingsLP(Q::QP{T}; kwargs...) where {T}
    Simplex.Settings{T}(; kwargs...)
end

struct Settings{T<:AbstractFloat}
    maxIter::Int64  #7777
    tol::T          #2^-26
    tolNorm::T      #2^-26
end

Settings(; kwargs...) = Settings{Float64}(; kwargs...)

function Settings{Float64}(; maxIter=7777,
    tol=2^-26,
    tolNorm=2^-26)
    Settings{Float64}(maxIter, tol, tolNorm)
end

function Settings{BigFloat}(; maxIter=7777,
    tol=BigFloat(2)^-76,
    tolNorm=BigFloat(2)^-76)
    Settings{BigFloat}(maxIter, tol, tolNorm)
end

function Settings(Q::QP{T}; kwargs...) where {T}
    Settings{T}(; kwargs...)
end



function getRows(A::Matrix{T}, tolNorm=sqrt(eps(T))) where {T}
    #indicate the non-redundant rows, the begaining few rows can be zeros (redundant)
    M, N = size(A)
    if N == 0
        @warn "zero columns" size(A)
        return collect(axes(A, 1))
    end
    R = falses(M)
    if M == 0
        return findall(R)
    end

    r1 = M + 1
    #find the 1st non-zero row
    for r in 1:M
        v = @view A[r, :]
        if norm(v) <= tolNorm
            continue
        else
            R[r] = true
            r1 = r + 1
            break
        end
    end

    #rows after the 1st non-zero row
    H = @view A[R, :]
    for r in r1:M
        #v = @view A[r:r, :]
        v = @view A[r, :]
        if norm(v) > tolNorm && norm(v - H' * (H' \ v)) > tolNorm
            R[r] = true
            H = @view A[R, :]
        end
    end
    return findall(R)
end


function solveQP(V::Matrix{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T}, G::Matrix{T}, g::Vector{T},
    d::Vector{T}, u::Vector{T}; settings=Settings{T}(), settingsLP=Simplex.Settings{T}()) where {T}

    #N::Int32 = length(q)
    N = length(q)
    M = length(b)
    J = length(g)
    Q = QP(V, A, G, q, b, g, d, u, N, M, J)
    solveQP(Q; settings=settings, settingsLP=settingsLP)
end


function aStep!(p, z::Vector{T}, S, F, Og, alpha, G, g, d, u, fu, fd, N, J, tol) where {T}
    #compute step
    Lo = Vector{Event{T}}(undef, 0)
    ik = findall(F)
    for k in eachindex(alpha)
        j = ik[k]
        t = p[k]
        h = z[j]
        dL = (d[j] - h) / t
        uL = (u[j] - h) / t
        if t > tol && fu[j]
            push!(Lo, Event{T}(IN, UP, j, uL))
        elseif t < -tol && fd[j]
            push!(Lo, Event{T}(IN, DN, j, dL))
        end
    end

    if J > 0
        zo = g[Og] - G[Og, :] * z
        po = G[Og, F] * p
        ik = findall(Og)
        for k in eachindex(zo)
            j = ik[k]
            t = po[k]
            if t > tol
                push!(Lo, Event{T}(OE, EO, j, zo[k] / t))
            end
        end
    end

    L1::T = 1.0
    nL = length(Lo)
    if nL > 0
        sort!(Lo, by=x -> x.L)
        L1 = Lo[1].L
    end

    if L1 < 1.0
        #adding blocking constraints
        z[F] .+= L1 * p
        #multi blocking
        for i in 1:lastindex(Lo)
            Lt = Lo[i]
            k = Lt.id
            To = Lt.To
            if Lt.L - L1 < tol
                if To == EO
                    k += N
                end
                S[k] = To
                if k <= N
                    z[k] = To == DN ? d[k] : u[k]
                end
            end
        end
        return -1
    else
        # if step size L1 == 1.0, some z_i hit bounds
        z[F] = alpha
        return 1
    end

end

function KKTchk!(S, B, Eg, gamma, alphaL, GE, idAE, ra, M, tol::T) where {T}
    ib = findall(B)
    Li = Vector{Event{T}}(undef, 0)
    for k in eachindex(gamma)
        j = ib[k]
        t = gamma[k]
        if S[j] == UP && t > tol
            push!(Li, Event{T}(UP, IN, j, -t))
        elseif S[j] == DN && t < -tol
            push!(Li, Event{T}(DN, IN, j, t))
        end
    end

    JE = size(GE, 1)
    if JE > 0
        iE = zeros(Int, JE)
        iR = findall(ra .> M)
        iE[idAE[ra[iR]]] = iR
        Lda = zeros(JE)
        for j in 1:JE
            k = iE[j]
            if k == 0
                x = AE' \ GE[j, F]
                Lda[j] = alphaL' * x
            else
                Lda[j] = alphaL[k]
            end
        end

        ib = findall(Eg)
        for k in 1:JE
            t = Lda[k]
            if t < -tol
                push!(Li, Event{T}(EO, OE, ib[k], t))
            end
        end
    end

    nL = length(Li)
    if nL > 0   #entering one only, the tighest one
        sort!(Li, by=x -> x.L)
        Lt = Li[1]
        k = Lt.id
        To = Lt.To
        if To == OE
            k += N
        end
        S[k] = To
        return -1
    else
        return 1
    end
end

function solveQP(Q::QP{T}; settings=Settings(Q), settingsLP=SettingsLP(Q), x0=nothing, S=nothing) where {T}
    #function solveQP(Q::QP{T}; settings=Settings(Q), settingsLP=SettingsLP(Q)) where {T}
    #function solveQP(Q::QP{T}; settings=Settings{T}(), settingsLP=SettingsLP(Q)) where {T}
    (; V, A, G, q, b, g, d, u, N, M, J) = Q
    (; maxIter, tol, tolNorm) = settings

    if isnothing(S)
        x0, S = initSSQP(Q, settingsLP)
    end

    fu = u .< Inf   #finite upper bound
    fd = d .> -Inf   #finite lower bound

    Sz = @view S[1:N]
    Se = @view S[(N.+(1:J))]
    z = copy(x0)


    iter = 0
    #@inbounds     while true
    while true
        iter += 1
        if iter >= maxIter
            return z, S, -iter
        end

        #=
        if iter == 6 #mod(iter, 10) == 0
            display((iter, S))
        end =#

        F = (Sz .== IN)
        B = .!F
        Eg = (Se .== EO)
        Og = (Se .== OE)
        GE = @view G[Eg, :]
        AE = vcat(A[:, F], GE[:, F])
        idAE = vcat(axes(A, 1), axes(GE, 1)) # id of each constraint
        zB = z[B]
        AB = vcat(A[:, B], GE[:, B])
        bE = vcat(b, g[Eg]) - AB * zB

        ra = getRows(AE, tolNorm)
        W = length(ra)
        if W < length(bE)
            rb = getRows([AE bE], tolNorm)
            if W != length(rb)
                return z, S, 0    #infeasible
            end
            AE = AE[ra, :]
            bE = bE[ra]
            AB = AB[ra, :]
        end

        K = sum(F)          #to do: K == 0
        #display((S, iter, K))
        if K == 0
            @warn "Hit a degenerate point, moving variables are all on the bounds"
            return z, S, -iter
        end


        iV = inv(cholesky(V[F, F]))
        VBF = V[B, F]
        c = VBF' * zB + q[F]
        mT = iV * AE'   #T=V_{I}⁻¹A_{I}′
        C = AE * mT
        C = (C + C') / 2
        C = inv(cholesky(C))
        TC = mT * C
        VQ = iV - mT * TC'    #Q=I-A_{I}′CT′   V_{I}⁻¹Q=V_{I}⁻¹-TCT′
        alpha = TC * bE - VQ * c    #α=TCb_{E}-V_{I}⁻¹Qc
        p = alpha - z[F]

        #direction p ≠ 0
        if norm(p) > tolNorm
            status = aStep!(p, z, S, F, Og, alpha, G, g, d, u, fu, fd, N, J, tol)
            if status < 0
                continue
            end
        end

        #direction p = 0
        #α_{λ}=-C(T′c+b_{E})
        alphaL = -(TC' * c + C * bE)
        gamma = VBF * alpha + V[B, B] * zB + q[B] + AB' * alphaL
        #z[F] = alpha
        status = KKTchk!(S, B, Eg, gamma, alphaL, GE, idAE, ra, M, tol)
        if status > 0
            return z, S, iter
        end
    end
end


function initSSQP(Q::QP{T}, settingsLP) where {T}
    (; A, G, q, b, g, d, u, N, M, J) = Q
    (; tol, rule) = settingsLP

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    end

    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    Ms = M + J  #convert Gz<=g to equality
    Ns = N + J
    N1 = Ms + Ns
    S = fill(DN, N1)
    B = collect(Ns .+ (1:Ms))
    S[B] .= IN

    As = [A zeros(T, M, J)
        G Matrix{T}(I, J, J)]
    invB = Matrix{T}(I, Ms, Ms)
    bs = [b; g]
    ds = [d; zeros(T, J)]
    us = [u; fill(Inf, J)]

    q = As * ds
    for j in 1:Ms
        invB[j, j] = bs[j] >= q[j] ? one(T) : -one(T)
    end
    #q = abs.(As * ds - bs)
    q = abs.(q - bs)
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #我的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    f, x, _q, _invB, _iH = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    if f > tol
        error("feasible region is empty")
    end

    x0 = x[1:N]
    S = S[1:N+J]
    for k in N+1:N+J
        S[k] = S[k] == IN ? OE : EO
    end
    return x0, S

end

end

