"Status Switching Method for QP"
module SSQP

using LinearAlgebra

using EfficientFrontier: EfficientFrontier, Problem, Status, Event, IN, DN, UP, OE, EO
import EfficientFrontier: SettingsLP

using EfficientFrontier.Simplex: cDantzigLP, maxImprvLP, Simplex

export solveQP, QP

"""

        QP(V::Matrix{T}; q, u, d, G, g, A, b) where T
        QP(P::Problem{T}, L::T) where T
        QP(mu::T, P::Problem{T}) where T

Setup a quadratic programming model:

```math
    min   (1/2)z′Vz+q′z
    s.t.   Az=b ∈ R^{M}
           Gz≤g ∈ R^{J}
           d≤z≤u ∈ R^{N}
```

some variable may be free, say -Inf < zi < +Inf. No equalities if M=0. Default values: q = 0, u = +∞, d = 0, G = [], g = [], A = ones(1,N), b = [1]

See also [`Problem`](@ref), [`solveQP`](@ref)

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

    #check feasibility and redundancy of Ax=b
    rb = rank([A b])
    @assert rb == rank(A) "infeasible: Ax=b"
    @assert M == rb "redundant rows in Ax=b"       #full row rank

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
    #Am = [E'; A]
    #bm = [mu; b]
    Am = [A; E']
    bm = [b; mu]
    return QP(V, Am, G, q, bm, g, d, u, N, M, J)
end


function SettingsLP(Q::QP{T}; kwargs...) where {T}
    Simplex.Settings{T}(; kwargs...)
end


"""

        Settings(Q::QP{T}; kwargs...)        The default Settings to given quadratic programming Q
        Settings(; kwargs...)       The default Settings is set by Float64 type
        Settings{T<:AbstractFloat}(; kwargs...)

kwargs are from the fields of Settings{T<:AbstractFloat} for Float64 and BigFloat

            maxIter::Int        #7777
            tol::T              #2^-26 ≈ 1.5e-8   general scalar
            tolN::T             #2^-26 ≈ 1.5e-8   for norms
            tolG::T             #2^-33 ≈ 1.2e-10  for Greeks (beta and gamma) "portfolio weights"

See also [`QP`](@ref), [`solveQP`](@ref)
"""
struct Settings{T<:AbstractFloat}
    maxIter::Int    #7777
    tol::T          #2^-26
    tolN::T         #2^-26
    tolG::T         #2^-33 for Greeks (beta and gamma)
end

Settings(; kwargs...) = Settings{Float64}(; kwargs...)

function Settings{Float64}(; maxIter=7777,
    tol=2^-26,
    tolN=2^-26,
    tolG=2^-33)
    Settings{Float64}(maxIter, tol, tolN, tolG)
end

function Settings{BigFloat}(; maxIter=7777,
    tol=BigFloat(2)^-76,
    tolN=BigFloat(2)^-76,
    tolG=BigFloat(2)^-87)
    Settings{BigFloat}(maxIter, tol, tolN, tolG)
end

function Settings(Q::QP{T}; kwargs...) where {T}
    Settings{T}(; kwargs...)
end



function getRows(A::Matrix{T}, tol=sqrt(eps(T))) where {T}
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
        if norm(v, Inf) <= tol
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
        if norm(v, Inf) > tol && norm(v - H' * (H' \ v), Inf) > tol
            R[r] = true
            H = @view A[R, :]
        end
    end
    return findall(R)
end

function freeK!(S, z, V, q, N, tol)  #for K=0
    #modify: S
    p = V * z + q
    #= if norm(p, Inf) < tol
        return 1
    end =#
    S0 = copy(S)

    t = true   #hit optimal
    @inbounds for k in 1:N
        if (p[k] >= -tol && S[k] == UP) || (p[k] <= tol && S[k] == DN)
            #if (p[k] > tol && S[k] == UP) || (p[k] < -tol && S[k] == DN)
            S[k] = IN
            t = false
        end
    end
    if t
        return 1
    else
        ip = findall(S .== IN)
        if length(ip) > 0 && norm(p[ip], Inf) <= tol  #all movable are optimal
            S[ip] = S0[ip]  #restore the status
            return 1
        end
        return -1
    end
    #= if !t && norm(p[S .== IN], Inf) <= tol  #all movable are optimal
        return 1
    end =#
    #return t ? 1 : -1

end

function aStep!(p, z::Vector{T}, S, F, Og, alpha, G, g, d, u, fu, fd, N, J, tol) where {T}
    #compute step
    Lo = Vector{Event{T}}(undef, 0)
    ik = findall(F)
    @inbounds for k in eachindex(alpha)
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

    @inbounds if J > 0
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

    @inbounds if L1 < 1.0
        #adding blocking constraints
        z[F] .+= L1 * p
        #multi blocking
        for i in 1:lastindex(Lo)
            Lt = Lo[i]
            if Lt.L - L1 > tol
                break
            end
            k = Lt.id
            To = Lt.To
            if To == EO
                k += N
            end
            S[k] = To
            if k <= N
                z[k] = To == DN ? d[k] : u[k]
            end
        end
        #= ik = findall(F)
        for k in ik
            if abs(z[k] - d[k]) < tol
                z[k] = d[k]
                S[k] = DN
            elseif abs(z[k] - u[k]) < tol
                z[k] = u[k]
                S[k] = UP
            end
        end =#
        return -1
    else
        # if step size L1 == 1.0, some z_i hit bounds
        z[F] = alpha
        return 1
    end

end

function KKTchk!(S, F, B, Eg, gamma, alphaL, AE, GE, idAE, ra, N, M, tolG::T) where {T}
    ib = findall(B)
    Li = Vector{Event{T}}(undef, 0)
    @inbounds for k in eachindex(gamma)
        j = ib[k]
        t = gamma[k]
        if S[j] == UP && t > tolG
            push!(Li, Event{T}(UP, IN, j, -t))
        elseif S[j] == DN && t < -tolG
            push!(Li, Event{T}(DN, IN, j, t))
        end
    end

    JE = size(GE, 1)
    @inbounds if JE > 0
        iE = zeros(Int, JE)
        iR = findall(ra .> M)
        iE[idAE[ra[iR]]] = iR
        Lda = zeros(T, JE)
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
            if t < -tolG
                push!(Li, Event{T}(EO, OE, ib[k], t))
            end
        end
    end

    nL = length(Li)
    @inbounds if nL > 0   #entering one only, the tighest one
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

"""

        solveQP(Q::QP{T}; settings, settingsLP) where T
        solveQP(V::Matrix{T}, q::Vector{T}, A, b, G, g, d, u; settings, settingsLP) where T
        solveQP(Q::QP{T}, S::Vector{Status}, x0; settings) where T

for quadratic programming problems: the initial feasible point (S, x0) if given,  can be obtained from SimplexLP

```math
    min   (1/2)z′Vz+q′z
    s.t.   Az=b ∈ R^{M}
           Gz≤g ∈ R^{J}
           d≤z≤u ∈ R^{N}
```

Outputs

    z               : solution,  N x 1 vector
    S               : Vector{Status}, (N+J)x1
    status          : > 0 if successful (=iter_count), 0 if infeasibility detected, < 0 if not converged (=-iter_count)

See also [`QP`](@ref), [`EfficientFrontier.Simplex.SimplexLP`](@ref), [`EfficientFrontier.SSQP.Settings`](@ref), [`EfficientFrontier.Simplex.Settings`](@ref), [`initQP`](@ref), [`initSSQP`](@ref)
"""
function solveQP(V::Matrix{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T}, G::Matrix{T}, g::Vector{T},
    d::Vector{T}, u::Vector{T}; settings=Settings{T}(), settingsLP=Simplex.Settings{T}()) where {T}

    #N::Int32 = length(q)
    N = length(q)
    M = length(b)
    J = length(g)
    Q = QP(V, A, G, q, b, g, d, u, N, M, J)
    solveQP(Q; settings=settings, settingsLP=settingsLP)
end

function solveQP(Q::QP{T}; settings=Settings(Q), settingsLP=SettingsLP(Q)) where {T}
    S, x0 = initSSQP(Q, settingsLP)
    solveQP(Q, S, x0; settings=settings)
end


function solveQP(Q::QP{T}, S, x0; settings=Settings(Q)) where {T}
    (; V, A, G, q, b, g, d, u, N, M, J) = Q
    #(; maxIter, tol, tolN, tolG) = settings
    (; maxIter, tol, tolG) = settings

    fu = u .< Inf   #finite upper bound
    fd = d .> -Inf   #finite lower bound

    Sz = @view S[1:N]
    Se = @view S[(N.+(1:J))]
    z = copy(x0)

    iter = 0
    @inbounds while true
        #while true
        iter += 1
        if iter > maxIter
            return z, S, -iter
        end

        F = (Sz .== IN)
        K = sum(F)
        if K == 0
            status = freeK!(S, z, V, q, N, tol)
            if status > 0
                return z, S, iter
            else
                continue
            end
        end

        B = .!F
        Eg = (Se .== EO)
        Og = (Se .== OE)
        GE = @view G[Eg, :]
        AE = vcat(A[:, F], GE[:, F])
        idAE = vcat(axes(A, 1), axes(GE, 1)) # id of each constraint
        zB = z[B]
        AB = vcat(A[:, B], GE[:, B])
        bE = vcat(b, g[Eg]) - AB * zB

        ra = getRows(AE, tol)
        W = length(ra)
        if W < length(bE)
            rb = getRows([AE bE], tol)
            if W != length(rb)
                return z, S, 0    #infeasible
            end
            AE = AE[ra, :]
            bE = bE[ra]
            AB = AB[ra, :]
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
        if norm(p, Inf) > tolG  #= IMPORTANT: in theory, norm(p, Inf)=0  == norm(p)=0, but in numerical computation, NO!  p = (e/2)*1,  norm(p) = (e/2)*sqrt(N) > e if N>4 =#
            #if norm(p) > tolN
            status = aStep!(p, z, S, F, Og, alpha, G, g, d, u, fu, fd, N, J, tol)
            if status < 0
                continue
            end
        #else
            #z[F] = alpha
        end

        #= if iter >= 61
            display(iter)
        end =#

        #direction p = 0
        #α_{λ}=-C(T′c+b_{E})
        alphaL = -(TC' * c + C * bE)
        gamma = VBF * alpha + V[B, B] * zB + q[B] + AB' * alphaL
        #z[F] = alpha
        status = KKTchk!(S, F, B, Eg, gamma, alphaL, AE, GE, idAE, ra, N, M, tolG)
        if status > 0
            ik = findall(F)
            for k in ik #check fake IN
                if abs(z[k] - d[k]) < tol
                    z[k] = d[k]
                    S[k] = DN
                elseif abs(z[k] - u[k]) < tol
                    z[k] = u[k]
                    S[k] = UP
                end
            end

            #should we compute the final analytical z[I]?
            #z1 = copy(z)
            #z1[F] = alphaCal(F, z, Se, V, A, G, q, b, g, tol)
            #z[F] = (z[F]+z1[F])/2

            return z, S, iter
        end
    end
end

#=
function alphaCal(F, z, Se, V, A, G, q, b, g, tol)
    B = .!F
    Eg = (Se .== EO)
    GE = @view G[Eg, :]
    AE = vcat(A[:, F], GE[:, F])

    zB = z[B]
    AB = vcat(A[:, B], GE[:, B])
    bE = vcat(b, g[Eg]) - AB * zB

    ra = getRows(AE, tol)
    W = length(ra)
    if W < length(bE)
        rb = getRows([AE bE], tol)
        if W != length(rb)
            return z, S, 0    #infeasible
        end
        AE = AE[ra, :]
        bE = bE[ra]
        AB = AB[ra, :]
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
    #p = alpha - z[F]
    return alpha
end
=#


"""
        S, x0 = initSSQP(Q::QP{T}, settingsLP)


do not handle free variables such that -∞ < x < +∞. and d shoue be finite. OK for EfficientFrontier
"""
function initSSQP(Q::QP{T}, settingsLP) where {T}
    (; A, G, b, g, d, u, N, M, J) = Q
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
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #灯塔的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    _iH, x, _invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    f = sum(x[Ns+1:end])
    if f > tol
        error("feasible region is empty")
    end

    x0 = x[1:N]
    S = S[1:N+J]
    for k in N+1:N+J
        S[k] = S[k] == IN ? OE : EO
    end
    return S, x0

end

"""
        S, x0 = initQP(Q::QP{T}, settingsLP)

performing Phase-I Simplex on the polyhedron {Az=b, Gz≤g, d≤z≤u} to find an initial feasible point
allowing free variables such that -∞ < z < +∞
"""
function initQP(Q::QP{T}, settingsLP) where {T}
    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    (; A, G, b, g, d, u, N, M, J) = Q
    (; tol, rule) = settingsLP

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    end

    #convert free variable: -∞ < x < +∞
    fu = u .== Inf   #no upper bound
    fd = d .== -Inf   #no lower bound
    fv = fu .&& fd  #free variable
    iv = findall(fv)
    n = length(iv)
    id = findall(fd .&& .!fv)   # (-∞, u]

    #add slack variables for Gz<=g , and 2nd part of free variables
    M0 = M + J
    N0 = N + J + n
    A0 = [A zeros(T, M, J) -A[:, iv]
        G Matrix{T}(I, J, J) -G[:, iv]]
    b0 = [b; g]
    d0 = [d; zeros(T, J + n)]
    u0 = [u; fill(Inf, J + n)]


    #(-∞, +∞)  -> two copy of [0, +∞)
    d0[iv] .= 0     #the 1st part of free variable
    #u0[iv] .= Inf

    #(-∞, u]  -> [-u, +∞)
    d0[id] .= -u0[id]
    u0[id] .= Inf
    A0[:, id] .= -A0[:, id]

    N1 = M0 + N0
    S = fill(DN, N1)
    B = collect(N0 .+ (1:M0))
    S[B] .= IN

    invB = Matrix{T}(I, M0, M0)
    q = A0 * d0
    for j in 1:M0
        invB[j, j] = b0[j] >= q[j] ? one(T) : -one(T)
    end
    q = abs.(q - b0)
    c1 = [zeros(T, N0); fill(one(T), M0)]   #灯塔的　模型　是　min
    A1 = [A0 invB]
    b1 = b0
    d1 = [d0; zeros(T, M0)]
    u1 = [u0; fill(Inf, M0)]


    #f, x, _q, _invB, _iH = EfficientFrontier.Simplex.cDantzigLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    _iH, x0, _invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    f = sum(x[N0+1:end])
    if f > tol
        #if abs(f) > tol
        #display(f)
        error("feasible region is empty")
    end

    x = x0[1:N]
    S = S[1:N+J]

    for k in N+1:N+J    #inequalities
        S[k] = S[k] == IN ? OE : EO
    end

    if n > 0    #free variable
        x[iv] .-= x0[N+J+1:N+J+n]
        S[iv] .= IN
    end

    m = length(id)
    if m > 0   # flip u d
        x[id] = -x[id]
        for k in 1:m
            #S[k] = S[k] == DN ? UP : DN
            if S[k] == DN
                S[k] == UP
            end
        end
    end
    return S, x
end

end

