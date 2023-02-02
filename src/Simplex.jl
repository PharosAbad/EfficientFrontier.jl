"Simplex algorithm for EfficientFrontier"
module Simplex
using LinearAlgebra, Combinatorics
using EfficientFrontier: EfficientFrontier, Problem, Status, Event, sCL, IN, DN, UP, OE, EO, computeCL!
import EfficientFrontier: Settings as SettingsEF    #if by `using`, the compiler will complain conflict when import Simplex.Settings in to EfficientFrontier as SettingsLP
export SimplexLP #,SimplexCL!

#=
REMARK: writing the next basis as a product of the current basis times an easily invertible matrix can be extended over several iterations. We do not adopt this for accuracy
    If this method is adopt, how often should one recompute an inverse of the current basis?
=#

"""

        Settings(P::Problem; kwargs...)        The default Settings to given Problem
        Settings(; kwargs...)       The default Settings is set by Float64 type
        Settings{T<:AbstractFloat}(; kwargs...)

kwargs are from the fields of Settings{T<:AbstractFloat} for Float64 and BigFloat

            tol::T          #general scalar
            rule::Symbol    #rule for Simplex {:Dantzig, :maxImprovement}

"""
struct Settings{T<:AbstractFloat}
    tol::T         #general scalar
    rule::Symbol    #rule for Simplex {:Dantzig, :maxImprovement}
end

Settings(; kwargs...) = Settings{Float64}(; kwargs...)
function Settings{Float64}(; tol=2^-26, rule=:Dantzig)
    Settings{Float64}(tol, rule)
end

function Settings{BigFloat}(; tol=BigFloat(2)^-76,
    rule=:Dantzig)
    Settings{BigFloat}(tol, rule)
end

function Settings(P::Problem{T}; kwargs...) where {T}
    Settings{T}(; kwargs...)
end


"""
        cDantzigLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)

using  combined Dantzig's pivot rule to solve LP (combine the largest-distance rule with Dantzig's pivot rule, and switch to Bland's rule if iters > N)
invB: inv(B)
q:    x[B]

```math
        min   z=c′x
        s.t.  Ax=b
              d≤x≤u
```
"""
function cDantzigLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)
    #B is always sorted. S in caller is change, compute invB and q=xB each step, switch Dantzig to Bland rule when iter > N 

    T = typeof(c[1])
    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)
    ip = zeros(Int64, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    cA = zeros(T, N)    #norm of each column of A
    for k in 1:N
        cA[k] = norm(A[:, k])
    end

    x = zeros(T, N)
    x[S.==UP] .= u[S.==UP]
    x[S.==DN] .= d[S.==DN]

    Y = invB * A[:, F]
    h = c[F] - Y' * c[B]
    ih = S[F] .== DN
    h[ih] .= -h[ih]
    ih = h .> tol
    hp = h[ih]
    iH = findall(F)[ih]
    nH = length(iH)
    Bland = false
    loop = 0
    while nH > 0
        #Dantzig rule or Bland rule
        #k0 = Bland ? 1 : argmax(hp)
        k0 = Bland ? 1 : argmax(hp ./ cA[iH]) #Largest-Distance Rule
        k = iH[k0]
        p = invB * A[:, k]
        kd = S[k] == DN
        m = 0
        if kd
            for j in 1:M
                i = B[j]
                if p[j] > tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                elseif p[j] < -tol && u[i] < Inf
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                end
            end

            #= if m == 0
                error("infeasible or unbounded")
            end
            g0 = u[k] - d[k]
            (gl, l) = findmin(gt[1:m])
            Sl = Sb[l]
            l = ip[l]
            if gl > g0
                S[k] = UP
                x[k] = u[k]
                q -= g0 * p
                l = -k
            end =#
            g0 = u[k] - d[k]
            if m == 0
                if isfinite(g0) #DN -> UP
                    S[k] = UP
                    x[k] = u[k]
                    q -= g0 * p
                    l = -k
                else
                    error("infeasible or unbounded")
                end
            else
                (gl, l) = findmin(gt[1:m])
                Sl = Sb[l]
                l = ip[l]
                if gl > g0
                    S[k] = UP
                    x[k] = u[k]
                    q -= g0 * p
                    l = -k
                end
            end

        else
            for j in 1:M
                i = B[j]
                if p[j] > tol && u[i] < Inf
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                elseif p[j] < -tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                end
            end

            #= if m == 0
                error("infeasible or unbounded")
            end
            g0 = d[k] - u[k]
            (gl, l) = findmax(gt[1:m])
            Sl = Sb[l]
            l = ip[l]
            if gl < g0
                S[k] = DN
                x[k] = d[k]
                q -= g0 * p
                l = -k
            end =#

            g0 = d[k] - u[k]
            if m == 0
                if isfinite(g0) #UP -> DN
                    S[k] = DN
                    x[k] = d[k]
                    q -= g0 * p
                    l = -k
                else
                    error("infeasible or unbounded")
                end
            else
                (gl, l) = findmax(gt[1:m])
                Sl = Sb[l]
                l = ip[l]
                if gl < g0
                    S[k] = DN
                    x[k] = d[k]
                    q -= g0 * p
                    l = -k
                end
            end

        end
        loop += 1
        if loop > N
            Bland = true    #Dantzig rule switch to Bland rule
        end
        if l > 0
            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            B = sort(B)
            invB = inv(A[:, B])

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            q = invB * b - Y * x[F]
        end
        h = c[F] - Y' * c[B]
        ih = S[F] .== DN
        h[ih] .= -h[ih]
        ih = h .> tol
        hp = h[ih]
        iH = findall(F)[ih]
        nH = length(iH)
    end

    #@. ih = abs(h) < tol   # h==0
    ih = abs.(h) .< tol   # h==0
    iH = findall(F)[ih]
    x[B] = q
    #return q, B, invB, iH, c' * x
    return c' * x, x, q, B, invB, iH
end


"""
        maxImprvLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)

using  max improvement pivot rule to solve LP (no cycle since it becomes Bland's rule if the improvement is 0)
invB: inv(B)
q:    x[B]

```math
        min   z=c′x
        s.t.  Ax=b
              d≤x≤u
```
"""
function maxImprvLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)
    #greatest improvement, B is always sorted. S in caller is change, compute invB and q=xB each step
    T = typeof(c[1])
    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)    #theta
    ip = zeros(Int64, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    x = zeros(T, N)
    x[S.==UP] .= u[S.==UP]
    x[S.==DN] .= d[S.==DN]


    #h = c[F] - (invB * mF)' * c[B]
    Y = invB * A[:, F]
    h = c[F] - Y' * c[B]
    vS = S[F]   #State of leaving to be, for each candidate k
    g = zeros(N - M)    #theta for each candidate k
    gi = zeros(Int64, N - M)    #min subscrip for theta for each candidate k
    vD = falses(N - M)  #DN or not, for each candidate k
    ud = u - d

    ih = S[F] .== DN
    h[ih] .= -h[ih]
    ih = h .> tol
    #hp = h[ih]
    iH = findall(F)[ih]
    nH = length(iH)
    while nH > 0
        P = @view Y[:, ih]
        for n in 1:nH
            k = iH[n]
            p = P[:, n]
            kd = S[k] == DN
            vD[n] = kd
            m = 0
            if kd
                for j in 1:M
                    i = B[j]
                    if p[j] > tol
                        m += 1
                        gt[m] = (q[j] - d[i]) / p[j]
                        ip[m] = j
                        Sb[m] = DN
                    elseif p[j] < -tol && u[i] < Inf
                        m += 1
                        gt[m] = (q[j] - u[i]) / p[j]
                        ip[m] = j
                        Sb[m] = UP
                    end
                end

                if m == 0
                    #error("infeasible or unbounded")
                    #g[n] = isfinite(u[k]) ? ud[k] : 0
                    g[n] = ud[k]
                    l = 1
                    ip[1] = 0   #isfinite(u[k]) to flip state, isinf(u[k]) unbounded
                else
                    (g[n], l) = findmin(gt[1:m])
                end
                #(g[n], l) = findmin(gt[1:m])
                if g[n] > ud[k]
                    g[n] = ud[k]
                end
            else
                for j in 1:M
                    i = B[j]
                    if p[j] > tol && u[i] < Inf
                        m += 1
                        gt[m] = (q[j] - u[i]) / p[j]
                        ip[m] = j
                        Sb[m] = UP
                    elseif p[j] < -tol
                        m += 1
                        gt[m] = (q[j] - d[i]) / p[j]
                        ip[m] = j
                        Sb[m] = DN
                    end
                end

                if m == 0
                    #error("infeasible or unbounded")
                    #g[n] = isfinite(ud[k]) ? ud[k] : 0
                    g[n] = -ud[k]
                    l = 1
                    ip[1] = 0   #isfinite(u[k]) to flip state, isinf(u[k]) unbounded
                else
                    (g[n], l) = findmax(gt[1:m])
                end
                #(g[n], l) = findmax(gt[1:m])
                if g[n] < -ud[k]
                    g[n] = -ud[k]
                end
            end
            gi[n] = ip[l]
            vS[n] = Sb[l]
        end
        k = getfield(findmax(abs.(h[ih] .* g[1:nH])), 2)
        l = gi[k]
        if l == 0 && isinf(u[k])
            error("infeasible or unbounded")
        end
        # if l == 0 && isfinite(u[k]),  then gl == ±ud[k], will handle later
        kd = vD[k]
        p = P[:, k]
        gl = g[k]
        Sl = vS[k]
        k = iH[k]   #entering index
        if kd
            if gl == ud[k]  # l>=0, maybe l==0
                S[k] = UP
                x[k] = u[k]
                q -= gl * p
                l = -k
            end
        else
            if gl == -ud[k]
                S[k] = DN
                x[k] = d[k]
                q -= gl * p
                l = -k
            end
        end
        #remove all l==0 case, now, if l<0, we have flip the sate of k
        if l > 0
            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            B = sort(B)
            invB = inv(A[:, B])

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            q = invB * b - Y * x[F]
        end

        h = c[F] - Y' * c[B]
        ih = S[F] .== DN
        h[ih] .= -h[ih]
        ih = h .> tol
        #hp = h[ih]
        iH = findall(F)[ih]
        nH = length(iH)
    end

    ih = abs.(h) .< tol   # h==0
    iH = findall(F)[ih]
    x[B] = q
    #return q, B, invB, iH, c' * x
    return c' * x, x, q, B, invB, iH
end



"""

        SimplexLP(PS::Problem; settings=Simplex.Settings(PS), min=true)

find the `Status` for assets by simplex method. If `min=false`, we maximize the objective function

See also [`Status`](@ref), [`Problem`](@ref), [`EfficientFrontier.Simplex.Settings`](@ref), [`EfficientFrontier.Simplex.cDantzigLP`](@ref), [`EfficientFrontier.Simplex.maxImprvLP`](@ref)
"""
function SimplexLP(PS::Problem{T}; settings=Settings(PS), min=true) where {T}
    (; E, u, d, G, g, A, b, N, M, J) = PS
    (; tol, rule) = settings

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    end

    Ms = M + J
    Ns = N + J
    S = fill(DN, Ms + Ns)
    B = collect(Ns .+ (1:Ms))
    S[B] .= IN

    As = [A zeros(T, M, J)
        G Matrix{T}(I, J, J)]
    invB = Matrix{T}(I, Ms, Ms)
    Es = [E; zeros(T, J)]
    bs = [b; g]
    ds = [d; zeros(T, J)]
    us = [u; fill(Inf, J)]

    #display((N, M, J, Ms, Ns))

    q = As * ds
    for j in 1:Ms
        invB[j, j] = bs[j] >= q[j] ? one(T) : -one(T)
    end
    q = abs.(As * ds - bs)
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #我的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    f, x, q, B, invB, iH = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    if f > tol
        error("feasible region is empty")
    end

    #display("--- --- phase 2 --- ---")
    S = S[1:Ns]
    #sgn = min == true ? 1 : -1
    if !min
        Es = -Es
    end
    #f, x, q, B, invB, iH = solveLP(-Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)
    f, x, q, B, invB, iH = solveLP(Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)
    if !min
        f = -f
    end

    for k in N+1:N+J
        S[k] = S[k] == IN ? OE : EO
    end

    return S, iH, q, f
end


function WolfeLP!(L, c, A, b, d, B, S; invB, q, tol=2^-26)
    #Complementary slackness LP, #no upper bound, only lower bound
    #B is always sorted. S in caller is change, compute invB and q=xB each step, using Bland rule to speed up

    T = typeof(c[1])
    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)
    ip = zeros(Int64, M)    #tracking the rows of p

    C = trues(N)
    x = zeros(T, N)
    x[S.==DN] .= d[S.==DN]

    Y = invB * A[:, F]
    h = c[F] - Y' * c[B]
    ih = (h .< -tol) .&& C[F]

    iH = findall(F)[ih]
    nH = length(iH)
    while nH > 0
        k = iH[1]
        p = invB * A[:, k]
        m = 0
        for j in 1:M
            i = B[j]
            if p[j] > tol
                m += 1
                gt[m] = (q[j] - d[i]) / p[j]
                ip[m] = j
            end
        end

        if m == 0
            error("infeasible or unbounded")
        end

        l = getfield(findmin(gt[1:m]), 2)
        m = ip[l]
        l = B[m]       #leaving index
        F[k] = false
        F[l] = true
        B[m] = k

        #Complementary slackness
        if k <= 2 * L
            C[k] = false
            C[k <= L ? k + L : k - L] = false
        end
        if l <= 2 * L
            C[l] = true
            C[l <= L ? l + L : l - L] = true
        end


        #B .= sort(B)
        sort!(B)
        invB = inv(A[:, B])
        S[k] = IN
        S[l] = DN
        x[l] = d[l]
        Y = invB * A[:, F]
        q = invB * b - Y * x[F]
        h = c[F] - Y' * c[B]
        ih = (h .< -tol) .&& C[F]
        iH = findall(F)[ih]
        nH = length(iH)
    end
    return nothing
end

function SimplexQP!(mu, aCL::Vector{sCL{T}}, PS::Problem{T}; nS=SettingsEF(PS), settings=Settings(PS)) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = PS


    #slack variable to Gz≤g, convert to Gz+s=g
    Vs = [V zeros(T, N, J)
        zeros(T, J, N) zeros(T, J, J)]
    As = [A zeros(T, M, J)
        G Matrix{T}(I, J, J)
        E' zeros(T, 1, J)]
    bs = [b; g; mu]
    ds = [d; zeros(T, J)]
    us = [u; fill(Inf, J)]
    Ms = M + J + 1
    Ns = N + J

    #convert upper bound to equality
    iu = findall(us .< Inf)
    nu = length(iu)
    Vu = [Vs zeros(T, Ns, nu)
        zeros(T, nu, Ns) zeros(T, nu, nu)]
    Au = [As zeros(T, Ms, nu)
        Matrix{T}(I, Ns, Ns)[iu, :] Matrix{T}(I, nu, nu)]
    bu = [bs; us[iu]]
    du = [ds; zeros(T, nu)]
    Mu = Ms + nu
    Nu = Ns + nu

    #KKT conditions
    A0 = [Au zeros(T, Mu, Nu + 2 * Mu)
        Vu -Matrix{T}(I, Nu, Nu) -Au' Au']
    b0 = [bu; zeros(T, Nu)]
    d0 = [du; zeros(T, Nu + 2 * Mu)]
    M0 = Mu + Nu
    N0 = 2 * M0 #(Nu + Mu)

    #introduction of artificial variables
    b1 = b0
    d1 = [d0; zeros(T, M0)]
    S1 = fill(DN, M0 + N0)
    B = collect(N0 .+ (1:M0))
    S1[B] .= IN
    invB = Matrix{T}(I, M0, M0)
    A1 = [A0 invB]
    #q = A1 * d1
    q = A0 * d0
    for j in 1:M0
        invB[j, j] = b0[j] >= q[j] ? one(T) : -one(T)
    end
    q = abs.(q - b0)
    c1 = [zeros(T, N0); fill(one(T), M0)]
    #M1 = M0
    #N1 = 3*M0

    WolfeLP!(Nu, c1, A1, b1, d1, B, S1; invB=invB, q=q, tol=settings.tol)


    S = S1[1:Ns]

    #display((N0, B, S))
    #=  Remark: if the only non-zero AV (artificial variable) is the one  for z′μ=μ, it may not be a problem. The KKT are hold, but for other mu value. 
    very slow, see N=263 in EF-dev-0108.jl
    WolfeLP! has a good chance only if mu is very close to the highest expected return, as in our case.  It doese not work 100%.
      It does not work for general QP. Since the K (IN, not on boundary) in a QP is variable, but LP has a fixed B (on boundary if degenerated).
    =#

    m = 1
    for k in iu
        if S[k] == IN
            S[k] = S1[Ns+m] == IN ? IN : UP
        end
        m += 1
    end

    for k in N+1:N+J
        S[k] = S[k] == IN ? OE : EO
    end
    computeCL!(aCL, S, PS, nS)

end




"""

        SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where T

compute the Critical Line Segments by Simplex method, for the highest expected return. Save the CL to aCL if done

Since `SimplexQP!` may fail, we give it up. Now using `EfficientFrontier.SimplexCL!` where `LightenQP`'s  QP solver is employed.

"""
function SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=SettingsEF(PS), settingsLP=Settings(PS), kwargs...) where {T}
    #function SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=SettingsEF(PS), settings=Settings(PS), settingsLP=Settings(PS)) where {T}
    #(; u, d, N, M, J) = PS
    #(; muShft, tolS, rule) = nS

    (S, iH, q, f) = SimplexLP(PS; settings=settingsLP, min=false)
    if computeCL!(aCL, S, PS, nS)
        return true
    end
    #display("------- SimplexQP!  -------")
    #mu = f * (nS.muShft - 1)    #assume Highest Mean > 0
    # mu = -f
    #sgn = (mu >= 0 ? -1 : 1)
    shft = nS.muShft
    if mu < -1 || mu > 1
        #sgn *= mu
        shft *= abs(mu)
    end
    #mu += sgn * nS.muShft
    mu -= shft

    return SimplexQP!(mu, aCL, PS; nS=nS, settings=settingsLP)
end


function oneCL!(aCL::Vector{sCL{T}}, PS::Problem{T}, S0::Vector{Status}, iH; nS=SettingsEF(PS)) where {T}
    #enumerating the states in iH for S0
    (; u, N) = PS
    iv = iH[iH.<=N]
    nk = length(iv)
    nj = length(iH[iH.>N])
    cbk = 1:nk
    cbj = 1:nj
    B = trues(nk)
    S = fill(IN, nk + nj)
    Sv = @view S[1:nk]
    Se = @view S[nk.+cbj]
    Se .= EO
    uk = u[cbk]
    for k in 1:nk
        cbK = combinations(cbk, k)
        for ii in cbK       #for IN
            S[ii] .= IN
            B[ii] .= false
            ic = (B .& (uk .< Inf))
            cbb = @view cbk[ic]
            nb = length(cbb)

            if nb == 0
                Sv[B] .= DN
                if nj == 0
                    S0[iH] = S
                    # display((iH, ii, S0))
                    if computeCL!(aCL, S0, PS, nS)
                        return true
                    end
                    B[ii] .= true
                    continue
                end
                # nj > 0
                for j in 0:nj
                    cbJ = combinations(cbj, j)
                    for ij in cbJ
                        Se[ij] .= OE
                        S0[iH] = S
                        # display((iH, ii, S0))
                        if computeCL!(aCL, S0, PS, nS)
                            return true
                        end
                        Se[ij] .= EO
                    end
                end
                B[ii] .= true
                continue
            end

            # nb > 0
            for m in nb:-1:0    #bounded
                cbB = combinations(cbb, m)
                for iu in cbB   #for UP
                    B[iu] .= false  #D
                    Sv[B] .= DN
                    S[iu] .= UP
                    if nj == 0
                        S0[iH] = S
                        # display((iH, ii, ij, S0))
                        if computeCL!(aCL, S0, PS, nS)
                            return true
                        end
                        B[iu] .= true
                        continue
                    end
                    # nj > 0
                    for j in 0:nj
                        cbJ = combinations(cbj, j)
                        for ij in cbJ
                            Se[ij] .= OE
                            S0[iH] = S
                            # display((iH, ii, ij, S0))
                            if computeCL!(aCL, S0, PS, nS)
                                return true
                            end
                            Se[ij] .= EO
                        end
                    end
                    B[iu] .= true
                end
            end
            B[ii] .= true
        end
    end
    return false
end

function LPcbCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=SettingsEF(PS), settings=Settings(PS), settingsLP=Settings(PS)) where {T}
    #Simplex LP + combinations
    (; u, d, N, M, J) = PS
    #(; tolS, rule) = nS
    tolS = nS.tolS

    (S, iH, q, f) = SimplexLP(PS; settings=settingsLP, min=false)
    j0 = sum(S[N+1:N+J] .== EO)

    #iH = iH[iH .<= N]   #It seems shoud be do so
    nH = length(iH)
    if nH == 0  #unique solution
        if !computeCL!(aCL, S, PS, nS)
            #check if hit a boundary point
            ip = findall(S .== IN)
            n = 1
            for i in ip
                t = q[n]
                if abs(t - d[i]) < tolS     #isapprox(t, d[i])
                    S[i] = DN
                elseif abs(t - u[i]) < tolS
                    S[i] = UP
                end
                n += 1
            end
            ip = findall(S .== IN)
            n = length(ip)
            if n == 0   #a boundary point
                cbk = 1:N
                cbK = combinations(cbk, max(2, M + j0))   #this works for j0==0, if j0>0, it will be complicated, the EO may turn OE when DN/UP go IN
                #cbK = combinations(cbk, 1)
                Sb = copy(S)
                for ii in cbK       #for IN
                    S[ii] .= IN
                    #display((ii, S))
                    if computeCL!(aCL, S, PS, nS)
                        return true
                    end
                    S[ii] = Sb[ii]
                end
            else    #NOT a boundary point
                return computeCL!(aCL, S, PS, nS)
            end
        end

    else
        aH = sort([iH; findall(S .== IN)])
        return oneCL!(aCL, PS, S, aH; nS=nS)
    end
    return false
end

end