"Simplex algorithm for EfficientFrontier"
module Simplex
using LinearAlgebra, Combinatorics
using EfficientFrontier: EfficientFrontier, Problem, Status, Event, sCL, IN, DN, UP, OE, EO, computeCL!
import EfficientFrontier: Settings as SettingsEF    #if by `using`, the compiler will complain conflict when import Simplex.Settings in to EfficientFrontier as SettingsLP
export SimplexLP, cDantzigLP, maxImprvLP

#=
REMARK: writing the next basis as a product of the current basis times an easily invertible matrix can be extended over several iterations. We do not adopt this for accuracy
    If this method is adopt, how often should one recompute an inverse of the current basis?
=#

"""

        Settings(P::Problem; kwargs...)        The default Settings to given Problem
        Settings(; kwargs...)       The default Settings is set by Float64 type
        Settings{T<:AbstractFloat}(; kwargs...)

kwargs are from the fields of Settings{T<:AbstractFloat} for Float64 and BigFloat

            tol::T          #2^-26 ≈ 1.5e-8  general scalar
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
the native model requires that `d` is finite
"""
function cDantzigLP(c::Vector{T}, A, b, d, u, B, S; invB, q, tol=2^-26) where {T}
    #B is always sorted. B and S in caller is change, compute invB and q=xB each step, switch Dantzig to Bland rule when iter > N

    #T = typeof(c[1])
    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)
    ip = zeros(Int, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    cA = zeros(T, N)    #norm of each column of A
    for k in 1:N
        cA[k] = norm(A[:, k])
    end
    #ldr = all(cA .> tol)    #some cols maybe zero. EfficientFrontier do not have this prob

    x = zeros(T, N)
    x[S.==UP] .= u[S.==UP]
    x[S.==DN] .= d[S.==DN]
    ud = u - d
    du = -ud
    fu = u .< Inf   #finite upper bound

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
    @inbounds while nH > 0
        loop += 1
        if loop > N
            Bland = true    #Dantzig rule switch to Bland rule
        end

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
                elseif p[j] < -tol #&& fu[i]
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                end
            end

            if m == 0   # p=0 => A[:,k]=0
                if fu[k]    #DN -> UP
                    l = -1
                else    # unbounded
                    #return c' * x, x, invB, 3
                    return 3, x, invB
                end
            else
                (gl, l) = findmin(gt[1:m])  #gl>0
                if fu[k]
                    if gl >= ud[k]   #DN -> UP
                        l = -1
                    else
                        Sl = Sb[l]
                        l = ip[l]
                    end
                else
                    if isinf(gl) #unbounded
                        #return c' * x, x, invB, 3
                        return 3, x, invB
                    end
                    Sl = Sb[l]
                    l = ip[l]
                end
            end

        else    #UP
            for j in 1:M
                i = B[j]
                if p[j] > tol #&& fu[i]
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

            if m == 0   # p=0 => A[:,k]=0
                l = -2  #UP -> DN
            else
                (gl, l) = findmax(gt[1:m])  #gl<0
                if gl <= du[k]  #du[k] is finite, for S[k]==UP and d is finite
                    l = -2  #UP -> DN
                else
                    Sl = Sb[l]
                    l = ip[l]
                end
            end
        end


        if l == -1  #flip the sate
            S[k] = UP
            x[k] = u[k]
            #q -= g0 * p
        elseif l == -2  #flip the sate
            S[k] = DN
            x[k] = d[k]
            #q -= g0 * p
        elseif l > 0
            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            #B = sort(B)
            sort!(B)
            invB = inv(A[:, B])

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            #q = invB * b - Y * x[F]
        end

        q = invB * b - Y * x[F]
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
    #iH = findall(F)[ih]
    x[B] = q

    status = length(ih) > 0 ? 2 : 1
    #return c' * x, x, invB, 2  #infinitely many solutions
    #return c' * x, x, invB, status
    return status, x, invB
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
the native model requires that `d` is finite
"""
function maxImprvLP(c::Vector{T}, A, b, d, u, B, S; invB, q, tol=2^-26) where {T}
    #greatest improvement, B is always sorted. B and S in caller is change, compute invB and q=xB each step
    #T = typeof(c[1])
    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)    #theta
    ip = zeros(Int, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    x = zeros(T, N)
    x[S.==UP] .= u[S.==UP]
    x[S.==DN] .= d[S.==DN]

    Y = invB * A[:, F]
    h = c[F] - Y' * c[B]
    vS = S[F]   #State of leaving to be, for each candidate k
    g = zeros(T, N - M)    #theta for each candidate k
    ig = zeros(Int, N - M)    #min subscrip for theta for each candidate k
    #vD = falses(N - M)  #DN or not, for each candidate k
    ud = u - d
    du = -ud
    fu = u .< Inf   #finite upper bound

    ih = S[F] .== DN
    h[ih] .= -h[ih]
    ih = h .> tol
    iH = findall(F)[ih]
    nH = length(iH)
    @inbounds while nH > 0
        P = @view Y[:, ih]
        for n in 1:nH
            k = iH[n]
            p = P[:, n]
            kd = S[k] == DN
            #vD[n] = kd
            m = 0
            if kd
                for j in 1:M
                    i = B[j]
                    if p[j] > tol
                        m += 1
                        gt[m] = (q[j] - d[i]) / p[j]
                        ip[m] = j
                        Sb[m] = DN
                    elseif p[j] < -tol #&& u[i] < Inf
                        m += 1
                        gt[m] = (q[j] - u[i]) / p[j]
                        ip[m] = j
                        Sb[m] = UP
                    end
                end

                g0 = ud[k]
                if m == 0
                    if fu[k]    #DN -> UP
                        g[n] = g0
                        l = 1
                        ip[1] = -1
                    else    # unbounded
                        #return c' * x, x, invB, 3
                        return 3, x, invB
                    end
                else
                    (g[n], l) = findmin(gt[1:m])

                    if fu[k]
                        if g[n] >= g0 #DN -> UP
                            g[n] = g0
                            ip[l] = -1
                        end
                    else
                        if isinf(g[n])  #unbounded
                            #return c' * x, x, invB, 3
                            return 3, x, invB
                        end
                    end
                end

            else    #UP
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
                g0 = du[k]  #finite
                if m == 0
                    g[n] = g0
                    l = 1
                    ip[1] = -2  #UP -> DN
                else
                    (g[n], l) = findmax(gt[1:m])
                    if g[n] <= g0
                        g[n] = g0
                        ip[l] = -2
                    end
                end
            end
            ig[n] = ip[l]
            vS[n] = Sb[l]
        end
        k = getfield(findmax(abs.(h[ih] .* g[1:nH])), 2)
        l = ig[k]
        #= if l < 0 && isinf(u[k])
            error("infeasible or unbounded")
        end =#

        #kd = vD[k]
        p = P[:, k]
        #gl = g[k]
        Sl = vS[k]
        k = iH[k]   #entering index
        if l == -1  #flip the sate
            S[k] = UP
            x[k] = u[k]
            #q -= gl * p
        elseif l == -2  #flip the sate
            S[k] = DN
            x[k] = d[k]
            #q -= gl * p
        elseif l > 0
            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            #B = sort(B)
            sort!(B)
            invB = inv(A[:, B])

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            #q = invB * b - Y * x[F]
        end

        q = invB * b - Y * x[F]
        h = c[F] - Y' * c[B]
        ih = S[F] .== DN
        h[ih] .= -h[ih]
        ih = h .> tol
        iH = findall(F)[ih]
        nH = length(iH)
    end

    ih = abs.(h) .< tol   # h==0
    #iH = findall(F)[ih]
    x[B] = q

    status = length(ih) > 0 ? 2 : 1
    #return c' * x, x, invB, 2  #infinitely many solutions
    #return c' * x, x, invB, status
    return status, x, invB

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

    Ms = M + J  #convert Gz<=g to equality
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
    #q = abs.(As * ds - bs)
    q = abs.(q - bs)
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #灯塔的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    #f, x, q, B, invB, iH = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    #f, x, invB, iH = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    iH, x, invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    #f = x' * c1
    f = sum(x[Ns+1:end])
    if f > tol
        #error("empty feasible region")
        return S, 0, x  #, f
    end

    #display("--- --- phase 2 --- ---")
    S = S[1:Ns]
    #sgn = min == true ? 1 : -1
    if !min
        Es = -Es
    end
    q = x[B]
    #f, x, q, B, invB, iH = solveLP(-Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)
    #f, x, q, B, invB, iH = solveLP(Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)
    #f, x, invB, iH = solveLP(Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)
    iH, x, invB = solveLP(Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)
    #= if !min
        f = -f
    end =#

    for k in N+1:N+J
        S[k] = S[k] == IN ? OE : EO
    end

    x = x[1:N]
    return S, iH, x #, f
end




end
