"Simplex algorithm"
module uSimplex
using LinearAlgebra, Combinatorics, EfficientFrontier
export SimplexCL!

#=
REMARK: writing the next basis as a product of the current basis times an easily invertible matrix can be extended over several iterations. We do not adopt this for accuracy
    If this method is adopt, how often should one recompute an inverse of the current basis?
=#

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

            if m == 0
                error("infeasible or degenerate")
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
                error("infeasible or degenerate")
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
    return q, B, invB, iH, c' * x
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
    vD = falses(N - M)  #DN  or not for each candidate k
    ud = u - d

    ih = S[F] .== DN
    h[ih] .= -h[ih]
    ih = h .> tol
    hp = h[ih]
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
                    error("infeasible or degenerate")
                end
                (g[n], l) = findmin(gt[1:m])
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
                    error("infeasible or degenerate")
                end
                (g[n], l) = findmax(gt[1:m])
                if g[n] < -ud[k]
                    g[n] = -ud[k]
                end
            end
            gi[n] = ip[l]
            vS[n] = Sb[l]
        end
        k = getfield(findmax(abs.(h[ih] .* g[1:nH])), 2)
        l = gi[k]
        kd = vD[k]
        p = P[:, k]
        gl = g[k]
        Sl = vS[k]
        k = iH[k]   #entering index
        if kd
            if gl == ud[k]
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

    ih = abs.(h) .< tol   # h==0
    iH = findall(F)[ih]
    x[B] = q
    return q, B, invB, iH, c' * x
end

function SimplexLP(PS::Problem{T}, tol=sqrt(eps(T)); rule=:Dantzig) where {T}
    (; E, u, d, G, g, A, b, N, M, J) = PS

    Solver = cDantzigLP
    if rule == :maxImprovement
        Solver = maxImprvLP
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

    #(q, B, invB, iH) = DantzigBlandLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    (q, B, invB, iH, f) = Solver(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    if f > tol
        error("feasible region is empty")
    end

    #display("--- --- phase 2 --- ---")
    # q, B, invB, h
    S = S[1:Ns]


    (q, B, invB, iH, f) = Solver(-Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)

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
    ih = S[F] .== DN
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
            error("infeasible or degenerate")
        end

        l = getfield(findmin(gt[1:m]), 2)
        m = ip[l]
        l = B[m]       #leaving index
        F[k] = false
        F[l] = true
        B[m] = k

        #Complementary slackness
        if k <= 2 * L
            m = mod(k, L)
            m = m == 0 ? L : m
            C[m] = false
            C[m+L] = false
        end
        if l <= 2 * L
            m = mod(l, L)
            m = m == 0 ? L : m
            C[m] = true
            C[m+L] = true
        end


        B = sort(B)
        invB = inv(A[:, B])
        S[k] = IN
        S[l] = DN
        x[l] = d[l]
        Y = invB * A[:, F]
        q = invB * b - Y * x[F]
        h = c[F] - Y' * c[B]
        ih = S[F] .== DN
        ih = (h .< -tol) .&& C[F]
        iH = findall(F)[ih]
        nH = length(iH)
    end

    #ih = abs.(h) .< tol   # h==0
    #iH = findall(F)[ih]
    #x[B] = q
    #return q, B, invB, iH, x
    return nothing
end

function SimplexQP!(mu, aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
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
    N0 = 2 * (Nu + Mu)

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
    #N1 = N0 + M0

    WolfeLP!(Nu, c1, A1, b1, d1, B, S1; invB=invB, q=q, tol=nS.tol)


    S = S1[1:Ns]
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


"""
function SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    #(; u, d, N, M, J) = PS
    (; muShft, tolS, rule) = nS

    (S, iH, q, f) = SimplexLP(PS, tolS; rule=rule)
    if computeCL!(aCL, S, PS, nS)
        return true
    end

    mu = f * (muShft - 1)
    return SimplexQP!(mu, aCL, PS; nS=nS)
    #display((f,mu))
    #= nH = length(iH)
    if nH == 0
        if !computeCL!(aCL, S, PS, nS)
            return SimplexQP(mu, aCL, PS; nS=nS)
        end
        return true
    else
        return SimplexQP(mu, aCL, PS; nS=nS)
    end =#

end


function oneCL!(aCL::Vector{sCL{T}}, PS::Problem{T}, S0::Vector{Status}, iH; nS=Settings(PS)) where {T}
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

function LPcbCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    #Simplex LP + combinations
    (; u, d, N, M, J) = PS
    (; tolS, rule) = nS


    #(S, iH, q) = SimplexLP(PS, nS.tolS; rule=:maxImprovement)
    (S, iH, q, f) = SimplexLP(PS, tolS; rule=rule)
    #(S, iH, q) = SimplexLP(PS, nS.tolS)
    #display((S, iH))
    j0 = 0
    for i in 1:J
        n = N + i
        #S[n] = S[n] == IN ? OE : EO
        if S[n] == IN
            S[n] = OE
        else
            S[n] = EO
            j0 += 1
        end
    end
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
                #return (S, iH)
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
        return oneCL!(aCL, PS, S, aH)
    end
    return false
end

end