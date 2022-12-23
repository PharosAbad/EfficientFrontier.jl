"Simplex algorithm"
module uSimplex
using LinearAlgebra, Combinatorics, EfficientFrontier
export SimplexCL!


#function LP0d1(c, A, b, d, u, B, S; invB=inv(A[:, sort(vec(B))]), q=invB * b, tol=2^-26)
function DantzigBlandLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)
    #B is always sorted. S in caller is change, compute invB and q=xB each step, switch Dantzig to Bland rule when iter > N 

    T = typeof(c[1])
    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)
    ip = zeros(Int64, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

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
        k0 = Bland ? 1 : argmax(hp)
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


#function LP0g1(c, A, b, d, u, B, S; invB=inv(A[:, sort(vec(B))]), q=invB * b, tol=2^-26)
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

    Solver = DantzigBlandLP
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
    #= for i in 1:Ms   #make sure bs>0
        if bs[i] < 0
            bs[i] = -bs[i]
            As[i, :] = -As[i, :]
        end
    end =#

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

    return S, iH, q

end

function oneCL!(aCL::Vector{sCL{T}}, PS::Problem{T}, S0::Vector{Status}, iH; nS=Settings(PS)) where {T}

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

function SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    (; u, d, N, M, J) = PS
    (; tolS, rule) = nS
    

    #(S, iH, q) = SimplexLP(PS, nS.tolS; rule=:maxImprovement)
    (S, iH, q) = SimplexLP(PS, tolS; rule=rule)
    #(S, iH, q) = SimplexLP(PS, nS.tolS)
    #display((S, iH))
    for i in 1:J
        n = N + i
        S[n] = S[n] == IN ? OE : EO
    end
    #return S, iH
    # display((S, iH, q))
    #iH = iH[iH .<= N]   #It seems shoud be do so
    nH = length(iH)
    if nH == 0
        #check if hit a boundary point
        ip = findall(S .== IN)
        n = 1
        for i in ip
            t = q[n]
            if abs(t-d[i]) < tolS
            #if isapprox(t, d[i])
                S[i] = DN
            elseif abs(t-u[i]) < tolS   #isapprox(t, u[i])
                S[i] = UP
            end
            n += 1
        end

        ip = findall(S .== IN)
        n = length(ip)

        if n == 0
            #return (S, iH)
            cbk = 1:N
            cbK = combinations(cbk, max(2, M))
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
            #display("NOT a boundary point")
            #display((n, S))
            return computeCL!(aCL, S, PS, nS)
        end

    else
        #= S = Status[DN, DN, DN, DN, DN, DN, DN, DN, DN, DN, DN, DN, DN, IN]
        iH = [9, 11]
        iH = sort([iH;findall(S .== IN)])
        nH = length(iH)
        Sb = copy(S)
        for k in M:nH
        #for k in 1:nH
            cbK = combinations(iH, k)
            #cbK = combinations(cbk, 1)
            for ii in cbK       #for IN
                S[ii] .= IN
                display((-ii, S))
                if computeCL!(aCL, S, PS, nS)
                    return true
                end
                S[ii] = Sb[ii]
            end
        end =#
        aH = sort([iH; findall(S .== IN)])        
        return oneCL!(aCL, PS, S, aH)
    end
    return false
    #return S, iH
end

end