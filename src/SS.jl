#imported functions
using StatusSwitchingQP.Simplex: cDantzigLP, maxImprvLP

"""

        SimplexLP(PS::Problem; settings=Simplex.Settings(PS), min=true)

find the `Status` for assets by simplex method. If `min=false`, we maximize the objective function

See also [`StatusSwitchingQP.Status`](@ref), [`Problem`](@ref), [`StatusSwitchingQP.Settings`](@ref), [`StatusSwitchingQP.Simplex.cDantzigLP`](@ref), [`StatusSwitchingQP.Simplex.maxImprvLP`](@ref)
"""
function SimplexLP(PS::Problem{T}; settings=SettingsLP(PS), min=true) where {T}
    #d is finite
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

    q = As * ds
    for j in 1:Ms
        invB[j, j] = bs[j] >= q[j] ? one(T) : -one(T)
    end
    q = abs.(q - bs)
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #灯塔的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    iH, x, invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    f = sum(x[Ns+1:end])
    if f > tol
        error("empty feasible region")
    end

    #display("--- --- phase 2 --- ---")
    q = x[B]
    ia = findall(B .> Ns)
    m = length(ia)

    while m > 0
        F = trues(Ns)
        F[B[B.<=Ns]] .= false
        Y = invB * As[:, F]
        l = B[end] - Ns
        if rank(Y) < Ms     # purge redundant row
            ir = trues(Ms)
            ir[l] = false
            Ms -= 1
            As = As[ir, :]
            bs = bs[ir]
            A1 = A1[ir, :]
            B = B[1:end-1]
            invB = inv(lu(A1[:, B]))
            q = q[1:end-1]
        else    #AV go out, replace by x[k]
            r = findfirst(abs.(Y[l, :]) .>= tol)
            k = findall(F)[r]
            B[end] = k
            ib = sortperm(B)
            B = B[ib]
            invB = inv(lu(A1[:, B]))
            q[end] = x[k]
            q = q[ib]
            S[k] = IN
        end
        ia = findall(B .> Ns)
        m = length(ia)
    end

    S = S[1:Ns]
    if !min
        Es = -Es
    end
    iH, x, invB = solveLP(Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)


    #= n = length(ia)
    if n == 0
        S = S[1:Ns]
        if !min
            Es = -Es
        end
        iH, x, invB = solveLP(Es, As, bs, ds, us, B, S; invB=invB, q=q, tol=tol)
    else
        #keep AV as BV (basic variable)
        xv = [collect(1:Ns); B[ia]]
        S = S[1:Ns+n]
        a1 = collect(Ns+1:Ns+n)
        S[a1] .= IN
        B[ia] .= a1

        c1 = [Es; fill(zero(T), n)]
        A1 = A1[:, xv]
        d1 = d1[xv]
        u1 = u1[xv]
        if !min
            c1 = -c1
        end
        iH, x, invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    end =#

    for k in N+1:N+J
        S[k] = S[k] == IN ? OE : EO
    end

    x = x[1:N]
    return x, S, iH

end


function QP(P::Problem{T}, L::T=0.0) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = P
    q = -L * E
    mc = 1
    return QP(V, A, G, q, b, g, d, u, N, M, J, mc)
end

function QP(mu::T, P::Problem{T}) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = P
    q = zeros(T, N)
    M += 1
    Am = [A; E']
    bm = [b; mu]
    mc = 1
    return QP(V, Am, G, q, bm, g, d, u, N, M, J, mc)
end