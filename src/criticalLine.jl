
function Int(v::Vector{Status})
    x = zeros(Int, length(v))
    for k in eachindex(v)
        x[k] = Int(v[k])
    end
    return x
end

#=
function ClarabelQP(E::Vector{T}, V::Matrix{T}, mu::T, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, Ae::Matrix{T}, be::Vector{T}) where {T}
    N = length(E)
    iu = u .!= Inf
    P = sparse(V)
    q = zeros(T, N) #Clarabel need P, q, A, b to be in type T
    A = sparse([E'; Ae; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    b = [mu; be; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(1 + length(be)), Clarabel.NonnegativeConeT(length(g) + N + sum(iu))]
    settings = Clarabel.Settings{T}()
    settings.verbose = false
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    Clarabel.solve!(solver)
end

function ClarabelLP(E::Vector{T}, u::Vector{T}, d::Vector{T}, G::Matrix{T}, g::Vector{T}, Ae::Matrix{T}, be::Vector{T}) where {T}
    N = length(E)
    iu = u .!= Inf  #Float64(Inf) == BigFloat(Inf)
    P = sparse(zeros(T, N, N))
    q = -E
    A = sparse([Ae; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]])
    b = [be; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(length(be)), Clarabel.NonnegativeConeT(length(g) + N + sum(iu))]
    settings = Clarabel.Settings{T}()
    settings.verbose = false
    solver = Clarabel.Solver{T}()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    Clarabel.solve!(solver)
end
=#

function getRows(A, nS)
    #indicate the non-redundant rows, the 1st row should be non-zeros (vector one here)
    H = @view A[1:1, :]
    R = falses(size(A, 1))
    R[1] = true
    for r in axes(A, 1)[(begin+1):end]
        v = @view A[r, :]
        if norm(v - H' * (H' \ v)) > nS.tolNorm
            R[r] = true
            H = @view A[R, :]
        end
    end
    return findall(R)
end

"""
        computeCL!(aCL::Vector{sCL{T}}, S::Vector{Status}, PS::Problem{T}, nS::Settings{T}) where T

compute the Critical Line Segment for S::Vector{Status}, save to aCL[end]. Return value: true if done.

"""
function computeCL!(aCL::Vector{sCL{T}}, S::Vector{Status}, PS::Problem{T}, nS::Settings{T}) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = PS
    (; tolL, tolG) = nS

    Sv = @view S[1:N]
    F = (Sv .== IN)
    B = .!F
    D = (Sv .== DN)
    U = (Sv .== UP)
    Se = @view S[(N.+(1:J))]
    Eg = (Se .== EO)
    Og = (Se .== OE)
    GE = @view G[Eg, :]
    AE = vcat(A[:, F], GE[:, F])
    idAE = vcat(axes(A, 1), axes(GE, 1)) # id of each constraint

    zB = d[B]
    zB[U[B]] = u[U]
    AB = vcat(A[:, B], GE[:, B])
    bE = vcat(b, g[Eg]) - AB * zB
    rb = getRows([AE bE], nS)
    if length(getRows(AE, nS)) != length(rb)
        return false    #infeasible
    else
        AE = AE[rb, :]
        bE = bE[rb]
        AB = AB[rb, :]
    end
    K = sum(F)
    if K < size(AE, 1)    #infeasible
        return false
    end

    EF = @view E[F]


    VF = @view V[F, F]
    iV = inv(cholesky(VF))  #make sure iV is symmetric
    #=
    https://discourse.julialang.org/t/inverse-of-a-symmetric-matrix-is-not-symmetric/10132/5
    https://groups.google.com/g/julia-users/c/QBXmLGQp3Mw/m/IPXaxh-6Fk0J	you are telling Julia that your matrix is positive definite and Julia will exploit that to give you a fast inverse which is also positive definite.
    https://scicomp.stackexchange.com/a/35645	best method is really application dependent

    VF = @view V[F, F]
    #iV = inv(VF)   #ERROR: MethodError: no method matching factorize(::SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Vector{Int64}}, false})    
    iV = inv(Symmetric(VF))
    =#
    VBF = @view V[B, F]
    #c = V[F, B] * zB  #c=V_{IB}z_{B}
    c = VBF' * zB  #c=V_{IB}z_{B}
    mT = iV * AE'   #T=V_{I}⁻¹A_{I}′
    #C = inv(AE * mT)   #C=(A_{I}T)⁻¹
    #C = inv(Symmetric(AE * mT))   #C=(A_{I}T)⁻¹    
    C = AE * mT
    C = (C + C') / 2
    C = inv(cholesky(C))
    #AL = AE*cholesky(iV).L
    #C = inv(cholesky(AL*AL'))  #ERROR: PosDefException: matrix is not Hermitian

    TC = mT * C
    VQ = iV - mT * TC'    #Q=I-A_{I}′CT′   V_{I}⁻¹Q=V_{I}⁻¹-TCT′

    alpha = TC * bE - VQ * c    #α=TCb_{E}-V_{I}⁻¹Qc
    beta = VQ * EF    #β=V_{I}⁻¹Qμ_{I}

    #α_{λ}=-C(T′c+b_{E})
    alphaL = -(TC' * c + C * bE)
    #β_{λ}=CT′μ_{I}
    betaL = TC' * EF

    #η_{B}==V_{BI}α+V_{B}z_{B}+A_{B}′α_{λ}+L(V_{BI}β-μ_{B}+A_{B}′β_{λ})
    gamma = VBF * alpha + V[B, B] * zB + AB' * alphaL
    delta = VBF * beta - E[B] + AB' * betaL

    LE0 = Vector{Event{T}}(undef, 0)
    LE1 = Vector{Event{T}}(undef, 0)
    ik = findall(F)
    for k in eachindex(alpha)
        j = ik[k]
        t = beta[k]
        h = alpha[k]
        dL = (d[j] - h) / t
        uL = (u[j] - h) / t
        if t > tolG
            push!(LE0, Event{T}(IN, DN, j, dL))
            if u[j] < Inf
                push!(LE1, Event{T}(UP, IN, j, uL))
            end
        elseif t < -tolG
            push!(LE1, Event{T}(DN, IN, j, dL))
            if u[j] < Inf
                push!(LE0, Event{T}(IN, UP, j, uL))
            end
        else
            if !(d[j] < h < u[j])
                return false
            end

        end
    end

    ib = findall(B)
    for k in eachindex(gamma)
        j = ib[k]
        t = delta[k]
        h = gamma[k]
        dL = -h / t
        if t > tolG
            if D[j]
                push!(LE0, Event{T}(DN, IN, j, dL))
            else
                push!(LE1, Event{T}(IN, UP, j, dL))
            end
        elseif t < -tolG
            if D[j]
                push!(LE1, Event{T}(IN, DN, j, dL))
            else
                push!(LE0, Event{T}(UP, IN, j, dL))
            end
        else
            if (D[j] && h <= 0) || (U[j] && h >= 0)
                return false
            end
        end
    end

    #check the IO of inequalities
    #z_{o} =(g_{O}-G_{OB}z_{B}-G_{OI}α)-LG_{OI}β
    zoA = g[Og] - G[Og, B] * zB - G[Og, F] * alpha
    zoB = -G[Og, F] * beta

    ik = findall(Og)
    for k in eachindex(zoA)
        j = ik[k]
        t = zoB[k]
        h = zoA[k]
        dL = -h / t
        if t > tolG
            push!(LE0, Event{T}(OE, EO, j, dL))
        elseif t < -tolG
            push!(LE1, Event{T}(EO, OE, j, dL))
        else
            if h <= 0
                return false
            end
        end
    end

    #restore full ldE if refined    
    JE = size(GE, 1)
    if JE > 0
        iE = zeros(Int, JE)
        iR = findall(rb .> M)
        iE[idAE[rb[iR]]] = iR
        Lda = zeros(JE)
        Ldb = zeros(JE)

        for j in 1:JE
            k = iE[j]
            if k == 0
                x = AE' \ GE[j, F]
                Lda[j] = alphaL' * x
                Ldb[j] = betaL' * x
            else
                Lda[j] = alphaL[k]
                Ldb[j] = betaL[k]
            end
        end

        ib = findall(Eg)
        for k in 1:JE
            j = ib[k]
            t = Ldb[k]
            h = Lda[k]
            dL = -h / t
            if t > tolG
                push!(LE0, Event{T}(EO, OE, j, dL))
            elseif t < -tolG
                push!(LE1, Event{T}(OE, EO, j, dL))
            else
                if h <= 0
                    return false
                end
            end
        end
    end

    L0::T = -Inf
    I0 = Vector{Event{T}}(undef, 0)
    nL0 = length(LE0)
    if nL0 > 0
        LE0 = sort(LE0, by=x -> x.L, rev=true)
        ic0 = LE0[1]
        L0 = ic0.L
        I0 = [ic0]
    end

    if L0 < -tolL
        L0 = 0.0
    end

    L1::T = Inf
    I1 = Vector{Event{T}}(undef, 0)
    nL1 = length(LE1)
    if nL1 > 0
        LE1 = sort(LE1, by=x -> x.L)
        ic1 = LE1[1]
        L1 = ic1.L
        I1 = [ic1]
    end
    if L1 < L0 + tolL
        return false
    end

    if nL0 > 1 && L0 - LE0[2].L < tolL
        I0 = [I0; LE0[2]]
        for k in 3:nL0
            if L0 - LE0[k].L < tolL
                I0 = [I0; LE0[k]]
            else
                break
            end
        end
    end

    if nL1 > 1 && LE1[2].L - L1 < tolL
        I1 = [I1; LE1[2]]
        for k in 3:nL1
            if LE1[k].L - L1 < tolL
                I1 = [I1; LE1[k]]
            else
                break
            end
        end
    end

    push!(aCL, sCL{T}(copy(S), K, L1, I1, L0, I0, alpha, beta))
    return true
end



#=
function BoundaryCP!(aCL::Vector{sCL{T}}, S::Vector{Status}, PS::Problem{T}, nS::Settings{T}) where {T}
    #S is modified too
    (; E, V, u, d, G, g, A, b, N, M, J) = PS
    (; tolS, muShft) = nS

    Sv = @view S[1:N]
    D = (Sv .== DN)
    U = (Sv .== UP)
    z = zeros(T, N)
    z[D] = d[D]
    z[U] = u[U]
    mu = z' * E * (1 - muShft)

    y = ClarabelQP(E, V, mu, u, d, G, g, A, b)
    if Int(y.status) != 1   #SOLVED
        error("Not able to find the efficient portfolio at mean: " * sprint(show, mu))
    end

    Y = y.s
    iu = u .!= Inf
    D .= Y[(1:N).+(M+J+1)] .< tolS
    U .= false
    U[iu] = Y[(1:sum(iu)).+(M+J+N+1)] .< tolS
    #S = fill(IN, N + J)
    S .= IN
    Sv = @view S[1:N]
    Sv[D] .= DN
    Sv[U] .= UP
    if J > 0
        Se = @view S[(1:J).+N]
        Se .= EO
        Og = Y[(1:J).+(M+1)] .> tolS
        Se[Og] .= OE
    end
    computeCL!(aCL, S, PS, nS)
end
=#

"""

        aCL = ECL(PS::Problem, init=:Combinatorics; numSettings = Settings(PS))

compute all the Critical Line Segments. Init the CL by combinatorial search (`init=:Clarabel` by Clarabel.jl, an interior point numerical solver)

the return aCL has the follwing structure

        struct sCL{T<:AbstractFloat}    #critical line segment
            S::Vector{Status}       # (N+J)x1
            K::Integer              #number of IN assets
            L1::T                   #higher lambda
            I1::Vector{Event{T}}    #go in/out events at L1
            L0::T                   #lower lambda
            I0::Vector{Event{T}}    #go in/out events at L0
            alpha::Vector{T}        # K x 1
            beta::Vector{T}         # K x 1
        end

"""
function ECL(PS::Problem, init=:Combinatorics; numSettings=Settings(PS))
    aCL = Vector{sCL{typeof(PS).parameters[1]}}(undef, 0)
    if init == :Clarabel
        #return cECL(PS; numSettings=numSettings)

        #=if !initCL!(aCL, PS, numSettings)
            error("Not able to find the first Critical Line")
        end =#
        initCL!(aCL, PS, numSettings)
    else
        cbCL!(aCL, PS; nS=numSettings)
    end
    if lastindex(aCL) == 0
        error("Not able to find the first Critical Line")
    end

    ECL!(aCL, PS; numSettings=numSettings, incL=true)
    ECL!(aCL, PS; numSettings=numSettings)
    return aCL
end


function initCL!(aCL, PS, nS)
    (; E, V, u, d, G, g, A, b, N, M, J) = PS
    (; tolS, muShft) = nS

    x = ClarabelLP(E, u, d, G, g, A, b)
    if Int(x.status) != 1   #SOLVED
        error("Not able to find the max expected return of efficient frontier")
    end
    y = ClarabelQP(E, V, x.obj_val * (muShft - 1), u, d, G, g, A, b)
    if Int(y.status) != 1   #SOLVED
        error("Not able to find the max expected return efficient portfolio")
    end

    Y = y.s
    iu = u .!= Inf
    D = Y[(1:N).+(M+J+1)] .< tolS
    U = falses(N)
    U[iu] = Y[(1:sum(iu)).+(M+J+N+1)] .< tolS
    S = fill(IN, N + J)
    Sv = @view S[1:N]
    Sv[D] .= DN
    Sv[U] .= UP
    if J > 0
        Se = @view S[(1:J).+N]
        Se .= EO
        Og = Y[(1:J).+(M+1)] .> tolS
        Se[Og] .= OE
    end
    computeCL!(aCL, S, PS, nS)
end

#=
function cECL(PS::Problem; numSettings=Settings(PS))
#using ClarabelLP thr initCL!
    aCL = Vector{sCL{typeof(PS).parameters[1]}}(undef, 0)
    if !initCL!(aCL, PS, numSettings)
        error("Not able to find the first Critical Line")
    end
    t = aCL[end]
    while true
        S = copy(t.S)
        q = t.I0
        for k in eachindex(q)
            c = q[k]
            To = c.To
            id = c.id
            if To == EO || To == OE
                id += PS.N
            end
            S[id] = To
        end
        #display(Int(S)')
        if computeCL!(aCL, S, PS, numSettings)
            t = aCL[end]
            if t.L0 <= 0
                break
            end
        else
            display(Int(S)')
            error("Critical Line Not Finished")
        end
        #display(t)
        #error("Debug ...")
    end
    return aCL
end
=#

"""

       ECL!(aCL::Vector{sCL{T}}, PS::Problem{T}; numSettings = Settings(PS), incL=false)  where T

compute all the Critical Line Segments as L decreasing to 0 (increasing to +Inf if incL=true), given the first one in aCL[end]. 
Return value: true if done

"""
function ECL!(aCL::Vector{sCL{T}}, PS::Problem{T}; numSettings=Settings(PS), incL=false) where {T}
    (; N, M) = PS
    t = aCL[end]
    S = copy(t.S)
    while true
        #S = copy(t.S)
        if incL
            q = t.I1
            if isempty(q)
                break
            end
            for k in eachindex(q)
                c = q[k]
                From = c.From
                id = c.id
                if From == EO || From == OE
                    id += N
                end
                S[id] = From
            end
            if sum(S .== IN) < M + sum(S .== EO)
                break
            end
        else
            q = t.I0
            for k in eachindex(q)
                c = q[k]
                To = c.To
                id = c.id
                if To == EO || To == OE
                    id += N
                end
                S[id] = To
            end
        end

        if computeCL!(aCL, S, PS, numSettings)
            t = aCL[end]
            #= if incL
                if isinf(t.L1)
                    break
                end
                if aCL[end-1].S == t.S
                    display(Int(aCL[end-1].S)')
                    display(Int(t.S)')
                    pop!(aCL)
                    break
                end
            else
                if t.L0 <= 0
                    break
                end
            end =#
            if isinf(t.L1) || t.L0 <= 0
                break
            end

        else
            #return false
            display(Int(S)')
            error("Critical Line Not Finished")
        end
        S .= t.S
    end
    sort!(aCL, by=x -> x.L1, rev=true)
    return true
end


"""

        cbCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), oneCL=true, K=PS.M+PS.J+1) where T

compute one or all (oneCL=false, K=PS.M) the Critical Line Segments by enumerating (combinations of Status)


"""
function cbCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), oneCL=true, K=PS.M + PS.J + 1) where {T}
    #(3^N -2^N -N*2^(N-1)) * 2^J cases  ( K= 0 -> 2^N ; K=1 -> N*2^(N-1))
    (; u, N, J) = PS
    cbk = 1:N
    cbj = 1:J
    B = trues(N)
    S = fill(IN, N + J)
    Sv = @view S[1:N]
    Se = @view S[N.+cbj]
    Se .= EO

    for k in max(2, K):N    # k<PS.M+PS.J no degree of freedom
        cbK = combinations(cbk, k)
        for ii in cbK       #for IN
            S[ii] .= IN
            B[ii] .= false
            ic = (B .& (u .< Inf))
            cbb = @view cbk[ic]
            nb = length(cbb)

            if nb == 0
                Sv[B] .= DN
                if J == 0
                    if computeCL!(aCL, S, PS, nS) && oneCL
                        return nothing
                    end
                    B[ii] .= true
                    continue
                end
                # J > 0
                for j in 0:J
                    cbJ = combinations(cbj, j)
                    for ij in cbJ
                        Se[ij] .= OE
                        if computeCL!(aCL, S, PS, nS) && oneCL
                            return nothing
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
                    if J == 0
                        if computeCL!(aCL, S, PS, nS) && oneCL
                            return nothing
                        end
                        B[iu] .= true
                        continue
                    end
                    # J > 0
                    for j in 0:J
                        cbJ = combinations(cbj, j)
                        for ij in cbJ
                            Se[ij] .= OE
                            if computeCL!(aCL, S, PS, nS) && oneCL
                                return nothing
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
    return nothing
end
