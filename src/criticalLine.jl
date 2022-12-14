
function Int(v::Vector{Status})
    x = zeros(Int, length(v))
    for k in eachindex(v)
        x[k] = Int(v[k])
    end
    return x
end


function getRows(A, tolNorm)
    #indicate the non-redundant rows, the 1st row should be non-zeros (vector one here)
    H = @view A[1:1, :]
    R = falses(size(A, 1))
    R[1] = true
    for r in axes(A, 1)[(begin+1):end]
        v = @view A[r, :]
        if norm(v - H' * (H' \ v)) > tolNorm
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
    (; tolNorm, tolL, tolG) = nS

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
    rb = getRows([AE bE], tolNorm)
    if length(getRows(AE, tolNorm)) != length(rb)
        return false    #infeasible
    else
        AE = AE[rb, :]
        bE = bE[rb]
        AB = AB[rb, :]
    end
    K = sum(F)
    if K < size(AE, 1)    # <:  infeasible,  or =: no freedom
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
    mT = iV * AE'   #T=V_{I}?????A_{I}???
    #C = inv(AE * mT)   #C=(A_{I}T)?????
    #C = inv(Symmetric(AE * mT))   #C=(A_{I}T)?????    
    C = AE * mT
    C = (C + C') / 2
    C = inv(cholesky(C))
    #AL = AE*cholesky(iV).L
    #C = inv(cholesky(AL*AL'))  #ERROR: PosDefException: matrix is not Hermitian

    TC = mT * C
    VQ = iV - mT * TC'    #Q=I-A_{I}???CT???   V_{I}?????Q=V_{I}?????-TCT???

    alpha = TC * bE - VQ * c    #??=TCb_{E}-V_{I}?????Qc
    beta = VQ * EF    #??=V_{I}?????Q??_{I}

    ik = findall(F)
    #= if K == 1
        k = ik[1]
        if abs(alpha[1]-d[k]) < tolG || abs(alpha[1]-u[k]) < tolG   #a boundary point
            return false
        end
    end =#

    #??_{??}=-C(T???c+b_{E})
    alphaL = -(TC' * c + C * bE)
    #??_{??}=CT?????_{I}
    betaL = TC' * EF

    #??_{B}==V_{BI}??+V_{B}z_{B}+A_{B}?????_{??}+L(V_{BI}??-??_{B}+A_{B}?????_{??})
    gamma = VBF * alpha + V[B, B] * zB + AB' * alphaL
    delta = VBF * beta - E[B] + AB' * betaL

    LE0 = Vector{Event{T}}(undef, 0)
    LE1 = Vector{Event{T}}(undef, 0)
    #ik = findall(F)
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
            if !(d[j]+tolG < h < u[j]-tolG) #d[j]-tolG < h < d[j]+tolG means it is on the boundary!, so to u[j]
            #if !(d[j]-tolG < h < u[j]+tolG)
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
            if (D[j] && h <= -tolG) || (U[j] && h >= tolG)
            #if (D[j] && h < tolG) || (U[j] && h > -tolG)
                return false
            end
        end
    end

    #check the IO of inequalities
    #z_{o} =(g_{O}-G_{OB}z_{B}-G_{OI}??)-LG_{OI}??
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
            if h <= -tolG   #z_{oA}???0
            #if h < tolG
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
                #if h < tolG
                if h <= -tolG
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







"""

        aCL = ECL(PS::Problem; numSettings = Settings(PS), init::Function=SimplexCL!)

compute all the Critical Line Segments. Init the CL by Simplex Method (combinatorial search `init=cbCL!`)

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

See [`Documentation for EfficientFrontier.jl`](https://github.com/PharosAbad/EfficientFrontier.jl/wiki)

See also [`Problem`](@ref), [`Settings`](@ref)

"""
function ECL(PS::Problem; numSettings=Settings(PS), init::Function=SimplexCL!)
#function ECL(PS::Problem; numSettings=Settings(PS), init::Function=cbCL!)
    #aCL = Vector{sCL{typeof(PS).parameters[1]}}(undef, 0)
    aCL = sCL(PS)
    init(aCL, PS; nS=numSettings)
    if lastindex(aCL) == 0
        error("Not able to find the first Critical Line")
    end

    ECL!(aCL, PS; numSettings=numSettings, incL=true)
    ECL!(aCL, PS; numSettings=numSettings)
    return aCL
end



function badK(S, PS, tolNorm)    #infeasible to compute CL
    (; u, d, G, g, A, b, N, J) = PS

    Sv = @view S[1:N]
    F = (Sv .== IN)
    B = .!F
    U = (Sv .== UP)
    Se = @view S[(N.+(1:J))]
    Eg = (Se .== EO)
    GE = @view G[Eg, :]
    AE = vcat(A[:, F], GE[:, F])

    zB = d[B]
    zB[U[B]] = u[U]
    AB = vcat(A[:, B], GE[:, B])
    bE = vcat(b, g[Eg]) - AB * zB
    rb = getRows([AE bE], tolNorm)
    if length(getRows(AE, tolNorm)) != length(rb)
        return true    #infeasible
    else
        AE = AE[rb, :]
    end
    K = sum(F)
    if K < size(AE, 1)    # <:  infeasible,  or =: no freedom
        return true 
    end
    return false
end

"""

       ECL!(aCL::Vector{sCL{T}}, PS::Problem{T}; numSettings = Settings(PS), incL=false)  where T

compute all the Critical Line Segments as L decreasing to 0 (increasing to +Inf if incL=true), given the first one in aCL[end]. 
Return value: true if done

"""
function ECL!(aCL::Vector{sCL{T}}, PS::Problem{T}; numSettings=Settings(PS), incL=false) where {T}
    #(; N, M) = PS
    N = PS.N
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
        else
            if t.L0 <= 0
                break
            end
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
        if badK(S, PS, numSettings.tolNorm)
        #if sum(S .== IN) <= M + sum(S .== EO)
            break
        end

        if computeCL!(aCL, S, PS, numSettings)
            t = aCL[end]
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

compute one or all (oneCL=false, K=PS.M) the Critical Line Segments by enumerating (combinations of Status). Save the CL(s) to aCL if done


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
