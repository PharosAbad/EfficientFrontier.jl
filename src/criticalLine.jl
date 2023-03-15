
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
    if K < size(AE, 1)    # impossible
        return false    # in our setting, size(AE, 1)>=1, hence skip K=0
    end
    # K >= size(AE, 1)  since the feasible set is not empty


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

    ik = findall(F)
    #= if K == 1
        k = ik[1]
        if abs(alpha[1]-d[k]) < tolG || abs(alpha[1]-u[k]) < tolG   #a boundary point
            return false
        end
    end =#

    #α_{λ}=-C(T′c+b_{E})
    alphaL = -(TC' * c + C * bE)
    #β_{λ}=CT′μ_{I}
    betaL = TC' * EF

    #η_{B}==V_{BI}α+V_{B}z_{B}+A_{B}′α_{λ}+L(V_{BI}β-μ_{B}+A_{B}′β_{λ})
    gamma = VBF * alpha + V[B, B] * zB + AB' * alphaL
    delta = VBF * beta - E[B] + AB' * betaL

    LE0 = Vector{Event{T}}(undef, 0)
    LE1 = Vector{Event{T}}(undef, 0)
    #ik = findall(F)
    @inbounds for k in eachindex(alpha)
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
            if !(d[j] + tolG < h < u[j] - tolG) #d[j]-tolG < h < d[j]+tolG means it is on the boundary!, so to u[j]
                #if !(d[j]-tolG < h < u[j]+tolG)
                return false
            end

        end
    end

    ib = findall(B)
    @inbounds for k in eachindex(gamma)
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
            if h <= -tolG   #z_{oA}≥0
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
        #LE0 = sort(LE0, by=x -> x.L, rev=true)
        sort!(LE0, by=x -> x.L, rev=true)
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
        #LE1 = sort(LE1, by=x -> x.L)
        sort!(LE1, by=x -> x.L)
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

        aCL = ECL(PS::Problem; init::Function=SimplexCL!; numSettings = Settings(PS), settings=SettingsQP(PS), settingsLP=SettingsLP(PS))

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
function ECL(PS::Problem; init::Function=SimplexCL!, numSettings=Settings(PS), settings=SettingsQP(PS), settingsLP=SettingsLP(PS))
    aCL = sCL(PS)
    init(aCL, PS; nS=numSettings, settings=settings, settingsLP=settingsLP)
    if lastindex(aCL) == 0
        error("Not able to find the first Critical Line")
    end

    ECL!(aCL, PS; numSettings=numSettings, settings=settings, settingsLP=settingsLP, incL=true)
    ECL!(aCL, PS; numSettings=numSettings, settings=settings, settingsLP=settingsLP)
    return aCL
end



function badK(S, PS, tolNorm)    #infeasible to compute CL
    (; u, d, G, g, A, b, N, J) = PS

    Sv = @view S[1:N]
    F = (Sv .== IN)
    K = sum(F)
    if K == 0   #degenerated
        return true, 0  # Status `S` from nextS, no need to check feasiblity of equality constraints
    end

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
        return true, -1    #infeasible
    else
        AE = AE[rb, :]
    end

    if K >= size(AE, 1) #good
        return false, K
    else
        @warn "Ignoring impossible: K < size(AE, 1)"
        return true, -2    #impossible
    end
end



"""

       ECL!(aCL::Vector{sCL{T}}, PS::Problem{T}; numSettings = Settings(PS), incL=false, settings=SettingsQP(PS), settingsLP=SettingsLP(PS))  where T

compute all the Critical Line Segments as L decreasing to 0 (increasing to +Inf if incL=true), given the first one in aCL[end].  The `settings` is the `LightenQP.Settings` setting from LightenQP
Return value: true if done

"""
function ECL!(aCL::Vector{sCL{T}}, PS::Problem{T}; numSettings=Settings(PS), incL=false, settings=SettingsQP(PS), settingsLP=SettingsLP(PS), endCL=true, f=Inf) where {T}
    tolNorm = numSettings.tolNorm
    N = PS.N
    t = aCL[end]
    #f = Inf
    while true
        if incL && isempty(t.I1)
            break
        end
        if !incL && t.L0 <= 0
            break
        end
        S = nextS(t, incL, N)

        bad, K = badK(S, PS, tolNorm)
        if bad
            if K != 0 || !endCL  #infeasible or impossible if `bad && K != 0`
                break
            end

            #degenerated
            if incL
                if isinf(f)
                    f = getfield(SimplexLP(PS; settings=settingsLP, min=false), 4)
                end
                mu = getMu(PS, t, incL)
                if abs(mu - f) < numSettings.tol  #hit HMFP=HVEP
                    break
                end
            end
            L1 = incL ? Inf : t.L0 / 2
            L0 = incL ? t.L1 * 2 : 0
            p = sCL{T}(Vector{Status}(undef, 0), 0, L1, Vector{Event{T}}(undef, 0), L0, Vector{Event{T}}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0))
            #p =  sCL(PS)
            #push!(p, sCL{T}(Vector{Status}(undef, 0), 0, L1, Vector{Event{T}}(undef, 0), L0, Vector{Event{T}}(undef, 0), Vector{T}(undef, 0), Vector{T}(undef, 0)))

            p = adjoinCL(p, t, S, incL, PS, numSettings, settings, settingsLP, f)
            S .= p.S
        end

        if computeCL!(aCL, S, PS, numSettings)
            t = aCL[end]
            if isinf(t.L1) || t.L0 <= 0
                break
            end
        else
            display(Int(S)')
            error("Critical Line Not Finished")
        end
        S .= t.S
    end
    sort!(aCL, by=x -> x.L1, rev=true)
    return true
end
#keep the argument `settings=SettingsQP`, we may use it for asQP in future



function nextS(t, incL, N)
    S = copy(t.S)
    if incL
        q = t.I1
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
    return S
end

function adjoinCL(p, t, S, incL, PS, nS, settings, settingsLP, f)
    L1 = incL ? p.L0 : t.L0
    L0 = incL ? t.L1 : p.L1
    L = (L1 + L0) / 2
    xS = copy(S)
    while xS == S   #in the region of degenerated
        x = asQP(PS; settingsLP=settingsLP, L=L)
        xS = getSx(x, PS, nS)
        K = sum(xS .== IN)
        if K > 0
            break
        end
        #now K = 0
        if xS != S   #hit another degenerated
            xS = copy(S)
            L1 = incL ? L : L1  #move in, go L close to t
            L0 = incL ? L0 : L
        else    #K=0 && xS == S, still in the region of degenerated
            L1 = incL ? 2 * L1 : L      #move out, go L away from t
            L0 = incL ? L : 0.5 * L0
        end
        L = (L1 + L0) / 2
    end
    xCL = sCL(PS)
    computeCL!(xCL, xS, PS, nS)
    ECL!(xCL, PS; numSettings=nS, settings=settings, settingsLP=settingsLP, incL=!incL, endCL=false, f=f)
    p = incL ? xCL[end] : xCL[1]
    xS = nextS(p, !incL, PS.N)
    if xS != S      # not connected
        p = adjoinCL(p, t, S, incL, PS, nS, settings, settingsLP, f)
    end
    return p
end

function getMu(PS::Problem{T}, t, incL) where {T}
    (; E, u, d, N) = PS
    z = zeros(T, N)
    S = @view t.S[1:N]
    F = (S .== IN)
    U = (S .== UP)
    D = (S .== DN)
    z[D] = d[D]
    z[U] = u[U]
    z[F] = t.alpha + t.beta * (incL ? t.L1 : t.L0)
    return z' * E
end



"""

        cbCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), oneCL=true, K=PS.M+PS.J+1, kwargs...) where T

compute one or all (oneCL=false, K=PS.M) the Critical Line Segments by enumerating (combinations of Status). Save the CL(s) to aCL if done

See also [`SimplexCL!`](@ref), [`Problem`](@ref), [`Settings`](@ref)
"""
function cbCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), oneCL=true, K=PS.M + PS.J + 1, kwargs...) where {T}
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



"""

        SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settingsLP=SettingsLP(PS)) where T

compute the Critical Line Segments by Simplex method, for the highest expected return. Save the CL to aCL if done

    settingsLP      : `SimplexLP.Settings`, for SimplexLP solver, we always first try the SimplexLP, and a chance of >=99.99% we find the critical line
                       when the remaining <=0.01%  happens (when the corner portfolio is degenerated, or infinite many solutions to SimplexLP encounted),
                       we use ASQP, which use LP to obtain an initial point, so no SettingsQP needed


See also [`cbCL!`](@ref), [`Problem`](@ref), [`Settings`](@ref), [`SettingsQP`](@ref)
"""
function SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settingsLP=SettingsLP(PS), kwargs...) where {T}
    #HMEP, HMFP (Highest Mean Frontier Portfolio), or HVEP (Highest Variance Efficient Portfolio), >=99.9% hit
    S, iH, q, f = SimplexLP(PS; settings=settingsLP, min=false)
    #display((f,S))
    if computeCL!(aCL, S, PS, nS)
        return true     #>=99.9% done
    end

    #test the LMFP (inefficient branch) >=0.09%
    rP = deepcopy(PS)
    E = rP.E
    E .= -E #negative
    S, iH, q, f = SimplexLP(rP; settings=settingsLP, min=false)
    xCL = sCL(rP)
    if computeCL!(xCL, S, rP, nS)
        ECL!(xCL, rP; numSettings=nS, settingsLP=settingsLP)
        #ECL!(xCL, rP; numSettings=nS, settings=settings, settingsLP=settingsLP)
        t = xCL[end]
        if t.L0 == 0    #GMVP
            return computeCL!(aCL, t.S, PS, nS)
            #return true
        end
    end

    #for the remaining <=0.01%
    asCL!(aCL, PS; nS=nS, settingsLP=settingsLP)
    #display("asCL!, buy lottery")

end
#=
function SimplexCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settingsLP=SettingsLP(PS), kwargs...) where {T}
    #HMEP, HMFP (Highest Mean Frontier Portfolio), or HVEP (Highest Variance Efficient Portfolio), >=99.9% hit
    S, iH, q, f = SimplexLP(PS; settings=settingsLP, min=false)
    if computeCL!(aCL, S, PS, nS)
        return true     #>=99.9% done
    end

    #x, aS = asQP(PS; settingsLP=settingsLP)   #GMVP
    #S = activeS(aS, PS)

    #asCL!  for the remaining <=0.1%
    #= x = asQP(PS; settingsLP=settingsLP)   #GMVP
    S = getSx(x, PS, nS)
    computeCL!(aCL, S, PS, nS)  =#

    asCL!(aCL, PS; nS=nS, settingsLP=settingsLP)    #for the remaining <=0.1%
end
=#



"""

        LightenCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsQP(PS), settingsLP=SettingsLP(PS)) where T

compute the Critical Line Segments by Simplex method, for the highest expected return. Save the CL to aCL if done

    settings        : for LightenQP solver, chance for a QP solver needed is <=0.1% (when the corner portfolio is degenerated, or infinite many solutions to SimplexLP encounted)
    settingsLP      : `SimplexLP.Settings`, for SimplexLP solver, we always first try the SimplexLP, and a chance of >=99.9% we find the critical line

See also [`cbCL!`](@ref), [`Problem`](@ref), [`Settings`](@ref), [`SettingsQP`](@ref)
"""
function LightenCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settings=SettingsQP(PS), settingsLP=SettingsLP(PS)) where {T}
    #using LightenQP

    S, iH, q, f = SimplexLP(PS; settings=settingsLP, min=false)  #HMFP (Highest Mean Frontier Portfolio), or HVEP (Highest Variance Efficient Portfolio), >=99.9% hit
    if computeCL!(aCL, S, PS, nS)
        return true
    end

    x, status = LightenQP.fPortfolio(OOQP(PS); settings=settings)   #GMVP, LVEP (Lowest Variance Efficient Portfolio),  may fail, leave it to computeCL!
    S = getS(x.s, PS, nS.tolS)

    computeCL!(aCL, S, PS, nS)
end



function getS(Y, PS, tolS)
    #from slack variables
    (; u, N, J) = PS
    iu = findall(u .< Inf)
    D = Y[(1:N).+J] .< tolS
    U = falses(N)
    U[iu] = Y[(1:length(iu)).+(J+N)] .< tolS
    S = fill(IN, N + J)
    Sv = @view S[1:N]
    Sv[D] .= DN
    Sv[U] .= UP
    if J > 0
        Se = @view S[(1:J).+N]
        Se .= EO
        Og = Y[1:J] .> tolS
        Se[Og] .= OE
    end
    return S
end

