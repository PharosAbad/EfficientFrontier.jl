
"""

        aEF = eFrontier(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where T

compute the Full Efficient Frontier by Status-Segment Method

the return aEF has the follwing structure

        struct sEF    #Efficient Frontier       Float64 is OK
            Ap::Matrix{Float64}     # v=a₂μ²+a₁μ+a₀, each row [a₂ a₁ a₀]
            mu::Vector{Float64}     #higher mean
            sgm::Vector{Float64}    #higher sigma
            Z::Matrix{Float64}      #weights, each corner portfolio in one row
            ic::Vector{Int64}       #id of related critical line
        end

See also [`EfficientFrontier.ECL`](@ref), [`ePortfolio`](@ref), [`Problem`](@ref), [`Settings`](@ref)
"""
function eFrontier(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    (; E, V, u, d, N, eE, eV) = PS
    nL = lastindex(aCL)
    W = trues(nL)
    for k in eachindex(aCL)
        if norm(aCL[k].beta) < nS.tolNorm
            W[k] = false
        end
    end

    idW = findall(W)
    nP = length(idW)
    Ap = zeros(nP + 1, 3)
    mu = zeros(nP + 1)
    sgm = zeros(nP + 1)
    Z = zeros(nP + 1, N)
    z1 = zeros(N)

    for k in 1:nP
        t = aCL[idW[k]]
        S = @view t.S[1:N]

        F = (S .== IN)
        B = .!F
        U = (S .== UP)
        D = (S .== DN)

        z1[D] = d[D]
        z1[U] = u[U]
        z1[F] = t.alpha + t.beta * t.L1
        zB = @view z1[B]
        EF = @view E[F]
        uo = EF' * t.alpha + E[B]' * zB #μ_{o}=z_{o}′μ=α′μ_{I}+z_{B}′μ_{B}
        vo = t.alpha' * (V[F, F] * t.alpha + V[F, B] * zB * 2) + zB' * V[B, B] * zB  #v_{o}=α′V_{I}α+2α′c+z_{B}′V_{B}z_{B}

        #q = E[F]' * t.beta    #q=μ_{I}′β
        #a2 = 1.0 / q
        a2 = 1.0 / (EF' * t.beta)
        a1 = -2 * a2 * uo  #a₁=-2a₂μ_{o}
        a0 = vo + a2 * uo^2  #a₀=v_{o}+a₂μ_{o}²


        Ap[k, :] = [a2 a1 a0]
        mu[k] = z1' * E
        sgm[k] = sqrt(z1' * V * z1)
        Z[k, :] = z1

        if k == nP
            Ap[k+1, :] = [a2; a1; a0]
            z1[F] = t.alpha
            mu[k+1] = z1' * E
            sgm[k+1] = sqrt(z1' * V * z1)
            Z[k+1, :] = z1
        end

        #μ=μ_{o}+qL
        #v=z′Vz=v_{o}+L²q
    end
    if PS.equilibrate
        s1 = eV / eE
        s2 = s1 / eE    #eV/(eE^2)
        mu *= eE
        sgm *= sqrt(eV)
        Ap[:, 1] *= s2
        Ap[:, 2] *= s1
        Ap[:, 3] *= eV
    end
    return sEF(Ap, mu, sgm, Z, [idW; idW[end]])
end




"""

        z = ePortfolio(aEF::sEF, mu)
        z = ePortfolio(P::Problem, mu; nS=Settings(P))                                  #one-step, given mu
        z = ePortfolio(aCL::Vector{sCL}, P::Problem, L)
        z = ePortfolio(P::Problem; nS=Settings(P), settings, settingsLP, L::T=0.0)      #one-step, given L

`ePortfolio(aEF, mu)` compute the efficient portfolio given mu when Efficient Frontier is known, returns portfolio weights z (Nx1 vector of NaN if mu is out of range)

`ePortfolio(P, mu; nS=Settings(P))` compute the efficient portfolio for `P` given mu, returns portfolio weights z (Nx1 vector of NaN if mu is out of range)
If mu is a vector, z is a matrix, its row k is the portfolio weights for mu[k]. Here you can not set `init::Function`

`ePortfolio(aCL, P, L)` compute the efficient portfolio for `P` with known Critical Lines `aCL`, given `L`

`ePortfolio(P::Problem; nS=Settings(P), settings, settingsLP, L::T=0.0)` compute the efficient portfolio for `P` given `L`

See also [`EfficientFrontier.ECL`](@ref), [`eFrontier`](@ref), [`Problem`](@ref), [`Settings`](@ref)

"""
function ePortfolio(aEF::sEF, mu::Float64)
    u = aEF.mu
    Z = aEF.Z
    z = zeros(size(Z, 2))
    if mu > u[1] || mu < u[end]
        fill!(z, NaN)
        return z
    end
    if mu == u[end]
        return Z[end, :]
    end
    if mu == u[1]
        return Z[1, :]
    end
    k = searchsortedlast(u, mu, rev=true)
    bt = (mu - u[k+1]) / (u[k] - u[k+1])
    return bt * Z[k, :] + (1 - bt) * Z[k+1, :]
end

function ePortfolio(aCL::Vector{sCL{T}}, P::Problem{T}, L::T) where {T}
    (; u, d, N) = P
    k = 1
    #L = L < 0 ? 0 : L
    if L < 0
        L = 0
        @warn "L < 0, reset L to zero, L = 0"
    end
    while aCL[k].L0 > L
        k += 1
    end
    t = aCL[k]
    S = @view t.S[1:N]
    z = zeros(N)
    F = (S .== IN)
    U = (S .== UP)
    D = (S .== DN)

    z[D] = d[D]
    z[U] = u[U]
    z[F] = t.alpha + t.beta * L
    return z
end

function ePortfolio(P::Problem{T}; nS=Settings(P), settings=SettingsQP(PS), settingsLP=SettingsLP(PS), L::T=0.0) where {T}
    aCL = EfficientFrontier.ECL(P; numSettings=nS, settings=settings, settingsLP=settingsLP)
    return ePortfolio(aCL, P, L)
end

function ePortfolio(P::Problem, mu::Float64; nS=Settings(P))
    aEF = eFrontier(ECL(P; numSettings=nS), P; nS=nS)
    return ePortfolio(aEF, mu)
end

function ePortfolio(P::Problem, mu::Vector{Float64}; nS=Settings(P))
    aEF = eFrontier(ECL(P; numSettings=nS), P; nS=nS)
    M = length(mu)
    Z = similar(mu, M, P.N)
    for k in 1:M
        Z[k, :] = ePortfolio(aEF, mu[k])
    end
    return Z
end




"""

        x         = fPortfolio(P::Problem; settingsLP=settingsLP, L=0)              # active-set numerical solver
        x, status = fPortfolio(P::Problem, mu; settings=SettingsQP, check=true)     # `LightenQP`'s numerical solver

compute frontier portfolio at given `L` or `mu`.

    L=-Inf          :FP(L=-Inf), LMFP (Lowest Mean Frontier Portfolio)
    L=+Inf          :FP(L=+Inf), HMFP (Highest Mean Frontier Portfolio) HVEP (Highest Variance Efficient Portfolio)
    L=L0            :FP(L=L0), the frontier (minimum variance) portfolio at L=L0. L=0, LVEP (Lowest Variance Efficient Portfolio, also called GMVP, Global Minimum Variance Portfolio)
    mu=-Inf         :FP(L=-Inf), LMFP (Lowest Mean Frontier Portfolio)
    mu=+Inf         :FP(L=+Inf), HMFP (Highest Mean Frontier Portfolio) == HVEP (Highest Variance Efficient Portfolio)
    mu=mu0          :FP(mu=mu0), the frontier (minimum variance) portfolio at mu=mu0

if `check=false`, we do not check if mu is feasible or not (between lowest and highest mean)

See also [`ePortfolio`](@ref), [`Problem`](@ref), [`SettingsQP`](@ref), [`LightenQP.fPortfolio`](@ref)

"""
function fPortfolio(P::Problem{T}; settingsLP=SettingsLP(P), L::T=0.0) where {T}
    asQP(P; settingsLP=settingsLP, L=L)
end

#=
function fPortfolio(P::Problem{T}; settings=SettingsQP(P), L=0.0) where {T}
    Q = OOQP(P)
    fPortfolio(Q; settings=settings, L=L)
end =#

function fPortfolio(P::Problem{T}, mu::T; settings=SettingsQP(P), check=true) where {T}
    Q = OOQP(P)
    fPortfolio(Q, mu; settings=settings, check=check)
end




