
"""

        aEF = eFrontier(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where T

compute the Full Efficient Frontier by Status-Segment Method

the return aEF has the follwing structure

        struct sEF    #Efficient Frontier       Float64 is OK
            Ap::Matrix{Float64}     # v=a₂μ²+a₁μ+a₀, each row [a₂ a₁ a₀]
            mu::Vector{Float64}     #higher mean
            sgm::Vector{Float64}    #higher sigma
            Z::Matrix{Float64}      #weights, each corner portfolio in one row
            ic::Vector{Int}         #id of related critical line
        end

See also [`EfficientFrontier.ECL`](@ref), [`ePortfolio`](@ref), [`Problem`](@ref), [`Settings`](@ref)
"""
function eFrontier(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS)) where {T}
    (; E, V, u, d, N, eE, eV) = PS
    (; tolNorm, tolL) = nS
    nL = lastindex(aCL)
    W = trues(nL)
    @inbounds for k in eachindex(aCL)
        t = aCL[k]
        if norm(t.beta) < tolNorm || t.L1 - t.L0 < tolL  #removable
            W[k] = false
        end
    end

    idW = findall(W)
    nP = length(idW)
    Ap = zeros(T, nP + 1, 3)
    mu = zeros(T, nP + 1)
    sgm = zeros(T, nP + 1)
    Z = zeros(T, nP + 1, N)
    z1 = zeros(T, N)

    @inbounds for k in 1:nP
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
    return sEF(convert(Matrix{Float64}, Ap), convert(Vector{Float64}, mu), convert(Vector{Float64}, sgm), convert(Matrix{Float64}, Z), [idW; idW[end]])
end




"""

        z = ePortfolio(mu, aEF::sEF)
        z = ePortfolio(mu, P::Problem; nS=Settings(P), settings, settingsLP)            #one-step, given mu
        z = ePortfolio(P::Problem, L, aCL::Vector{sCL})
        z = ePortfolio(P::Problem, L::T=0.0; nS=Settings(P), settings, settingsLP)      #one-step, given L

`ePortfolio(mu::Float64, aEF)` compute the efficient portfolio given mu when Efficient Frontier is known, returns portfolio weights z (Nx1 vector of NaN if mu is out of range)

`ePortfolio(mu, P; nS=Settings(P))` compute the efficient portfolio for `P` given mu::Float64, returns portfolio weights z (Nx1 vector of NaN if mu is out of range)
If mu is a ::Vector{Float64}, z is a matrix, its row k is the portfolio weights for mu[k]. Here you can not set `init::Function`

`ePortfolio(P, L::T, aCL)` compute the efficient portfolio for `P` with known Critical Lines `aCL`, given `L::T` (Float64, BigFloat)

`ePortfolio(P::Problem, L::T=0.0; nS=Settings(P), settings, settingsLP)` compute the efficient portfolio for `P` given `L`

See also [`EfficientFrontier.ECL`](@ref), [`eFrontier`](@ref), [`Problem`](@ref), [`Settings`](@ref)

"""
function ePortfolio(mu::Float64, aEF::sEF)
    u = aEF.mu
    Z = aEF.Z
    z = zeros(Float64, size(Z, 2))
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


function ePortfolio(mu::Float64, P::Problem; nS=Settings(P), settings=SettingsQP(P), settingsLP=SettingsLP(P))
    aEF = eFrontier(ECL(P; numSettings=nS, settings=settings, settingsLP=settingsLP), P; nS=nS)
    return ePortfolio(mu, aEF)
end

function ePortfolio(mu::Vector{Float64}, P::Problem; nS=Settings(P), settings=SettingsQP(P), settingsLP=SettingsLP(P))
    aEF = eFrontier(ECL(P; numSettings=nS, settings=settings, settingsLP=settingsLP), P; nS=nS)
    M = length(mu)
    Z = similar(mu, M, P.N)
    for k in 1:M
        Z[k, :] = ePortfolio(mu[k], aEF)
    end
    return Z
end


function ePortfolio(P::Problem{T}, L::T, aCL::Vector{sCL{T}}) where {T}
    (; u, d, N) = P
    k = 1
    #L = L < 0 ? 0 : L
    if L < 0
        L = 0
        @warn "L < 0, reset L to zero, L = 0"
    end

    L1 = aCL[1].L1
    if L > L1
        L = L1
    end

    while aCL[k].L0 > L
        k += 1
    end
    t = aCL[k]
    Kb = 1
    if k > 1 && L > t.L1
        Sb = nextS(t, true, N)
        Kb = sum(Sb .== IN)
    end

    if Kb == 0   # degenerate
        S = @view Sb[1:N]
    else
        S = @view t.S[1:N]
    end
    #S = @view t.S[1:N]
    z = zeros(T, N)
    U = (S .== UP)
    D = (S .== DN)

    z[D] = d[D]
    z[U] = u[U]
    if Kb > 0
        F = (S .== IN)
        z[F] = t.alpha + t.beta * L
    end
    return z
end

function ePortfolio(P::Problem{T}, L::T=0.0; nS=Settings(P), settings=SettingsQP(P), settingsLP=SettingsLP(P)) where {T}
    aCL = ECL(P; numSettings=nS, settings=settings, settingsLP=settingsLP)
    return ePortfolio(P, L, aCL)
end




"""

            mu2L(aEF::sEF, aCL::Vector{sCL{T}}, mu::T) where T

compute L from mu
"""
function mu2L(aEF::sEF, aCL::Vector{sCL{T}}, mu::T) where {T}
    #L=a₂(μ-μ_{o}) , thus we use the linearity L-L₀=β(L₁-L₀)   β=((μ-μ₀)/(μ₁-μ₀))
    u = aEF.mu
    if mu >= u[1]
        return aCL[aEF.ic[1]].L1
    elseif mu <= u[end]
        return 0.0
    end

    #k = searchsortedlast(u, mu, rev=true)
    k = 1
    while u[k] > mu
        k += 1
    end
    k -= 1
    t = aCL[aEF.ic[k]]  #  u[k] > mu >= u[k+1]
    bt = (mu - u[k+1]) / (u[k] - u[k+1])
    bt * (t.L1 - t.L0) + t.L0
end


"""

        L2mu(aEF::sEF, aCL::Vector{sCL{T}}, mu::T) where T

compute mu from L
"""
function L2mu(aEF::sEF, aCL::Vector{sCL{T}}, L::T) where {T}
    #μ=μ_{o}+qL => μ-μ₀=β(μ₁-μ₀)   β=((L-L₀)/(L₁-L₀))
    u = aEF.mu
    if L <= 0
        return u[end]
    end
    k = 1
    while aCL[aEF.ic[k]].L0 > L
        k += 1
    end
    t = aCL[aEF.ic[k]]  #L >= t.L0
    if L > t.L1     #degenerated or singular points
        L = t.L1
    end
    bt = (L - t.L0) / (t.L1 - t.L0)
    bt * (u[k] - u[k+1]) + u[k+1]

end