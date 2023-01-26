
"""

        aEF = eFrontier(aCL::Vector{sCL{T}}, PS::Problem{T}; tolNorm = 2^-26) where T

compute the Full Efficient Frontier by Status-Segment Method

the return aEF has the follwing structure

        struct sEF    #Efficient Frontier       Float64 is OK
            Ap::Matrix{Float64}     # v=a₂μ²+a₁μ+a₀, each row [a₂ a₁ a₀]
            mu::Vector{Float64}     #higher mean
            sgm::Vector{Float64}    #higher sigma
            Z::Matrix{Float64}      #weights, each corner portfolio in one row
            ic::Vector{Int64}       #id of related critical line
        end

"""
function eFrontier(aCL::Vector{sCL{T}}, PS::Problem{T}; tolNorm=2^-26) where {T}
    (; E, V, u, d, N, eE, eV) = PS
    nL = lastindex(aCL)
    W = trues(nL)
    for k in eachindex(aCL)
        if norm(aCL[k].beta) < tolNorm
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

compute the efficient portfolio given mu, returns portfolio weights z (Nx1 vector of NaN if mu is out of range)

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


"""

        z = ePortfolio(P::Problem, mu)

compute the efficient portfolio given mu, returns portfolio weights z (Nx1 vector of NaN if mu is out of range)
If mu is a vector, z is a matrix, its row k is the portfolio weights for mu[k]

"""
function ePortfolio(P::Problem, mu::Float64)
    aEF = eFrontier(ECL(P), P)
    return ePortfolio(mu, aEF)
end

function ePortfolio(mu::Vector{Float64}, P::Problem)
    aEF = eFrontier(ECL(P), P)
    J = length(mu)
    Z = similar(mu, J, P.N)
    for k in 1:J
        Z[k,:] = ePortfolio(mu[k], aEF)
    end
    return Z
end