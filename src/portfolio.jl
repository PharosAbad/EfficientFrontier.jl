
"""

        aEF = eFrontier(aCL::sCL, PS::Problem; tolNorm = 2^-26)

compute the Full Efficient Frontier by connecting Critical Line Segments

the return aEF has the follwing structure

        struct sEF    #Efficient Frontier       Float64 is OK
            Ap::Matrix{Float64}     # v=a₂μ²+a₁μ+a₀, each row [a₂ a₁ a₀]
            mu::Vector{Float64}     #higher mean
            sgm::Vector{Float64}    #higher sigma
            Z::Matrix{Float64}      #weights, each corner portfolio in one row
            ic::Vector{Int64}       #id of related critical line
        end

"""
function eFrontier(aCL, PS; tolNorm = 2^-26)
    (; E, V, u, d, N) = PS
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
        S = t.S[1:N]

        F = (S .== IN)
        B = .!F
        U = (S .== UP)
        D = (S .== DN)
        
        z1[D] = d[D]
        z1[U] = u[U]
        z1[F] = t.alpha + t.beta * t.L1
        zB = @view z1[B]
        uo = E[F]' * t.alpha + E[B]' * zB #μ_{o}=z_{o}′μ=α′μ_{I}+z_{B}′μ_{B}
        vo = t.alpha' * (V[F, F] * t.alpha + V[F, B] * zB * 2) + zB' * V[B, B] * zB  #v_{o}=α′V_{I}α+2α′c+z_{B}′V_{B}z_{B}
        
        #q = E[F]' * t.beta    #q=μ_{I}′β
        #a2 = 1.0 / q
        a2 = 1.0 / (E[F]' * t.beta)
        a1 = -2 * a2 * uo  #a₁=-2a₂μ_{o}
        a0 = vo + a2 * uo^2  #a₀=v_{o}+a₂μ_{o}²
        

        Ap[k, :] = [a2 a1 a0]
        mu[k] = z1' * E
        sgm[k] = sqrt(z1' * V * z1)
        Z[k, :] = z1

        if k == nP
            Ap[k+1, :] = [a2 a1 a0]
            z1[F] = t.alpha
            mu[k+1] = z1' * E
            sgm[k+1] = sqrt(z1' * V * z1)
            Z[k+1, :] = z1
        end

        #μ=μ_{o}+qL
        #v=z′Vz=v_{o}+L²q
        #push!(aEF, sEF(a2, a1, a0, z1'*E, sqrt(z1'*V*z1), z0'*E, sqrt(z0'*V*z0), copy(z1), copy(z0)))
    end
    return sEF(Ap, mu, sgm, Z, [idW; idW[end]])
end

#=
function CornerP()
    return sEF.Z
end


function CornerP()
    nP = lastindex(aCL)
    W = trues(nP)
    for k in eachindex(aCL)
        if norm(aCL[k].beta) < epsD
            W[k] = false
        end
    end

    eCL = aCL[W]
    nP = lastindex(eCL)
    Z = zeros(nP + 1, N)
    for k in 1:nP
        t = eCL[k]
        S = t.S[1:N]

        F = (S .== IN)
        D = (S .== DN)
        U = (S .== UP)

        Z[k, D] = d[D]
        Z[k, U] = u[U]
        Z[k, F] = t.alpha + t.beta * t.L1
        if k == nP
            Z[k+1, :] = Z[k, :]
            Z[k+1, F] = t.alpha
        end
    end
    return Z
end
=#

#=
function setup(E0, V0, A0, b0, d0, u0, G0, g0)
    #global aCL = aCL0
    global E = vec(E0)
    global N = length(E)
    global V = copy(V0)
    global A = copy(A0)
    global b = copy(b0)
    global d = copy(d0)
    global u = copy(u0)
    global G = copy(G0)
    global g = copy(g0)
    global M = length(b)
    global J = length(g)
end

function noShortsale(E0, V0)
    global E = vec(E0)
    global N = length(E)
    global V = copy(V0)
    global A = ones(Float64, 1, N)
    global b = ones(Float64,1)
    global d = zeros(Float64, N)
    global u = fill(Float64(Inf), N)
    global G = Matrix{Float64}(undef, 0, N)
    global g = Vector{Float64}(undef, 0)
    global M = length(b)
    global J = length(g)
end
=#