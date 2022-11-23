module EfficientFrontier
using LinearAlgebra, Clarabel, SparseArrays

#export Status, Event, sCL, aCL, computeCL, ECL, CornerP
export Status, Event, sCL, computeCL, IN, DN, UP, OE, EO


epsD::Float64 = 2^-26    #1.490116119384765625e-08

@enum Status begin
    IN
    DN  #down, lower bound
    UP  #upper bound
    OE  #original <= not binded
    EO  #edge, <= as =
end

struct Event
    From::Status
    To::Status
    id::Int
    L::Float64
end

struct sCL    #critical line segment
    S::Vector{Status}   # (N+J)x1
    K::Integer  #number of IN assets
    L1::Float64 #higher lambda
    I1::Vector{Event}   #go in/out events at L1
    L0::Float64 #lower lambda
    I0::Vector{Event}   #go in/out events at L0
    alpha::Vector{Float64}  # K x 1
    beta::Vector{Float64}   # K x 1
end

N::Int = 7
E = Vector{Float64}(undef, N)
V = Matrix{Float64}(undef, 0, N)
A = ones(Float64, 1, N)
b = ones(Float64, 1)
G = Matrix{Float64}(undef, 0, N)
g = Vector{Float64}(undef, 0)
d = zeros(Float64,N)
u = fill(Float64(Inf), N)
M = length(b)
J = length(g)
aCL = Vector{sCL}(undef, 0)

function Int(v::Vector{Status})
    x = zeros(Int, length(v))
    for k in eachindex(v)
        x[k] = Int(v[k])
    end
    return x
end


function ClarabelQP(E, V, mu, Ae, be, d, u, G, g)
    N = length(E)
    iu = u .!= Inf
    Nu = sum(iu)    
    P = sparse(V)
    q = zeros(N)
    A = sparse([E'; Ae; G; -Matrix{Float64}(I, N, N); Matrix{Float64}(I, N, N)[iu,:]])
    b = [mu; be; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(1 + length(be)), Clarabel.NonnegativeConeT(length(g) + N + Nu)]
    settings = Clarabel.Settings()
    settings.verbose = false
    solver = Clarabel.Solver()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    Clarabel.solve!(solver)
end



function ClarabelLP(E, Ae, be, d, u, G, g)
    N = length(E)
    iu = u .!= Inf
    Nu = sum(iu)
    P = sparse(zeros(Float64, N, N))
    q = -E
    A = sparse([Ae; G; -Matrix{Float64}(I, N, N); Matrix{Float64}(I, N, N)[iu,:]])
    b = [be; g; -d; u[iu]]
    cones = [Clarabel.ZeroConeT(length(be)), Clarabel.NonnegativeConeT(length(g) + N + Nu)]
    settings = Clarabel.Settings()
    settings.verbose = false
    solver = Clarabel.Solver()
    Clarabel.setup!(solver, P, q, A, b, cones, settings)
    Clarabel.solve!(solver)
end

function getRows(A)
    #indicate the non-redundant rows, the 1st row should be non-zeros (vector one here)
    H = @view(A[1:1, :])
    R = falses(size(A, 1))
    R[1] = true
    for r in axes(A, 1)[(begin+1):end]
        v = @view(A[r, :])
        if norm(v - H' * (H' \ v)) > epsD
            R[r] = true
            H = @view(A[R, :])
        end
    end
    return findall(R)
end


function initCL()
    d[isinf.(d)] .= -1.0    #replace -Inf in lower bound by -1.0
    iu = u .< d
    u[iu] = d[iu]   #make sure u >= d

    x = ClarabelLP(E, A, b, d, u, G, g)
    if Int(x.status) != 1   #SOLVED
        error("Not able to find the max expected return of efficient frontier")
    end
    #display(-x.obj_val)
    #y = ClarabelQP(E, V, -x.obj_val, A, b, d, u, G, g) #maybe nothing IN
    #y = ClarabelQP(E, V, -x.obj_val*0.999, A, b, d, u, G, g)    
    y = ClarabelQP(E, V, x.obj_val * (2^-14 - 1), A, b, d, u, G, g)
    if Int(y.status) != 1   #SOLVED
        error("Not able to find the max expected return efficient portfolio")
    end
    Y = y.s
    #display(Y')
    iu = u .!= Inf
    Nu = sum(iu)
    D = Y[(1:N).+(M+J+1)] .< epsD
    U = falses(N)
    U[iu] = Y[(1:Nu).+(M+J+N+1)] .< epsD
    #F = .!(U .| D)
    S = fill(IN, N)
    S[D] .= DN
    S[U] .= UP
    Se = fill(EO, J)
    Og = Y[(1:J).+(M+1)] .> epsD
    Se[Og] .= OE
    Sj = [S; Se]

    # to do 直到 IN 不为空，然后 往上追溯到 I1 为空 2022-11-21 09:46:41

    #display(Int(Status[IN, UP, DN, IN, DN, DN, UP, UP, DN, DN, DN, DN, UP, UP, OE, OE])' - Int(Sj)')
    #error("Debug ...")
    #computeCL(Sj, true)
    computeCL(Sj)

end

#function computeCL(Sj, snglOK=false)
function computeCL(Sj::Vector{Status})
    S = Sj[1:N]
    F = (S .== IN)
    if sum(F) < 1   #IN is empty
        return false
    end
    B = .!F
    D = (S .== DN)
    U = (S .== UP)
    Se = Sj[(N.+(1:J))]
    Eg = (Se .== EO)
    Og = (Se .== OE)
    GE = @view(G[Eg, :])
    AE = vcat(A[:, F], GE[:, F])
    idAE = vcat(axes(A, 1), axes(GE, 1)) # id of each constraint

    zB = d[B]
    zB[U[B]] = u[U]
    AB = vcat(A[:, B], GE[:, B])
    bE = vcat(b, g[Eg]) - AB * zB
    rb = getRows([AE bE])
    if length(getRows(AE)) != length(rb)
        return false    #infeasible
    else
        AE = AE[rb, :]
        bE = bE[rb]
        AB = AB[rb, :]
    end

    EF = E[F]
    VF = V[F, F]
    iV = inv(VF)
    K = length(EF)

    c = V[F, B] * zB  #c=V_{IB}z_{B}
    T = iV * AE'   #T=V_{I}⁻¹A_{I}′
    C = inv(AE * T)   #C=(A_{I}T)⁻¹
    TC = T * C
    VQ = iV - T * TC'    #Q=I-A_{I}′CT′   V_{I}⁻¹Q=V_{I}⁻¹-TCT′

    alpha = TC * bE - VQ * c    #α=TCb_{E}-V_{I}⁻¹Qc
    beta = VQ * EF    #β=V_{I}⁻¹Qμ_{I}

    #=
    if !snglOK && norm(beta) < epsD    #singula
        return false
    end
    =#

    #α_{λ}=-C(T′c+b_{E})
    alphaL = -(TC' * c + C * bE)
    #β_{λ}=CT′μ_{I}
    betaL = TC' * EF

    #η_{B}==V_{BI}α+V_{B}z_{B}+A_{B}′α_{λ}+L(V_{BI}β-μ_{B}+A_{B}′β_{λ})
    gamma = V[B, F] * alpha + V[B, B] * zB + AB' * alphaL
    delta = V[B, F] * beta - E[B] + AB' * betaL
    #display(delta)

    LE0 = Vector{Event}(undef, 0)
    LE1 = Vector{Event}(undef, 0)
    ik = findall(F)
    for k in eachindex(alpha)
        j = ik[k]
        t = beta[k]
        h = alpha[k]
        dL = (d[j] - h) / t
        uL = (u[j] - h) / t
        if t > epsD
            push!(LE0, Event(IN, DN, j, dL))
            if u[j] < Inf
                push!(LE1, Event(UP, IN, j, uL))
            end
        elseif t < -epsD
            push!(LE1, Event(DN, IN, j, dL))
            if u[j] < Inf
                push!(LE0, Event(IN, UP, j, uL))
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
        if t > epsD
            if D[j]
                push!(LE0, Event(DN, IN, j, dL))
            else
                push!(LE1, Event(IN, UP, j, dL))
            end
        elseif t < -epsD
            if D[j]
                push!(LE1, Event(IN, DN, j, dL))
            else
                push!(LE0, Event(UP, IN, j, dL))
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
        if t > epsD
            push!(LE0, Event(OE, EO, j, dL))
        elseif t < -epsD
            push!(LE1, Event(EO, OE, j, dL))
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
            if t > epsD
                push!(LE0, Event(EO, OE, j, dL))
            elseif t < -epsD
                push!(LE1, Event(OE, EO, j, dL))
            else
                if h <= 0
                    return false
                end
            end
        end
    end

    L0 = -Inf
    I0 = Vector{Event}(undef, 0)
    nL0 = length(LE0)
    if nL0 > 0
        LE0 = sort(LE0, by=x -> x.L, rev=true)
        ic0 = LE0[1]
        L0 = ic0.L
        I0 = [ic0]
    end

    if L0 < -epsD
        L0 = 0.0
    end

    L1 = Inf
    I1 = Vector{Event}(undef, 0)
    nL1 = length(LE1)
    if nL1 > 0
        LE1 = sort(LE1, by=x -> x.L)
        ic1 = LE1[1]
        L1 = ic1.L
        I1 = [ic1]
    end
    if L1 < L0 + epsD
        return false
    end

    if nL0 > 1 && L0 - LE0[2].L < epsD
        I0 = [I0; LE0[2]]
        for k in 3:nL0
            if L0 - LE0[k].L < epsD
                I0 = [I0; LE0[k]]
            else
                break
            end
        end
    end

    if nL1 > 1 && LE1[2].L - L1 < epsD
        I1 = [I1; LE1[2]]
        for k in 3:nL1
            if LE1[k].L - L1 < epsD
                I1 = [I1; LE1[k]]
            else
                break
            end
        end
    end

    push!(aCL, sCL(copy(Sj), K, L1, I1, L0, I0, alpha, beta))
    return true
end

function ECL()
    global aCL = Vector{sCL}(undef, 0)
    if !initCL()
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
                id += N
            end
            S[id] = To
        end

        #display(Int(S)')
        #if computeCL(S, true)
        if computeCL(S)
            t = aCL[end]
            if t.L0 <= 0
                break
            end
        else
            display(Int(S)')
            error("Critical Line Not Fineshed")
        end

        #display(t)
        #error("Debug ...")
    end
    #return true
    return aCL
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

end
