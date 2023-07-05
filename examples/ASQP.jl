"Numerically stable methods for quadratic programming: an inertia-controlling active-set solver for general (definite/indefinite) dense quadratic programs"
module ASQP
#https://github.com/oxfordcontrol/GeneralQP.jl
#to do: 把 Settings 单独出来
#why no asLP? cause we need a LP to init ASQP

using LinearAlgebra
using Polynomials
using EfficientFrontier: EfficientFrontier, Problem, Status, Event, sCL, IN, DN, UP, OE, EO, computeCL!, Settings, SettingsLP
export solveASQP, asQP, asCL!, getSx

using StatusSwitchingQP.Simplex: cDantzigLP, maxImprvLP, stpEdgeLP


mutable struct UpdatableQR{T} <: Factorization{T}
    """
    Gives the qr factorization an (n, m) matrix as Q1*R1
    Q2 is such that Q := [Q1 Q2] is orthogonal and R is an (n, n) matrix where R1 "views into".
    """
    Q::Matrix{T}
    R::Matrix{T}
    n::Int
    m::Int

    Q1::SubArray{T,2,Matrix{T},Tuple{Base.Slice{Base.OneTo{Int}},UnitRange{Int}},true}
    Q2::SubArray{T,2,Matrix{T},Tuple{Base.Slice{Base.OneTo{Int}},UnitRange{Int}},true}
    R1::UpperTriangular{T,SubArray{T,2,Matrix{T},Tuple{UnitRange{Int},UnitRange{Int}},false}}

    function UpdatableQR(A::AbstractMatrix{T}) where {T}
        n, m = size(A)
        @assert m <= n "Too many columns in the matrix."

        F = qr(A)
        Q = F.Q * Matrix(I, n, n)
        R = zeros(T, n, n)
        R[1:m, 1:m] .= F.R

        new{T}(Q, R, n, m,
            view(Q, :, 1:m), view(Q, :, m+1:n),
            UpperTriangular(view(R, 1:m, 1:m)))
    end
end

function addColumn!(F::UpdatableQR{T}, a::AbstractVector{T}) where {T}
    a1 = F.Q1' * a
    a2 = F.Q2' * a

    x = copy(a2)
    for i = length(x):-1:2
        G, r = givens(x[i-1], x[i], i - 1, i)
        lmul!(G, x)
        lmul!(G, F.Q2')
    end

    F.R[1:F.m, F.m+1] .= a1
    F.R[F.m+1, F.m+1] = x[1]

    F.m += 1
    updateViews!(F)

    return a2
end

function addColHouseholder!(F::UpdatableQR{T}, a::AbstractVector{T}) where {T}
    a1 = F.Q1' * a
    a2 = F.Q2' * a

    Z = qr(a2)
    LAPACK.gemqrt!('R', 'N', Z.factors, Z.T, F.Q2) # Q2 .= Q2*F.Q
    F.R[1:F.m, F.m+1] .= a1
    F.R[F.m+1, F.m+1] = Z.factors[1, 1]
    F.m += 1
    updateViews!(F)

    return Z
end

function removeColumn!(F::UpdatableQR{T}, idx::Int) where {T}
    Q12 = view(F.Q, :, idx:F.m)
    R12 = view(F.R, idx:F.m, idx+1:F.m)

    for i in 1:size(R12, 1)-1
        G, r = givens(R12[i, i], R12[i+1, i], i, i + 1)
        lmul!(G, R12)
        rmul!(Q12, G')
    end

    for i in 1:F.m, j in idx:F.m-1
        F.R[i, j] = F.R[i, j+1]
    end
    F.R[:, F.m] .= zero(T)

    F.m -= 1
    updateViews!(F)

    return nothing
end

function updateViews!(F::UpdatableQR{T}) where {T}
    F.R1 = UpperTriangular(view(F.R, 1:F.m, 1:F.m))
    F.Q1 = view(F.Q, :, 1:F.m)
    F.Q2 = view(F.Q, :, F.m+1:F.n)
end

mutable struct NullspaceHessianLDL{T}
    # Struct for the LDL factorization of the reduced Hessian
    # on the nullspace Z of the working constraints
    # U'*D*U = Z[:, end:-1:1]'*P*Z[:, end:-1:1]
    # The projected hessian can include at most one negative eigenvector

    n::Int # Dimension of the original space
    m::Int # Dimension of the nullspace
    artificial_constraints::Int

    P::Matrix{T}     # The full hessian
    # Z is the nullspace, i.e. QR.Q2 where QR is defined below
    Z::SubArray{T,2,Matrix{T},Tuple{Base.Slice{Base.OneTo{Int}},UnitRange{Int}},true}
    U::UpperTriangular{T,SubArray{T,2,Matrix{T},Tuple{UnitRange{Int},UnitRange{Int}},false}}
    D::Diagonal{T,SubArray{T,1,Vector{T},Tuple{UnitRange{Int}},true}}
    QR::UpdatableQR{T}
    data::Matrix{T}  # That's where U is viewing into
    d::Vector{T}     # That's where D.diag views into
    indefinite_tolerance::T

    function NullspaceHessianLDL(P::Matrix{T}, A::AbstractMatrix{T}) where {T}
        @assert size(A, 1) == size(P, 1) == size(P, 2) "Matrix dimensions do not match."

        F = UpdatableQR(A)
        n = F.n
        m = F.n - F.m

        data = zeros(T, n, n)
        F.Q2 .= F.Q2[:, m:-1:1]
        WPW = F.Q2' * P * F.Q2
        WPW .= (WPW .+ WPW') ./ 2
        indefinite_tolerance = 1e-12

        C = cholesky(WPW, Val(true); check=false)
        if all(diag(C.U) .> indefinite_tolerance)
            F.Q2 .= F.Q2[:, C.p]
            artificial_constraints = 0
        else
            idx = findfirst(diag(C.U) .<= indefinite_tolerance)
            F.Q2 .= [F.Q2[:, C.p[idx:end]] F.Q2[:, C.p[1:idx-1]]]
            artificial_constraints = m + 1 - idx
        end
        m -= artificial_constraints
        U = view(data, 1:m, 1:m)

        F.Q2[:, artificial_constraints+1:end] .= F.Q2[:, end:-1:artificial_constraints+1]
        Z = view(F.Q2, :, artificial_constraints+1:size(F.Q2, 2))
        U .= view(C.U, 1:m, 1:m)
        d = ones(T, n)
        D = Diagonal(view(d, 1:m))

        new{T}(n, m, artificial_constraints, P, Z, UpperTriangular(U), D, F, data, d, indefinite_tolerance)
    end

end

function addConstraint!(F::NullspaceHessianLDL{T}, a::AbstractVector{T}) where {T}
    a2 = addColumn!(F.QR, a)
    if F.m == 1 # Nothing to write
        F.m -= 1
        updateViews!(F)
        return nothing
    end

    l = length(a2)
    for i = l:-1:2
        if l - i + 2 <= F.m
            G, _ = givens(a2[i-1], a2[i], l - i + 1, l - i + 2)
            rmul!(F.U.data, G)
        end
        G, _ = givens(a2[i-1], a2[i], i - 1, i)
        lmul!(G, a2)
    end
    F.m -= 1
    updateViews!(F)
    lmul!(sqrt.(F.D), F.U.data)
    H2Triangular!(F.U.data)

    z = view(F.Z, :, 1)
    Pz = F.P * z
    if F.d[F.m+1] <= -10
        # Recompute the last column of U
        l = reverse!(F.Z' * Pz)
        F.U[:, end] = F.U' \ l
    end
    # Correct last element of U
    u1 = view(F.U, 1:F.m-1, F.m)
    d_new = dot(z, F.P * z) - dot(u1, u1)
    # Prevent F.U having zero columns
    F.U[end] = max(sqrt(abs(d_new)), F.indefinite_tolerance)

    # Scale matrices so that diag(F.U) = ones(..)
    F.D.diag .= one(T)
    F.D.diag[end] *= sign(d_new)

    return nothing
end

function removeConstraint!(F::NullspaceHessianLDL{T}, idx::Int) where {T}
    @assert F.m == 0 || F.D[end] > F.indefinite_tolerance
    "Constraints can be removed only when the reduced Hessian was already Positive Semidefinite."

    if F.artificial_constraints == 0
        removeColumn!(F.QR, idx)
    else
        idx != 0 && @warn "Ignoring non-zero index: Removing an artificial constraint."
        F.artificial_constraints -= 1
    end

    z = view(F.QR.Q2, :, F.artificial_constraints + 1)
    Pz = F.P * z
    u = F.D \ (F.U' \ reverse!(F.Z' * Pz))
    d_new = dot(z, Pz) - dot(u, F.D * u)

    F.m += 1
    updateViews!(F)

    F.U[1:end-1, end] .= u
    F.U[end, end] = one(T)
    F.D[end] = d_new

    return nothing
end

function updateViews!(F::NullspaceHessianLDL{T}) where {T}
    if size(F.U, 1) > F.m
        F.U.data[F.m+1, :] .= zero(T)
        F.U.data[:, F.m+1] .= zero(T)
    end
    F.U = UpperTriangular(view(F.data, 1:F.m, 1:F.m))
    F.D = Diagonal(view(F.d, 1:F.m))
    F.Z = view(F.QR.Q2, :, F.artificial_constraints+1:size(F.QR.Q2, 2))
end

function H2Triangular!(A::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    for i in 1:n-1
        G, _ = givens(A[i, i], A[i+1, i], i, i + 1)
        lmul!(G, A)
    end
    return A
end

mutable struct NullspaceHessian{T}
    n::Int # Dimension of the original space
    m::Int # Dimension of the nullspace

    P::Matrix{T}     # The full hessian
    # Z is the nullspace, i.e. QR.Q2 where QR is defined below
    Z::SubArray{T,2,Matrix{T},Tuple{Base.Slice{Base.OneTo{Int}},UnitRange{Int}},true}
    ZPZ::SubArray{T,2,Matrix{T},Tuple{UnitRange{Int},UnitRange{Int}},false}  # equal to Z'*P*Z
    QR::UpdatableQR{T}
    data::Matrix{T}  # That's where ZPZ is viewing into

    function NullspaceHessian{T}(P::Matrix{T}, A::Matrix{T}) where {T}
        @assert size(A, 1) == size(P, 1) == size(P, 2) "Matrix dimensions do not match."

        F = UpdatableQR(A)
        n = F.n
        m = F.n - F.m

        data = zeros(T, n, n)
        ZPZ = view(data, n-m+1:n, n-m+1:n)
        ZPZ .= F.Q2' * P * F.Q2
        ZPZ .= (ZPZ .+ ZPZ') ./ 2

        new{T}(n, m, P, F.Q2, ZPZ, F, data)
    end

    function NullspaceHessian{T}(P::Matrix{T}, F::UpdatableQR{T}) where {T}
        n = F.n
        m = F.n - F.m
        @assert n == size(P, 1) == size(P, 2) "Dimensions do not match."

        data = zeros(T, n, n)
        ZPZ = view(data, n-m+1:n, n-m+1:n)
        ZPZ .= F.Q2' * P * F.Q2
        ZPZ .= (ZPZ .+ ZPZ') ./ 2

        new{T}(n, m, P, F.Q2, ZPZ, F, data)
    end
end

function addConstraint!(H::NullspaceHessian{T}, a::Vector{T}) where {T}
    Z = addColHouseholder!(H.QR, a)

    LAPACK.gemqrt!('L', 'T', Z.factors, Z.T, H.ZPZ)
    LAPACK.gemqrt!('R', 'N', Z.factors, Z.T, H.ZPZ)
    # ToDo: Force symmetry? (i.e. H.ZPZ .= (H.ZPZ .+ H.ZPZ')./2)
    H.m -= 1
    updateViews!(H)
    # H.ZPZ .= (H.ZPZ .+ H.ZPZ')./2

    return nothing
end

function removeConstraint!(H::NullspaceHessian{T}, idx::Int) where {T}
    removeColumn!(H.QR, idx)
    H.m += 1
    updateViews!(H)

    Pz = H.P * view(H.Z, :, 1)  # ToDo: avoid memory allocation
    mul!(view(H.ZPZ, 1, :), H.Z', Pz)
    for i = 2:H.m
        H.ZPZ[i, 1] = H.ZPZ[1, i]
    end

    return nothing
end

function updateViews!(H::NullspaceHessian{T}) where {T}
    range = H.n-H.m+1:H.n
    H.ZPZ = view(H.data, range, range)
    H.Z = H.QR.Q2
end


function removeConstraint!(data, idx::Int)
    constraint_idx = data.workingSet[idx]
    ignoredConstraintsAddRow!(data, data.A[constraint_idx, :], data.b[constraint_idx])
    prepend!(data.ignoredSet, constraint_idx)
    deleteat!(data.workingSet, idx)
    removeConstraint!(data.F, idx)
    updateViews!(data) # Updates ignoredA, bIgnored views
end

function addConstraint!(data, idx::Int)
    constraint_idx = data.ignoredSet[idx]
    addConstraint!(data.F, data.A[constraint_idx, :])
    ignoredConstraintsRemoveRow!(data, idx)
    deleteat!(data.ignoredSet, idx)
    append!(data.workingSet, constraint_idx)
    updateViews!(data)
end

function ignoredConstraintsAddRow!(data, a::AbstractVector, b::AbstractFloat)
    l = length(data.ignoredSet)
    data.shuffledA[data.m-l, :] = a
    data.bShuffled[data.m-l] = b
end

function ignoredConstraintsRemoveRow!(data, idx::Int)
    l = length(data.ignoredSet)
    start = data.m - l + 1
    @inbounds for j in 1:data.n, i in start+idx-1:-1:start+1
        data.shuffledA[i, j] = data.shuffledA[i-1, j]
    end
    @inbounds for i in start+idx-1:-1:start+1
        data.bShuffled[i] = data.bShuffled[i-1]
    end
end

function updateViews!(data)
    range = data.m-length(data.ignoredSet)+1:data.m
    data.ignoredA = view(data.shuffledA, range, :)
    data.bIgnored = view(data.bShuffled, range)
end


mutable struct Data{T}
    """
    Data structure for the solution of
        min         ½x'Px + q'x
        subject to  Ax ≤ b
    with an active set algorithm.

    The most important element is F, which holds
    - F.QR:      an updatable QR factorization of the working constraints
    - F.Z:       a view on an orthonormal matrix spanning the nullspace of the working constraints
    - F.P:       the hessian of the problem
    - F.U & F.D: the ldl factorization of the projected hessian, i.e. F.U'*F.D*F.D = F.Z'*F.P*F.Z

    The Data structure also keeps matrices of the constraints not in the working set (ignoredA and bIgnored)
    stored in a continuous manner. This is done to increase use of BLAS.
    """
    x::Vector{T}

    n::Int
    m::Int

    F::NullspaceHessianLDL{T}
    q::Vector{T}
    A::Matrix{T}
    b::Vector{T}
    workingSet::Vector{Int}
    ignoredSet::Vector{Int}
    alpha::Vector{T}
    residual::T

    iter::Int
    done::Bool

    e::Vector{T} # Just an auxiliary vector used in computing step directions

    ignoredA::SubArray{T,2,Matrix{T},Tuple{UnitRange{Int},Base.Slice{Base.OneTo{Int}}},false}
    bIgnored::SubArray{T,1,Vector{T},Tuple{UnitRange{Int}},true}
    shuffledA::Matrix{T}
    bShuffled::Vector{T}
    tol::T
    rMin::T
    rMax::T

    function Data(P::Matrix{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T}, x::Vector{T};
        rMin::T=zero(T), rMax::T=Inf, tol=1e-9) where {T}
        if rMin < tol
            rMin = -one(T)
        end

        m, n = size(A)
        workingSet = zeros(Int, 0)
        @assert maximum(A * x - b) < tol "The initial point is infeasible!"
        @assert norm(x) - rMin > -tol "The initial point is infeasible!"
        @assert norm(x) - rMax < tol "The initial point is infeasible!"

        ignoredSet = setdiff(1:m, workingSet)

        F = NullspaceHessianLDL(P, Matrix(view(A, workingSet, :)'))
        if F.m == 0 # To many artificial constraints...
            removeConstraint!(F, 0)
        end
        shuffledA = zeros(T, m, n)
        l = length(ignoredSet)
        shuffledA[end-l+1:end, :] .= view(A, ignoredSet, :)
        bShuffled = zeros(T, m)
        bShuffled[end-l+1:end] .= view(b, ignoredSet)

        e = zeros(T, n)
        alpha = zeros(T, m)

        new{T}(x, n, m, F, q, A, b, workingSet, ignoredSet, alpha,
            NaN, 0, false, e,
            view(shuffledA, m-l+1:m, :),
            view(bShuffled, m-l+1:m),
            shuffledA, bShuffled,
            tol, rMin, rMax)
    end
end



"""

        solveASQP(P::Matrix{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T}, x::Vector{T}; maxIter::Int=77777, kwargs...) where T

for quadratic programming problems: the initial feasible point can be obtained by performing Phase-I Simplex on the polyhedron Ax ≤ b

```math
    min	(1/2)x′Px + x′q
    s.t.	Ax ≤ b ∈ R^{M}
```


See also [`asQP`](@ref)
"""
function solveASQP(P::Matrix{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T},
    x::Vector{T}; maxIter::Int=77777, kwargs...) where {T}

    data = Data(P, q, A, b, x; kwargs...)

    while !data.done && data.iter <= maxIter && norm(data.x) <= data.rMax - 1e-10 && norm(data.x) >= data.rMin + 1e-10
        iterate!(data)
    end

    #return x, data.iter, data.workingSet
    #return x, data.workingSet
    return x
end

function iterate!(data::Data{T}) where {T}
    direction, stepsize, new_constraints = calculateStep(data)
    data.x .+= stepsize * direction

    if !isempty(new_constraints)
        addConstraint!(data, new_constraints[1])
    end
    #ToDo: break next condition into simpler ones or a more representative flag
    if (isempty(new_constraints) || data.F.m == 0) && norm(data.x) <= data.rMax - 1e-10
        if data.F.artificial_constraints > 0
            removeConstraint!(data.F, 0)
        else
            idx = KKTcheck!(data)
            !data.done && removeConstraint!(data, idx)
        end
    end
    data.iter += 1
end

function calculateStep(data)
    gradient = data.F.P * data.x + data.q
    if data.F.D[end] >= data.F.indefinite_tolerance
        qw = data.F.Z' * (gradient)
        if norm(qw) <= 1e-10 # We are alread on the optimizer of the current subproblem
            return zeros(data.n), 0, []
        end
        # Gill & Murray (1978) calculate the direction/stepsize as:
        # direction = data.F.Z*reverse(data.F.U\e)
        # direction .*= -sign(dot(direction, gradient))
        # alpha_min = -dot(gradient, direction)/data.F.D[end]
        # This assumes that we start from a feasible vertex or a feasible stationary point.
        # This assumption is suitable e.g. when we use a Phase I algorithm to find a feasible point
        # and it allows for a faster/simpler calculation of the direction.
        # Nevertheless, for purposes of generality, we calculate the direction in a more general way.
        direction = -data.F.Z * reverse(data.F.U \ (data.F.D \ (data.F.U' \ reverse(qw))))
        alpha_min = 1
    else
        e = view(data.e, 1:data.F.m)
        e[end] = 1
        #=
        if norm(data.F.U[:, end]) <= 1e-11
            data.F.U[end] = 1e-11
        end
        =#
        direction = data.F.Z * reverse(data.F.U \ e)
        if dot(direction, gradient) >= 0
            direction .= -direction
        end
        e .= 0
        alpha_min = Inf
    end

    A_times_direction = data.ignoredA * direction
    ratios = abs.(data.bIgnored - data.ignoredA * data.x) ./ (A_times_direction)
    ratios[A_times_direction.<=1e-11] .= Inf

    alpha_constraint = minimum(ratios)
    idx = findlast(ratios .== alpha_constraint)
    # @show alpha_constraint, data.ignoredSet[idx]

    alpha = min(alpha_min, alpha_constraint)
    if alpha == Inf
        # variable "direction" holds the unbounded ray
        @info "detected the problem to be unbounded (unbounded ray found)."
    end

    alpha_max = Inf
    if isfinite(data.rMax)
        # Calculate the maximum allowable step alpha_max so that rMin ≤ norm(x) ≤ rMax
        roots_rmax = roots(Poly([norm(data.x)^2 - data.rMax^2, 2 * dot(direction, data.x), norm(direction)^2]))
        if data.rMin > 0
            roots_rmin = roots(Poly([norm(data.x)^2 - data.rMin^2, 2 * dot(direction, data.x), norm(direction)^2]))
            roots_all = [roots_rmin; roots_rmax]
        else
            roots_all = roots_rmax
        end
        # Discard complex and negative steps
        roots_all = real.(roots_all[isreal.(roots_all)])
        roots_all = roots_all[roots_all.>=0]
        if length(roots_all) > 0
            alpha_max = minimum(roots_all)
        end
    end
    stepsize = min(alpha, alpha_max)
    @assert isfinite(stepsize) "Set the keyword argument rMax to a finite value if you want to get bounded solutions."
    if alpha_constraint == stepsize
        new_constraints = [idx]
    else
        new_constraints = []
    end

    return direction, stepsize, new_constraints
end

function KKTcheck!(data)
    grad = data.F.P * data.x + data.q
    alpha = -data.F.QR.R1 \ data.F.QR.Q1' * grad
    data.alpha .= 0.0
    data.alpha[data.workingSet] .= alpha
    data.residual = norm(data.F.Z' * grad)
    # data.residual = norm(grad + data.A[data.workingSet, :]'*alpha)

    idx = NaN
    if all(alpha .>= -2^-37)
        #if all(alpha .>= -1e-8)
        data.done = true
    else
        data.done = false
        idx = argmin(alpha)
    end

    return idx
end


"""

        asQP(PS::Problem{T}, L::T=0.0; settingsLP=SettingsLP(PS)) where T       : L version
        asQP(mu::T, PS::Problem{T}; settingsLP=SettingsLP(PS)) where T          : mu version

QP numerical solver, to the portfolio selection problem define by `PS::Problem`

See [`Documentation for EfficientFrontier.jl`](https://github.com/PharosAbad/EfficientFrontier.jl/wiki)

See also [`Problem`](@ref), [`SettingsLP`](@ref), [`asCL!`](@ref)
"""
function asQP(PS::Problem{T}, L::T=0.0; settingsLP=SettingsLP(PS)) where {T}

    E = PS.E
    if isinf(L)
        min = L == Inf ? false : true
        #mu = getfield(SimplexLP(PS; settings=settingsLP, min=min), 4)
        mu = E'*getfield(SimplexLP(PS; settings=settingsLP, min=min), 1)
        return asQP(mu, PS; settingsLP=settingsLP)
    end

    #finite L
    #(; E, V, u, d, G, g, A, b, N, M, J) = PS
    (; V, u, d, G, g, A, b, N, M, J) = PS
    (; tol, rule) = settingsLP

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    elseif rule == :stpEdgeLP
        solveLP = stpEdgeLP
    end

    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    Ms = M + J  #convert Gz<=g to equality
    Ns = N + J
    N1 = Ms + Ns
    S = fill(DN, N1)
    B = collect(Ns .+ (1:Ms))
    S[B] .= IN

    As = [A zeros(T, M, J)
        G Matrix{T}(I, J, J)]
    invB = Matrix{T}(I, Ms, Ms)
    bs = [b; g]
    ds = [d; zeros(T, J)]
    us = [u; fill(Inf, J)]

    q = As * ds
    for j in 1:Ms
        invB[j, j] = bs[j] >= q[j] ? one(T) : -one(T)
    end
    #q = abs.(As * ds - bs)
    q = abs.(q - bs)
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #灯塔的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    iH, x, invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    f = sum(x[Ns+1:end])
    if f > tol
        error("feasible region is empty")
    end

    x = x[1:N]

    iu = findall(u .< Inf)
    Aq = [A; -A; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]]
    bq = [b; -b; g; -d; u[iu]]
    if L == 0
        qq = zeros(T, N)
    else
        qq = -L * E
    end
    return solveASQP(V, qq, Aq, bq, x)
end

function asQP(mu::T, PS::Problem{T}; settingsLP=SettingsLP(PS)) where {T}
    (; E, V, u, d, G, g, A, b, N, M, J) = PS
    (; tol, rule) = settingsLP

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    elseif rule == :stpEdgeLP
        solveLP = stpEdgeLP
    end

    if isinf(mu)
        min = mu == Inf ? false : true
        #mu = getfield(SimplexLP(PS; settings=settingsLP, min=min), 4)
        mu = E'*getfield(SimplexLP(PS; settings=settingsLP, min=min), 1)
    end

    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    Ms = M + 1 + J  #convert Gz<=g to equality
    Ns = N + J
    N1 = Ms + Ns
    S = fill(DN, N1)
    B = collect(Ns .+ (1:Ms))
    S[B] .= IN

    As = [A zeros(T, M, J)
        E' zeros(T, 1, J)
        G Matrix{T}(I, J, J)]
    invB = Matrix{T}(I, Ms, Ms)
    bs = [b; mu; g]
    ds = [d; zeros(T, J)]
    us = [u; fill(Inf, J)]

    q = As * ds
    for j in 1:Ms
        invB[j, j] = bs[j] >= q[j] ? one(T) : -one(T)
    end
    #q = abs.(As * ds - bs)
    q = abs.(q - bs)
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #灯塔的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    iH, x, invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    f = sum(x[Ns+1:end])
    if f > tol
        error("feasible region is empty")
    end
    x = x[1:N]
    iu = findall(u .< Inf)
    Aq = [A; -A; E'; -E'; G; -Matrix{T}(I, N, N); Matrix{T}(I, N, N)[iu, :]]
    bq = [b; -b; mu; -mu; g; -d; u[iu]]
    return solveASQP(V, zeros(T, N), Aq, bq, x)

end



"""

        asCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settingsLP=SettingsLP(PS), kwargs...) where T

compute the Critical Line Segments by an inertia-controlling active-set numerical solver, at LMEP. Save the CL to aCL if done

See also [`asQP`](@ref), [`sCL`](@ref), [`Problem`](@ref), [`Settings`](@ref), [`SettingsLP`](@ref)
"""
function asCL!(aCL::Vector{sCL{T}}, PS::Problem{T}; nS=Settings(PS), settingsLP=SettingsLP(PS), kwargs...) where {T}
    #x, aS = asQP(PS; settingsLP=settingsLP)   #GMVP
    #S = activeS(aS, PS) #if fully degenerated, using x to determine S may be more reliable

    x = asQP(PS; settingsLP=settingsLP)   #GMVP
    S = getSx(x, PS, nS)
    computeCL!(aCL, S, PS, nS)
end


"""

        getSx(x, P, nS)

extract the `Status` information from the portfolio weights `x`, given `P::Problem` and `nS::Settings`

See also [`Status`](@ref), [`Problem`](@ref), [`Settings`](@ref)
"""
function getSx(x, P, nS)
    #from  optimal solution directly. Hence, high accuracy needed
    (; u, d, G, g, N, J) = P
    tolS = nS.tolS

    S = fill(IN, N + J)
    Sv = @view S[1:N]
    Sv[abs.(x - d).<tolS] .= DN
    Sv[abs.(x - u).<tolS] .= UP

    for m = 1:J
        S[N+m] = abs(g[m] - G[m, :]' * x) < tolS ? EO : OE
    end
    return S
end

#=
function activeS(aS, P)
    #=
    not reliable, why?
    * for HMFP, it is a LP, we know there is degenerated case, even fully degenerated. that can belong to many working set
    * for corner portfolios, it is the joit of 2 CL, each has its working set
    =#
    (; u, M, N, J) = P
    S = fill(IN, N + J)
    S[(1:J).+N] .= OE
    k0 = 2 * M      #equalities
    kj = k0 + J     #inequalities
    kd = kj + N     #lower bounds
    iu = findall(u .< Inf)
    for k in eachindex(aS)
        n = aS[k]
        if n <= kj
            if J > 0 && n > k0
                S[N+n-k0] = EO
            end
            continue
        end

        if n <= kd  #z >= d
            S[n-kj] = DN
        else    #z <= u
            S[iu[n-kd]] = UP
        end
    end
    return S
end
=#

end

