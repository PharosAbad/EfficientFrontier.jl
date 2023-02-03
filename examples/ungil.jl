#source https://github.com/ungil/Markowitz.jl/blob/master/examples/frontier.jl
#Example: # lower and upper bounds, 2 inequality constraints, 2 equality constraints

using LinearAlgebra
using EfficientFrontier

if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    include("./Markowitz.jl")
    using .Markowitz
end

function main()

    assets = ["Bonds - US Government"
        "Bonds - US Corporate"
        "Bonds - International"
        "Bonds - High Yield"
        "Bonds - Bank Loans"
        "Bonds - Emerging USD Debt"
        "Bonds - Emerging Local Debt"
        "Alternative - Emerging FX"
        "Alternative - Commodities"
        "Alternative - REITs"
        "Stocks - US Large"
        "Stocks - US Small"
        "Stocks - International"
        "Stocks - Emerging"]

    E, V = EfficientFrontier.EVdata(:Ungil, false)
    println("\n--- Status-Segment Method vs Markowitz's CLA  ---\n")

    class = vec([:FI :FI :FI :FI :FI :FI :FI :ALT :ALT :ALT :EQ :EQ :EQ :EQ])
    m = markowitz(E, V, names=assets, # asset bounds by class: stocks -10/30, bonds 0/20, alt. 0/10
        lower=-0.1 * (class .== :EQ),
        upper=0.3 * (class .== :EQ) + 0.2 * (class .== :FI) + 0.1 * (class .== :ALT))
    unit_sum(m)
    add_constraint(m, 1 * (class .== :EQ), '>', 0.3) # net equity exposure between 30% and 60%
    add_constraint(m, 1 * (class .== :EQ), '<', 0.6)
    add_constraint(m, [1 1 0 0 0 0 0 0 0 0 0 0 0 0], '=', 0.25) # US govt + Investment Grade = 25%

    ts = @elapsed f = frontier(m)
    println("Markowitz CLA:  ", ts, "  seconds") #0.001356236  seconds
    display(f.weights)



    #EfficientFrontier
    A = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
        1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    b = [1.0; 0.25]
    G = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0]
    G[1, :] = -G[1, :]
    g = [-0.3; 0.6]
    d = vec([-0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.1 -0.1 -0.1 -0.1])
    u = vec([0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.3 0.3 0.3 0.3])


    P = Problem(E, V, u, d, G, g, A, b)
    ts = @elapsed aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)
    #display(aCL)    #there are 2  singular CL (beta is a zero vector, a line in R^(N+J+M+1) space, but a single point in mean-variance space )
    println("Status-Segment Method:  ", ts, "  seconds") #0.000833465  seconds
    display(aEF.Z)  # Kinks on the frontier due to singular CL

    Pt = Problem(E, V, u, d, G, g, A, b; equilibrate=true)
    ts = @elapsed aCLt = EfficientFrontier.ECL(Pt)
    aEFt = eFrontier(aCLt, Pt)
    println("equilibrate:  ", ts, "  seconds") #0.000916927  seconds


    #BigFloat
    Pb = Problem(convert(Vector{BigFloat}, vec(E)), V, u, d, G, g, A, b)
    ts = @elapsed aCLb = EfficientFrontier.ECL(Pb)
    aEFb = eFrontier(aCLb, Pb)
    println("BigFloat:  ", ts, "  seconds")   #0.006389718  seconds



    #println("improvements  ", round.([maximum(abs.(aEFb.Z-aEF.Z)), maximum(abs.(aEFt.Z-aEF.Z))], sigdigits=3))
    #println("BigFloat over Float64+equilibrate, improvements: ", round(maximum(abs.(aEFb.Z-aEFt.Z)), sigdigits=3))
    println("BigFloat over Float64+equilibrate: ", round(norm(aEFb.Z - aEFt.Z), sigdigits=3))
    println("Float64+equilibrate over Float64: ", round(norm(aEFt.Z - aEF.Z), sigdigits=3))
    println("BigFloat over Float64: ", round(norm(aEFb.Z - aEF.Z), sigdigits=3))
    #maximum(abs.(f.weights[end-15:end,:]-aEF.Z[end-15:end,:]))

    println("
    v0.3  `Float64` is enough thus the default
        `equilibrate` may improve precision if range of non-zeros in abs.(E) or abs.(V) is too large
        `BigFloat` can be tens (30 to 70), even more than 200 times slow.
        If precision is a problem, `Float64+equilibrate` is recommended. It is better (fast and stable) most of the time, but not always. (Due to the truncation of floating number, 
        e.g., 1/60*50-50/60 = 0, but if scaled down the data by 10, 1/6*5-5/6 = 1.1102230246251565e-16)
    ")

end


main()
nothing
