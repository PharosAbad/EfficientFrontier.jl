#Example: CLA (Critical Line Algorithm) fail, for it only handle either one IN or one OUT


using EfficientFrontier

if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    include("./Markowitz.jl")
    using .Markowitz
end


function main()

    println("\n--- Status-Segment Method vs Markowitz's CLA  ---\n")

    E, V = EfficientFrontier.EVdata(:Ungil, false)
    E = [7.2 0.7 0.8 2.3 2.2 1.9 5.6 5.6 7.2 1.3 0.7 -0.1 4.1 7.2] # 3 assets share the highest expected return

    m = markowitz(E, V)
    unit_sum(m) # total weight = 100%
    ts = @elapsed f = frontier(m)   #Warning: tweaking mu[8] mu[13] to ensure the solution is unique  => perturbed method
    println("Markowitz CLA:  ", ts, "  seconds")    #0.00059 seconds
    display(f.weights)


    P = Problem(E, V; equilibrate=false)
    ts = @elapsed aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)
    println("Status-Segment Method:  ", ts, "  seconds")   #0.00077 seconds
    display(aEF.Z)

    println("
    the first corner portfolio has THREE assets, not ONE asset.   CLA's perturbed method fail to work
    the first four corner portfolios by Markowitz's CLA are wrong, the corresponding critical lines violate the KKT conditions
    ")

    #=
    #find the first non-singular CL
    aCL1 = [aCL[1]]
    ECL!(aCL1, P; incL=true)
    =#

end

main()
nothing