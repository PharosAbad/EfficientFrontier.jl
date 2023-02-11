#Example: A non-endpoint corner portfolio is a degenerated portfolio (K=0, weights are on the lower or upper bounds)
#R: endpoint corner portfolios are: LVEP (Lowest Variance Efficient Portfolio), also called GMVP (Global Minimum Variance Portfolio) 
# and HVEP (Highest Variance Efficient Portfolio), also called HMFP (Highest Mean Frontier Portfolio) 

using EfficientFrontier, LinearAlgebra

if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    include("./Markowitz.jl")
    using .Markowitz
end

function main()
    E, V = EfficientFrontier.EVdata(:Ungil, false)
    u = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.3 0.3 0.3 0.3]
    u = vec(u)
    P = Problem(E, V, u)
    println("\n--- Status-Segment Method: A non-endpoint degenerated corner portfolio`  ---\n")
    ts = @elapsed aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)
    println("Status-Segment Method:  ", ts, "  seconds")    #0.000936215  seconds
    println("\n---  aCL[10:11]  are as follows")
    display(aCL[10:11])
    println("          which shows that from L=", aCL[10].L0, " down to L=", aCL[11].L1, ", we have K=0, a fully degenerated corner portfolio\n")

    println("\n Weigths for the first 9 corner portfolios")
    display(aEF.Z[1:9, :])
    
    println("\n Markowitz CLA:  starting ...")
    m = markowitz(E, V; upper=u)
    unit_sum(m) # total weight = 100%
    ts = @elapsed f = frontier(m)
    println("Markowitz CLA:  ", ts, "  seconds")     #0.14 seconds
    println("Weigths by CLA: first 12 corner portfolios")
    display(f.weights[1:12, :])
    #norm(aEF.Z[1:2,:]-f.weights[1:2,:], Inf)
    #(aEF.Z[3,:]-f.weights[3,:])'
    println("Mean (expected return) by CLA")
    display(f.ret[1:12]')
    #norm(aEF.Z[9:end,:]-f.weights[12:end,:], Inf)
    println("there are 3 kinks (see the repeated mu: 5.29, 5.1 and 4.42)")


    println("\n

    Remark: endpoint corner portfolios are 
    * LMEP (Lowest Mean Efficient Portfolio) == LVEP (Lowest Variance Efficient Portfolio) == GMVP (Global Minimum Variance Portfolio)
    * HMEP (Highest Mean Efficient Portfolio) == HVEP (Highest Variance Efficient Portfolio) == HMFP (Highest Mean Frontier Portfolio)

    ")
end

main()
nothing