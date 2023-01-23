
#source https://github.com/ungil/Markowitz.jl/blob/master/examples/frontier.jl
#Example: no short-sale

using EfficientFrontier

#https://stackoverflow.com/a/53030465 List of loaded/imported packages in Julia
#filter((x) -> typeof(eval(x)) <:  Module && x ≠ :Main, names(Main,imported=true))
if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    #https://stackoverflow.com/questions/71595632 How include a file module in Julia 1.7?   -> First option
    include("./Markowitz.jl")
    using .Markowitz
end

assets = [ "Bonds - US Government"
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
           "Stocks - Emerging" ]

E, V = EfficientFrontier.EVdata(:Ungil, false)

println("\n--- Status-Segment Method vs Markowitz's CLA  ---\n")

m = markowitz(E, V, names=assets)
unit_sum(m) # total weight = 100%
ts = @elapsed f = frontier(m)
println("Markowitz CLA:  ", ts, "  seconds")    #0.00049  seconds
display(f.weights)



P = Problem(E, V; equilibrate=false)
ts = @elapsed aCL = EfficientFrontier.ECL(P)
aEF = eFrontier(aCL, P)
println("Status-Segment Method:  ", ts, "  seconds")   #0.0006 seconds
display(aEF.Z)
println("NOTE: a kink in second row")

#f.weights[3:end,:]-aEF.Z[2:end,:]
