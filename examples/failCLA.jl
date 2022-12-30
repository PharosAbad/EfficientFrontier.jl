#Example: CLA (Critical Line Algorithm) fail, for it only handle either one IN or one OUT


using EfficientFrontier
#using Pkg; Pkg.add("EfficientFrontier")
#using Pkg; Pkg.add("https://github.com/PharosAbad/EfficientFrontier.jl.git")

if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    include("./Markowitz.jl")
    using .Markowitz
end

println("\n--- connecting Critical Line Segments vs Markowitz's CLA  ---\n")

V = [ 133.0  39.0  36.0  -17.0  -28.0   31.0   -5.0   -6.0  -34.0  -10.0  -46.0  -68.0  -52.0  -63.0
       39.0  16.0  17.0   10.0    0.0   22.0   13.0    7.0    1.0   19.0    3.0    0.0    5.0    9.0
       36.0  17.0  38.0   19.0    3.0   34.0   45.0   29.0   37.0   48.0   20.0   21.0   42.0   50.0
      -17.0  10.0  19.0   98.0   69.0   65.0   73.0   43.0   68.0  151.0   99.0  130.0  121.0  165.0
      -28.0   0.0   3.0   69.0   84.0   34.0   42.0   25.0   55.0   99.0   64.0   93.0   87.0  119.0
       31.0  22.0  34.0   65.0   34.0   93.0   83.0   48.0   57.0  115.0   76.0   85.0   99.0  145.0
       -5.0  13.0  45.0   73.0   42.0   83.0  154.0   87.0  111.0  168.0  113.0  143.0  165.0  225.0
       -6.0   7.0  29.0   43.0   25.0   48.0   87.0   57.0   75.0   95.0   71.0   89.0  104.0  142.0
      -34.0   1.0  37.0   68.0   55.0   57.0  111.0   75.0  286.0  117.0  101.0  125.0  164.0  239.0
      -10.0  19.0  48.0  151.0   99.0  115.0  168.0   95.0  117.0  491.0  231.0  327.0  259.0  322.0
      -46.0   3.0  20.0   99.0   64.0   76.0  113.0   71.0  101.0  231.0  214.0  256.0  217.0  269.0
      -68.0   0.0  21.0  130.0   93.0   85.0  143.0   89.0  125.0  327.0  256.0  380.0  265.0  342.0
      -52.0   5.0  42.0  121.0   87.0   99.0  165.0  104.0  164.0  259.0  217.0  265.0  297.0  359.0
      -63.0   9.0  50.0  165.0  119.0  145.0  225.0  142.0  239.0  322.0  269.0  342.0  359.0  556.0 ]

#E = [ 0.1 0.7 0.8 2.3 2.2 1.9 5.6 5.6 2.2 1.3 0.7 -0.1 4.1 7.2 ]
E = [7.2 0.7 0.8 2.3 2.2 1.9 5.6 5.6 7.2 1.3 0.7 -0.1 4.1 7.2] # 3 assets share the highest expected return

m = markowitz(E, V)
unit_sum(m) # total weight = 100%
ts = @elapsed f = frontier(m)   #Warning: tweaking mu[8] mu[13] to ensure the solution is unique  => perturbed method
println("Markowitz CLA:  ", ts, "  seconds")    #0.00059 seconds
display(f.weights)


P = Problem(E, V; equilibrate=false)
ts = @elapsed aCL = EfficientFrontier.ECL(P)
aEF = eFrontier(aCL, P)
println("connecting Critical Line Segments:  ", ts, "  seconds")   #0.00077 seconds
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