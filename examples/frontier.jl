
#source https://github.com/ungil/Markowitz.jl/blob/master/examples/frontier.jl

using EfficientFrontier
#using Pkg; Pkg.add("EfficientFrontier")
#using Pkg; Pkg.add("https://github.com/PharosAbad/EfficientFrontier.jl.git")

#https://stackoverflow.com/a/53030465 List of loaded/imported packages in Julia
#filter((x) -> typeof(eval(x)) <:  Module && x ≠ :Main, names(Main,imported=true))
if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
    #https://stackoverflow.com/questions/71595632 How include a file module in Julia 1.7?   -> First option
    include("./Markowitz.jl")
    using .Markowitz
end

#=
if length(filter((x) -> x == :EfficientFrontier, names(Main, imported=true))) == 0
    include("../src/EfficientFrontier.jl")
    using .EfficientFrontier
end
=#



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

E = [ 0.1 0.7 0.8 2.3 2.2 1.9 5.6 5.6 2.2 1.3 0.7 -0.1 4.1 7.2 ]

m = markowitz(E, V, names=assets)
unit_sum(m) # total weight = 100%
f = frontier(m)
display(f.weights)


N = length(E)

A = ones(1, N)
b = ones(1)
G = Matrix{Float64}(undef, 0, N)
g = Vector{Float64}(undef, 0)
d = zeros(Float64, N)
u = fill(Inf, N)

#aCL = Vector{sCL}(undef, 0)

EfficientFrontier.setup(E, V, A, b, d, u, G, g)
aCL = EfficientFrontier.ECL()
#display(EfficientFrontier.aCL)
display(aCL)

Z = EfficientFrontier.CornerP()

#f.weights[3:end,:]-Z[2:end,:]
