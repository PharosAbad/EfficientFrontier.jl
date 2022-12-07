#source https://github.com/ungil/Markowitz.jl/blob/master/examples/frontier.jl
#Example: # lower and upper bounds, 2 inequality constraints, 2 equality constraints

using EfficientFrontier
#using Pkg; Pkg.add("EfficientFrontier")
#using Pkg; Pkg.add("https://github.com/PharosAbad/EfficientFrontier.jl.git")

if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
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


class = vec([:FI :FI :FI :FI :FI :FI :FI :ALT :ALT :ALT :EQ :EQ :EQ :EQ])
m = markowitz(E, V, names=assets, # asset bounds by class: stocks -10/30, bonds 0/20, alt. 0/10
            lower = -0.1 * (class .== :EQ),
            upper = 0.3 * (class .== :EQ) + 0.2 * (class .== :FI) + 0.1 * (class .== :ALT))
unit_sum(m)
add_constraint(m, 1 * (class .== :EQ), '>', 0.3) # net equity exposure between 30% and 60%
add_constraint(m, 1 * (class .== :EQ), '<', 0.6)
add_constraint(m, [1 1 0 0 0 0 0 0 0 0 0 0 0 0], '=', 0.25) # US govt + Investment Grade = 25%

ts = @elapsed f = frontier(m)
println("Markowitz CLA:  ", ts, "  seconds") #0.001 seconds
display(f.weights)



#EfficientFrontier
E = vec(E)
A = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
b = [1.0; 0.25]
G = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0]
G[1, :] = -G[1, :]
g = [-0.3; 0.6]
d = vec([-0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.1 -0.1 -0.1 -0.1])
u = vec([0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.3 0.3 0.3 0.3])

#=  #v0.1.0
EfficientFrontier.setup(E, V, A, b, d, u, G, g)
aCL = EfficientFrontier.ECL()
display(aCL)    #there are 2  singular CL (beta is a zero vector, a line in R^(N+J+M+1) space, but a single point in mean-variance space )
Z = EfficientFrontier.CornerP() # Kinks on the frontier due to singular CL
=#


#v0.2.0
P = Problem(E, V, u, d, G, g, A, b)
ts = @elapsed aCL = EfficientFrontier.ECL(P)
aEF = eFrontier(aCL, P)
#display(aCL)    #there are 2  singular CL (beta is a zero vector, a line in R^(N+J+M+1) space, but a single point in mean-variance space )
println("connecting Critical Line Segments:  ", ts, "  seconds") #0.001 seconds
display(aEF.Z)  # Kinks on the frontier due to singular CL

Pt = Problem(E, V, u, d, G, g, A, b; equilibrate=true)
ts = @elapsed aCLt = EfficientFrontier.ECL(Pt)
aEFt = eFrontier(aCLt, Pt)
println("equilibrate:  ", ts, "  seconds") #0.001 seconds


#BigFloat
Pb = Problem(convert(Vector{BigFloat},E), V, u, d, G, g, A, b)
ts = @elapsed aCLb = EfficientFrontier.ECL(Pb)
aEFb = eFrontier(aCLb, Pb)
println("BigFloat:  ", ts, "  seconds")   #0.037 seconds, slower than Float64
println("improvements  ", round.([maximum(abs.(aEFb.Z-aEF.Z)), maximum(abs.(aEFt.Z-aEF.Z))], sigdigits=3))
#maximum(abs.(z[end-15:end,:]-aEF.Z[end-15:end,:]))     #5.495603971894525e-15

println("
For precision,  `Float64+equilibrate` is enough, it is approximate to BigFloat. However, BigFloat can be tens (30 to 70), even more than 200 times slow.
")
