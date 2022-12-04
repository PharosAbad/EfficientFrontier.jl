# lower and upper bounds, 2 inequality constraints

using EfficientFrontier
#using Pkg; Pkg.add("EfficientFrontier")
#using Pkg; Pkg.add("https://github.com/PharosAbad/EfficientFrontier.jl.git")

# from  CLA-Data.csv https://github.com/mdengler/cla
V = [0.40755159 0.03175842 0.05183923 0.05663904 0.0330226 0.00827775 0.02165938 0.01332419 0.0343476 0.02249903
   0.03175842 0.9063047 0.03136385 0.02687256 0.01917172 0.00934384 0.02495043 0.00761036 0.02874874 0.01336866
   0.05183923 0.03136385 0.19490901 0.04408485 0.03006772 0.01322738 0.03525971 0.0115493 0.0427563 0.02057303
   0.05663904 0.02687256 0.04408485 0.19528471 0.02777345 0.00526665 0.01375808 0.00780878 0.02914176 0.01640377
   0.0330226 0.01917172 0.03006772 0.02777345 0.34059105 0.00777055 0.02067844 0.00736409 0.02542657 0.01284075
   0.00827775 0.00934384 0.01322738 0.00526665 0.00777055 0.15983874 0.02105575 0.00518686 0.01723737 0.00723779
   0.02165938 0.02495043 0.03525971 0.01375808 0.02067844 0.02105575 0.68056711 0.01377882 0.04627027 0.01926088
   0.01332419 0.00761036 0.0115493 0.00780878 0.00736409 0.00518686 0.01377882 0.95526918 0.0106553 0.00760955
   0.0343476 0.02874874 0.0427563 0.02914176 0.02542657 0.01723737 0.04627027 0.0106553 0.31681584 0.01854318
   0.02249903 0.01336866 0.02057303 0.01640377 0.01284075 0.00723779 0.01926088 0.00760955 0.01854318 0.11079287]
E = [1.175, 1.19, 0.396, 1.12, 0.346, 0.679, 0.089, 0.73, 0.481, 1.08]

#EfficientFrontier
E0 = vec(E)
V0 = copy(V)
N = length(E)

A = ones(1, N)
b = ones(1)
d = zeros(N)
d[1] = 0.1
d[5] = 0.1
u = 0.3 * ones(N)
g = [-0.2; 0.5]
G = zeros(length(g), N)
G[1, 1:3] .= -1.0
G[1, 4] = -0.5
G[2, 4] = 0.5
G[2, 5:7] .= 1.0

#=  #v0.1.0
EfficientFrontier.setup(E0, V0, A, b, d, u, G, g)
aCL = EfficientFrontier.ECL()
display(aCL)
Z = EfficientFrontier.CornerP() #see Markowitz and Todd (2000), chapter 13, pp.337
=#

#v0.2.0
#P = Problem(E0, V0, u, d, G, g, A, b)
P = Problem(E0, V0, u, d, G, g)
ts = @elapsed aCL = EfficientFrontier.ECL(P)
aEF = EfficientFrontier.eFrontier(aCL, P)
display(ts)   #0.0007 seconds
display(aEF.Z)

# the CLA in Markowitz and Todd (2000) assumes either one IN or one OUT.
# Our code handles two or more IN and/or OUT.  The first CL in this example has: 
# Event[Event(UP, IN, 2, 2.494839663636366), Event(DN, IN, 10, 2.494839663636377)]. 
# Because the first Corner Portfolio is on the boundary, there is no IN assets