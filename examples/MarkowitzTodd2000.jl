# lower and upper bounds, 2 inequality constraints

using EfficientFrontier

#=
   if length(filter((x) -> x == :Markowitz, names(Main, imported=true))) == 0
      include("./Markowitz.jl")
      using .Markowitz
   end   
   =#

function main()

   # from  CLA-Data.csv https://github.com/mdengler/cla
   E, V = EfficientFrontier.EVdata(:Markowitz, false)
   N = length(E)
   #A = ones(1, N)
   #b = ones(1)
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

   println("--- Markowitz and Todd (2000), chapter 13, pp.337 --- ")

   #=
   m = markowitz(E, V, lower=d, upper=u)
   unit_sum(m)
   add_constraint(m, G[1,:], '<', g[1])
   add_constraint(m, G[2,:], '<', g[2])
   ts = @elapsed f = frontier(m)
   println("Markowitz CLA:  ", ts, "  seconds")    #0.0007 seconds
   display(f.weights)
   =#


   P = Problem(E, V, u, d, G, g)
   ts = @elapsed aCL = EfficientFrontier.ECL(P)
   #EfficientFrontier.ECL(P; init=cbCL!)
   aEF = eFrontier(aCL, P)
   println("Status-Segment Method:  ", ts, "  seconds")   #0.0004 seconds
   display(aEF.Z)

   # the CLA in Markowitz and Todd (2000) assumes either one IN or one OUT.
   # Our code handles two or more IN and/or OUT.  The first CL in this example has: 
   # Event[Event(UP, IN, 2, 2.494839663636366), Event(DN, IN, 10, 2.494839663636377)]. 
   # Because the first Corner Portfolio is on the boundary, there is no IN assets
   #aCL[1].I1[1].L - aCL[1].I1[2].L    #-1.1102230246251565e-14

   #=
   #BigFloat
   Pb = Problem(convert(Vector{BigFloat}, E), V, u, d, G, g)
   #ts = @elapsed aCLb = EfficientFrontier.ECL(Pb; numSettings = Settings{BigFloat}(tolL = BigFloat(2)^-51))
   ts = @elapsed aCLb = EfficientFrontier.ECL(Pb; numSettings = Settings(Pb; tolL = BigFloat(2)^-51))
   aEFb = eFrontier(aCLb, Pb)
   println("BigFloat:  ", ts, "  seconds")   #0.020 seconds
   #display(aEFb.Z)
   #maximum(abs.(aEFb.Z-aEF.Z))   #3.344546861683284e-15

   #=
   aCLb[1].I1
   2-element Vector{Event{BigFloat}}:
    Event{BigFloat}(UP, IN, 2, 2.494839663636366245529450347001271614280350227856738443120461658250300005883108)
    Event{BigFloat}(DN, IN, 10, 2.494839663636366495420530397802951875034152818858742312229266642281434186751967)
   =#
   #aCLb[1].I1[1].L - aCLb[1].I1[2].L  #-2.498910800508016802607538025910020038691088049840311341808673741032753977364765e-16
   #for BigFloat, default tolL = BigFloat(2)^-76 = 1.32348898008484427979425390731194056570529937744140625e-23
   =#

   #BigFloat from the source data: convert to Int (rationals) then to BigFloat to improve precision
   V0 = convert(Matrix{Int64}, round.(V * 10^8))
   E0 = convert(Vector{Int64}, round.(E * 10^3))
   d0 = convert(Vector{Int64}, round.(d * 10))
   u0 = convert(Vector{Int64}, round.(u * 10))
   g0 = convert(Vector{Int64}, round.(g * 10))
   G0 = convert(Matrix{Int64}, round.(G * 10))
   Vb = V0 / BigFloat(10^8)
   Eb = E0 / BigFloat(10^3)
   db = d0 / BigFloat(10)
   ub = u0 / BigFloat(10)
   gb = g0 / BigFloat(10)
   Gb = G0 / BigFloat(10)
   #println("\n xxx  ", maximum(abs.(Eb-E)))
   println("\n E[4]=1.12,  but BigFloat(\"1.12\")-E[4] is: ", BigFloat("1.12") - E[4])
   println("now convert the raw data to BigFloat to improve precision (raw data to BigFloat, not raw data to Float64 then to BigFloat)")


   #S = Status[UP, IN, DN, UP, DN, DN, DN, DN, DN, IN, OE, OE]
   S = aCL[1].S
   Pb = Problem(Eb, Vb, ub, db, Gb, gb)
   #aCLb = Vector{sCL{typeof(Pb).parameters[1]}}(undef, 0)
   aCLb = sCL(Pb)
   computeCL!(aCLb, S, Pb, Settings(Pb))
   #computeCL!(aCLb, S, Pb, Settings{BigFloat}(tolL = BigFloat(2)^-51)) #if the data is first truncated by Float64, adjust the TolL


   println("\n Assets 2 and 10 go IN at the same time, NOT the case that ONE at each time")
   display(aCLb[1].I1)
   println("using BigFloat, the difference on L: ", aCLb[1].I1[1].L - aCLb[1].I1[2].L)
   #aCLb[1].I1[1].L -aCLb[1].I1[2].L        #-1.450876317255866697064907112950467127947488061225295272683982182988323422931288e-75



   println("
   Because the first Corner Portfolio is on the boundary, there is no IN assets (K=0).
      Hence,  K>=2 assets go IN (since z'1=1, no way to change a single one while others fixed on boundary, thus, K=1 is ruled out)

   the CLA in Markowitz and Todd (2000) assumes either one IN or one OUT (if and only if one asset changes state).
   Our code handles two or more IN and/or OUT (allow any number of changes concurrently).
   ")

end

main()
nothing
