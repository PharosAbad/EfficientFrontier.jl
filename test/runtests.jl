using EfficientFrontier
using Test

@testset "EfficientFrontier.jl" begin
    # Write your tests here.
    V = [1/100 1/80 1/100
        1/80 1/16 1/40
        1/100 1/40 1/25]
    E = [109 / 100; 23 / 20; 119 / 100]

        
    #v0.3.0
    P = Problem(E, V; equilibrate=false)
    aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)

    @test length(aCL) == 4
    @test isapprox(aEF.Z[end,:]'*E, E[1])
end
