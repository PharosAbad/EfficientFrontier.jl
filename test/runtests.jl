using EfficientFrontier
using Test

@testset "EfficientFrontier.jl" begin
    # Write your tests here.
    V = [1/100 1/80 1/100
        1/80 1/16 1/40
        1/100 1/40 1/25]
    E = [109 / 100; 23 / 20; 119 / 100]

    #=  #v0.1.0
    EfficientFrontier.noShortsale(E, V)
    aCL = EfficientFrontier.ECL()
    Z = EfficientFrontier.CornerP()
    @test isapprox(Z[end,:]'*E, E[1])
    =#
    
    #v0.2.0
    P = Problem(E, V; equilibrate=false)
    aCL = EfficientFrontier.ECL(P)
    aEF = EfficientFrontier.eFrontier(aCL, P)

    @test length(aCL) == 3
    @test isapprox(aEF.Z[end,:]'*E, E[1])
end
