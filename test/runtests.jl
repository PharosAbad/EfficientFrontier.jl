using EfficientFrontier
using Test

@testset "EfficientFrontier.jl" begin
    # Write your tests here.
    V = [1/100 1/80 1/100
        1/80 1/16 1/40
        1/100 1/40 1/25]
    E = [109 / 100; 23 / 20; 119 / 100]

    EfficientFrontier.noShortsale(E, V)
    aCL = EfficientFrontier.ECL()
    Z = EfficientFrontier.CornerP()
end
