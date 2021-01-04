using Stochastics
using Test


@testset "Regression.jl" begin
    x₁ = [56.7, 66.6, 55.7, 63.5, 63.1, 57.6, 57.8, 60.4, 57.4, 59.2, 59.3] 
    y₁ = [7.5, -0.6,8.5, 4.0, 4.2, 7.6, 7.7, 5.8, 6.0, 7.0, 7.0]
    d = Stochastics.regression_linear_simple(x₁, y₁)
    @test round(get(d, "α", nothing), digits=4) == 48.9646
    @test round(get(d, "β", nothing), digits=4) == -0.7210
    @test round(√(get(d, "s²", nothing)), digits=4) == 0.9431
end