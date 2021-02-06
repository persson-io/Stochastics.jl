using Stochastics
using Test


@testset "Distributions.jl" begin
    X = Stochastics.NormalDistribution(0, 1)
    @test abs2(Stochastics.normal_cumulative_distribution_function(X, 1.96) - 0.975) < 1e-8
    @test abs2(Stochastics.gamma_function(1) - 1) < 1e-8
    @test abs2(Stochastics.gamma_function(1/2) - √(π)) < 1e-8
    @test_throws ErrorException("α must be in ℝ > 0") Stochastics.gamma_function(0)
end
