using Stochastics
using Test


@testset "Distributions.jl" begin
    X = Stochastics.NormalDistribution(0, 1)
    @test abs2(Stochastics.normal_cumulative_distribution_function(X, 1.96) - 0.975) < 1e-8
end
