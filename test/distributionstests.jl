using Stochastics
using Test


@testset "Distributions.jl" begin
    X = Stochastics.NormalDistribution(0, 1)
    @test abs(Stochastics.normal_cumulative_distribution_function(X, 1.96) - 0.975) < 1e-4
    @test abs(Stochastics.gamma_function(1) - 1) < 1e-4
    @test abs(Stochastics.gamma_function(1/2) - √(π)) < 1e-4
    @test_throws ErrorException("α must be in ℝ > 0") Stochastics.gamma_function(0)
    X_chi = Stochastics.ChiSquaredDistribution(5)
    @test abs(Stochastics.chi_squared_cumulative_distribution_function(X_chi, 11.070) - 0.95) < 1e-4
    @test abs(Stochastics.quantile_finder(X_chi, 0.9995) - 0.1581) < 1e-4
    @test abs(Stochastics.quantile_finder(X, 0.025) - 1.96) < 1e-4
    X_t = Stochastics.tDistribution(1)
    @test abs(Stochastics.quantile_finder(X_t, 0.1) - 3.0777) < 1e-4
end
