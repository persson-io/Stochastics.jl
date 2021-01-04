using Stochastics
using Test


@testset "DescriptiveStatistics.jl" begin
    x = [1, 2]
    y = [3, 4]
    @test Stochastics.sample_mean(x) == 1.5
    @test Stochastics.sample_variance(x) == 0.5
    @test Stochastics.sum_of_squares(x) == 0.5
    @test Stochastics.sum_of_squares(x, y) == 0.5
    @test Stochastics.two_sample_covariance(x, y) == 0.5
    @test abs2(Stochastics.two_sample_correlation_coefficient(x, y) - 1) < 10e-8
    @test Stochastics.two_sample_variance(x, y) == 0.5
end
