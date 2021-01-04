using Stochastics
using Test


@testset "Inference.jl" begin
    x = [1, 2]
    y = [3, 4]
    x₁ = [5.2, 5.8, 4.3, 5.2, 4.6, 4.7, 5.8, 5.5] 
    y₁ = [5.6, 6.3, 4.9, 5.8, 5.5, 5.7, 6.1, 5.4]
    @test Stochastics.two_sample_degrees_of_freedom_unknown_variance(x, y) == 2
    @test round.(Stochastics.normal_interval_two_sample_paired(x₁, y₁), digits=2) == [0.24, 0.81]
end
