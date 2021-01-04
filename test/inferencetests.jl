using Stochastics
using Test


@testset "Inference.jl" begin
    x = [1, 2]
    y = [3, 4]
    Stochastics.two_sample_degrees_of_freedom_unknown_variance(x, y) == 2
end
