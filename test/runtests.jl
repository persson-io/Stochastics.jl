using SafeTestsets


@safetestset "Inference Tests" begin include("inferencetests.jl") end
@safetestset "Descriptive Statistics Tests" begin include("descriptivestatisticstests.jl") end
