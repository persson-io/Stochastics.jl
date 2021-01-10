using SafeTestsets


@safetestset "Inference Tests" begin include("inferencetests.jl") end
@safetestset "Descriptive Statistics Tests" begin include("descriptivestatisticstests.jl") end
@safetestset "Regression Tests" begin include("regressiontests.jl") end
@safetestset "Distributions Tests" begin include("distributionstests.jl") end
