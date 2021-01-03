include("./DescriptiveStatistics.jl")


"""
"""
function normal_reference__mu(μ::AbstractFloat, x::Vector, σ²=nothing)
    n = length(x) 
    x̅ = sample_mean(x)
    if σ² === nothing
        s = √(sample_variance(x))
        Rμ = (x̅ - μ) / (s * √(n))
        return Rμ
    if σ² ≠ nothing
        σ = √(σ²)
        Rμ = (x̅ - μ) / (σ * √(n))
        return Rμ
    end 
end

"""
"""
function normal_reference_sigma(σ²::AbstractFloat, x::Vector)
    n = length(x)
    s² = sample_variance(x)
    Rσ² = (n - 1) * s² / σ²
end
