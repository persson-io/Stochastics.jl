include("DescriptiveStatistics.jl")


"""
"""
function normal_reference__mu(μ::AbstractFloat, x::Vector, σ²=nothing)
    n = length(x) 
    x̅ = sample_mean(x)
    if σ² === nothing
        s = √(sample_variance(x))
        Rμ = (x̅ - μ) / (s * √(n))
        return Rμ
    elseif σ² ≠ nothing
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
    return Rσ²
end


"""
"""
function normal_reference_two_sample_mu()
    return nothing    
end


"""
"""
function two_sample_degrees_of_freedom_unknown_variance(x::Vector, y::Vector)
    nᵪ = length(x)
    nᵧ = length(y)
    s²ᵪ = sample_variance(x)
    s²ᵧ = sample_variance(y)
    first = 1 / (nᵪ - 1)
    second = (nᵧ * s²ᵪ)^2 / (nᵧ * s²ᵪ + nᵪ * s²ᵧ)^2
    third = 1 / (nᵧ - 1)
    fourth = (nᵪ * s²ᵧ)^2 / (nᵧ * s²ᵪ + nᵪ * s²ᵧ)^2
    f = 1 / (first * second + third * fourth)
    return f
end
