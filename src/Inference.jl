include("DescriptiveStatistics.jl")


"""
"""
function normal_reference_mu_unknown_variance(μ::AbstractFloat, x::Vector)
    n = length(x) 
    x̅ = sample_mean(x)
    s = √(sample_variance(x))
    Rμ = (x̅ - μ) / (s * √(n))
    return Rμ
end


"""
"""
function normal_reference_mu_known_variance(μ::AbstractFloat, σ²::AbstractFloat, x::Vector)
    n = length(x) 
    x̅ = sample_mean(x)
    σ = √(σ²)
    Rμ = (x̅ - μ) / (σ * √(n))
    return Rμ 
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
function normal_reference_two_sample_mu_known_variance(μᵪ::AbstractFloat, μᵧ::AbstractFloat, σ²ᵪ::AbstractFloat, σ²ᵧ::AbstractFloat, x::Vector, y::Vector)
    nᵪ = length(x)
    nᵧ = length(y)
    x̅ = sample_mean(x) 
    y̅ = sample_mean(x)
    Rμᵪμᵧ = ((x̅ - y̅) - (μᵪ - μᵧ)) / √(σ²ᵪ / nᵪ + σ²ᵧ / nᵧ)
    return Rμᵪμᵧ
end


"""
"""
function normal_reference_two_sample_mu_unknown_equal_variance(μᵪ::AbstractFloat, μᵧ::AbstractFloat, x::Vector, y::Vector)
    nᵪ = length(x)
    nᵧ = length(y)
    x̅ = sample_mean(x) 
    y̅ = sample_mean(x)
    s = √(two_sample_variance(x, y))
    Rμᵪμᵧ = ((x̅ - y̅) - (μᵪ - μᵧ)) / (s * √(1 / nᵪ + 1 / nᵧ))  
    return Rμᵪμᵧ
end


"""
"""
function normal_reference_two_sample_mu_unknown_variance(μᵪ::AbstractFloat, μᵧ::AbstractFloat, x::Vector, y::Vector)
    nᵪ = length(x)
    nᵧ = length(y)
    x̅ = sample_mean(x) 
    y̅ = sample_mean(x)
    s²ᵪ = sample_variance(x)
    s²ᵧ = sample_variance(y)
    Rμᵪμᵧ = ((x̅ - y̅) - (μᵪ - μᵧ)) / √(s²ᵪ / nᵪ + s²ᵧ / nᵧ)
    return Rμᵪμᵧ   
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
