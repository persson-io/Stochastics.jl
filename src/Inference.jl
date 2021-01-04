include("DescriptiveStatistics.jl")


"""
Calculates the reference variable for ``μ`` when the elements in vector ``x`` comes from a normal distribution.
See also: [`sample_mean`](@ref), [`sample_variance`](@ref).
"""
function normal_reference_mu_unknown_variance(μ::AbstractFloat, x::Vector)
    n = length(x) 
    x̅ = sample_mean(x)
    s = √(sample_variance(x))
    Rμ = (x̅ - μ) / (s * √(n))
    return Rμ
end


"""
Calculates the reference variable for ``μ`` when the elements in vector ``x`` comes from a normal distribution.
This time the variance ``σ²`` has to be known for the distribution.
See also: [`sample_mean`](@ref).
"""
function normal_reference_mu_known_variance(μ::AbstractFloat, σ²::AbstractFloat, x::Vector)
    n = length(x) 
    x̅ = sample_mean(x)
    σ = √(σ²)
    Rμ = (x̅ - μ) / (σ * √(n))
    return Rμ 
end


"""
Reference variable for ``σ²`` when the elements in vector ``x`` comes from a normal distribution.
See also: [`sample_variance`](@ref).
"""
function normal_reference_sigma(σ²::AbstractFloat, x::Vector)
    n = length(x)
    s² = sample_variance(x)
    Rσ² = (n - 1) * s² / σ²
    return Rσ²
end


"""
Calculates the reference variable for tests of the difference between two means.
This variable is used when the variance is known and different for the two samples.
See also: [`sample_mean`](@ref).
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
Calculates the reference variable for tests of the difference between two means.
This variable is used when the variance is unknown but equal for the two samples.
See also: [`sample_mean`](@ref), [`two_sample_variance`](@ref).
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
Calculates the reference variable for tests of the difference between two means.
This variable is used when the variance is unknown but not necessarily equal for the two samples.
See also: [`sample_mean`](@ref), [`sample_variance`](@ref), [`two_sample_degrees_of_freedom_unknown_variance`](@ref).
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
Calculates the degrees of freedom of the t-distribution that the reference variable above is from.
See also: [`sample_variance`](@ref), [`normal_reference_two_sample_mu_unknown_variance`](@ref).
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
