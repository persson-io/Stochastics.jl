"""
Calculates the reference variable for ``μ`` when the elements in vector ``x`` comes from a normal distribution.
See also: [`sample_mean`](@ref), [`sample_variance`](@ref).
"""
function normal_reference_mu_unknown_variance(μ::AbstractFloat, x::Vector)
    n = length(x) 
    x̅ = sample_mean(x)
    s = √(sample_variance(x))
    T = (x̅ - μ) / (s / √(n))
    return T
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
    T = (x̅ - μ) / (σ / √(n))
    return T 
end


"""
Reference variable for ``σ²`` when the elements in vector ``x`` comes from a normal distribution.
See also: [`sample_variance`](@ref).
"""
function normal_reference_sigma(σ²::AbstractFloat, x::Vector)
    n = length(x)
    s² = sample_variance(x)
    T = (n - 1) * s² / σ²
    return T
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
    T = ((x̅ - y̅) - (μᵪ - μᵧ)) / √(σ²ᵪ / nᵪ + σ²ᵧ / nᵧ)
    return T
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
    T = ((x̅ - y̅) - (μᵪ - μᵧ)) / (s * √(1 / nᵪ + 1 / nᵧ))  
    return T
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
    T = ((x̅ - y̅) - (μᵪ - μᵧ)) / √(s²ᵪ / nᵪ + s²ᵧ / nᵧ)
    return T   
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


"""
Calculates the interval for the estimated mean difference of two paired samples.
See also: [`sample_mean`](@ref), [`sample_variance`](@ref).  
"""
function normal_interval_two_sample_paired(x::Vector, y::Vector, α::AbstractFloat)
    n = length(x)
    z = y .- x 
    z̅ = sample_mean(z)
    dist = tDistribution(n-1)  
    t = quantile_finder(dist, α)
    s = √(sample_variance(z))
    I = [z̅ - t * s / √(n), z̅ + t * s / √(n)]  
    return I
end


function normal_interval_mu_unknown_variance(μ::AbstractFloat, x::Vector, α::AbstractFloat)
    
end


function normal_interval_mu_known_variance(μ::AbstractFloat, σ²::AbstractFloat, x::Vector, α::AbstractFloat)
    
end


function normal_interval_sigma(x::Vector, α::AbstractFloat)
    
end


function normal_interval_two_sample_mu_known_variance(μᵪ::AbstractFloat, μᵧ::AbstractFloat, σ²ᵪ::AbstractFloat, σ²ᵧ::AbstractFloat, x::Vector, y::Vector, α::AbstractFloat)
    
end


function normal_interval_two_sample_mu_unknown_equal_variance(μᵪ::AbstractFloat, μᵧ::AbstractFloat, x::Vector, y::Vector, α::AbstractFloat)
    
end


function exp_reference_mu(μ::AbstractFloat, x::Vector)
    n = length(n)
    x̅ = sample_mean(x)
    T = (x̅ - μ) / (μ / √(n)) 
    return T
end


function exp_interval_mu(μ::AbstractFloat, x::Vector, α::AbstractFloat)
    
end


function binomial_reference_p()
    
end


function binomial_interval_p()
    
end


function hypergeometric_reference_p()
    
end


function hypergeometric_interval_p()
    
end


function binomial_interval_comparison_of_proportions()
    
end


function poisson_reference_mu()
    
end


function poisson_interval_mu()
    
end


function test()
    
end
