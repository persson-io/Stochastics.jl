"""
Calculates the average of the elements in a vector ``x``.
"""
function sample_mean(x::Vector)
    n = length(x)
    x̅ = 1 / n * sum(x)
    return x̅
end


"""
Calculates the sample variance of the elements in a vector ``x``.
"""
function sample_variance(x::Vector)
    n = length(x)
    s² = 1 / (n - 1) * (sum(x.^2) - 1 / n * sum(x)^2)
    return s² 
end


"""
Calculates the covariance of vector ``x``.
"""
function sample_covariance(x::Vector)
    n = length(x)
    cᵪᵪ = sum_of_squares(x) / (n - 1)
    return cᵪᵪ 
    
end


"""
Calculates the covariance of two vectors ``x`` and ``y``. 
"""
function two_sample_covariance(x::Vector, y::Vector)
    n = length(x)
    cᵪᵧ = 1 / (n -1) * (sum(x .* y) - 1 / n * sum(x) * sum(y))
    return cᵪᵧ
end


"""
Calculates the correlation coefficient of two vectors ``x`` and ``y``.
See also: [`sample_variance`](@ref), [`two_sample_covariance`](@ref). 
"""
function two_sample_correlation_coefficient(x::Vector, y::Vector)
    sᵪ = √(sample_variance(x))
    sᵧ = √(sample_variance(y))
    cᵪᵧ = two_sample_covariance(x, y)
    rᵪᵧ = cᵪᵧ / (sᵪ * sᵧ)
    return rᵪᵧ
end


"""
Calculates the sample variance of two vectors ``x`` and ``y``.
Used mainly for calculating reference variables for significance tests when the variance is unknown.
See also: [`sample_variance`](@ref). 
"""
function two_sample_variance(x::Vector, y::Vector)
    nᵪ = length(x)
    nᵧ = length(y)
    s²ᵪ = sample_variance(x)
    s²ᵧ = sample_variance(y)
    s²ᵪᵧ = ((nᵪ - 1) * s²ᵪ + (nᵧ - 1) * s²ᵧ) / (nᵪ + nᵧ - 2)
    return s²ᵪᵧ    
end


"""
This function calculates the sum of squares. 
This is used mainly for estimating the parameters in simple linear regression.
See also: [`not_implemented`](@ref).
"""
function sum_of_squares(x::Vector, y=x::Vector)
    n = length(x)
    Sᵪᵧ = (sum(x .* y) - 1 / n * sum(x) * sum(y))
    return Sᵪᵧ 
end
