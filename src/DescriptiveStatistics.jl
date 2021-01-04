"""
"""
function sample_mean(x::Vector)
    n = length(x)
    x̅ = 1 / n * sum(x)
    return x̅
end


"""
"""
function sample_variance(x::Vector)
    n = length(x)
    s² = 1 / (n - 1) * (sum(x.^2) - 1 / n * sum(x)^2)
    return s² 
end


"""
"""
function two_sample_covariance(x::Vector, y::Vector)
    n = length(x)
    cᵪᵧ = 1 / (n -1) * (sum(x .* y) - 1 / n * sum(x) * sum(y))
    return cᵪᵧ
end


"""
"""
function two_sample_correlation_coefficient(x::Vector, y::Vector)
    sᵪ = √(sample_variance(x))
    sᵧ = √(sample_variance(y))
    cᵪᵧ = two_sample_covariance(x, y)
    rᵪᵧ = cᵪᵧ / (sᵪ * sᵧ)
    return rᵪᵧ
end


"""
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
"""
function sum_of_squares(x::Vector, y=x::Vector)
    n = length(x)
    Sᵪᵧ = (sum(x .* y) - 1 / n * sum(x) * sum(y))
    return Sᵪᵧ 
end
