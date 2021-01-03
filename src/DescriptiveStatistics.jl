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
function sample_covariance(x::Vector, y::Vector)
    n = length(x)
    cᵪᵧ = 1 / (n -1) * (sum(x .* y) - 1 / n * sum(x) * sum(y))
    return cᵪᵧ
end


"""
"""
function correlation_coefficient(x::Vector, y::Vector)
    sᵪ = √(sample_variance(x))
    sᵧ = √(sample_variance(y))
    cᵪᵧ = sample_covariance(x, y)
    rᵪᵧ = cᵪᵧ / (sᵪ * sᵧ)
    return rᵪᵧ
end
