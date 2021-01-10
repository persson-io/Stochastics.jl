"""
Estimates all parameters for a simple linear regression. 
See also: [`sample_mean`](@ref), [`sum_of_squares`](@ref). 
"""
function regression_linear_simple(x::Vector, y::Vector)
    n = length(x)
    x̅ = sample_mean(x)
    y̅ = sample_mean(y)
    SSE = sum_of_squares(y) - sum_of_squares(x, y)^2 / sum_of_squares(x) 
    SSR = sum_of_squares(x, y)^2 / sum_of_squares(x) 
    SST = SSR + SSR
    R² = 1 - SSE / SST
    s² = SSE / (n - 2)
    β = sum_of_squares(x, y) / sum_of_squares(x) 
    α = y̅ - β * x̅
    return Dict([("α", α), ("β", β), ("s²", s²), ("R²", R²), ("SST", SST), ("SSR", SSR), ("SSE", SSE)]) 
end
