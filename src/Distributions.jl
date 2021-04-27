using QuadGK
using Roots


abstract type Distribution end


abstract type DiscreteDistribution <: Distribution end


"""
BinomialDistribution creates a struct for the binomial distribution. 
See also: [`binomial_mean`](@ref), [`binomial_variance`](@ref), [`binomial_probability_mass_funcion`](@ref).
"""
struct BinomialDistribution <: DiscreteDistribution
    n::Integer
    p::AbstractFloat
    q::AbstractFloat
    function BinomialDistribution(n::Integer, p::AbstractFloat)
        q = 1 - p
        new(n, p, q)
    end
end


"""
Calculates the binomial mean for a binomial distribution.
"""
function binomial_mean(distribution::BinomialDistribution)
    mean = distribution.n * distribution.p
    return mean
end


"""
Calculates the variance for a binomial distribution.
"""
function binomial_variance(distribution::BinomialDistribution)
    variance = distribution.n * distribution.p * distribution.q
    return variance
end


"""
Calculates the probability that ``P(Y = k)`` for a binomial distribution ``Y``. 
"""
function binomial_probability_mass_funcion(distribution::BinomialDistribution, k::Integer)
    pᵧ = binomial(distribution.n, k) * distribution.p^k * distribution.q^(distribution.n - k)
    return pᵧ 
end


abstract type ContiniousDistribution <: Distribution end


struct NormalDistribution <: ContiniousDistribution
    μ::AbstractFloat
    σ²::AbstractFloat
end


function normal_probability_density_function(distribution::NormalDistribution, x::AbstractFloat)
    μ = distribution.μ
    σ = √(distribution.σ²)
    f = 1 / (σ * √(2 * π)) * ℯ^(-1 / 2 * ((x - μ) / σ)^2)
    return f
end


function normal_cumulative_distribution_function(distribution::NormalDistribution, x::AbstractFloat)
    X = distribution
    F, error = quadgk(t -> normal_probability_density_function(X, t), -Inf, x)
    return F
end


abstract type GammaDistributions <: ContiniousDistribution end


struct GammaDistribution <: GammaDistributions
    α::AbstractFloat
    β::AbstractFloat
end


function gamma_probability_density_function(distribution::GammaDistribution, x::AbstractFloat)
    α = distribution.α
    β = distribution.β
    f = β^α / gamma_function(α) * x^(α - 1) * ℯ^(-β * x)
end


function gamma_cumulative_distribution_function(distribution::GammaDistribution, x::AbstractFloat)
    X = distribution
    F, error = quadgk(t -> gamma_probability_density_function(X, t), -Inf, x)
    return F
end


function gamma_function(α::AbstractFloat)
    α > 0 || error("α must be in ℝ > 0")
    Γ, error = quadgk(x -> x^(α - 1) * exp(-x), 0, Inf)
    return Γ
end


function gamma_function(α::Integer)
    α > 0 || error("α must be in ℝ > 0")
    Γ = factorial(α - 1)
    return Γ
end


struct ExponentialDistribution <: GammaDistributions
    α::Integer
    β::AbstractFloat
    function ExponentialDistribution(β::AbstractFloat)
        α = 1
        return new(α, β)
    end
end


function exponential_probability_density_function(distribution::ExponentialDistribution, x::AbstractFloat)
    f = gamma_probability_density_function(distribution, x)
end


function exponential_cumulative_distribution_function(distribution::ExponentialDistribution, x::AbstractFloat)
    F = gamma_cumulative_distribution_function(distribution, x)
end


struct ChiSquaredDistribution <: GammaDistributions
    df::Integer
    α::AbstractFloat
    β::AbstractFloat
    function ChiSquaredDistribution(df::Integer)
        α = df / 2
        β = 1 / 2
        return new(df, α, β)
    end
end


function chi_squared_probability_density_function(distribution::ChiSquaredDistribution, x::AbstractFloat)
    df = distribution.df
    C = 1 / (gamma_function(df / 2) * 2^(df / 2))
    f = C * x^(df / 2 - 1) * ℯ^(-x / 2)
end


function chi_squared_cumulative_distribution_function(distribution::ChiSquaredDistribution, x::AbstractFloat)
    X = distribution
    F, error = quadgk(t -> chi_squared_probability_density_function(X, t), 0, x)
    return F
end


abstract type tDistributions <: ContiniousDistribution end


struct tDistribution <: tDistributions
    df::Integer
end


function t_probability_density_function(distribution::tDistributions, x::AbstractFloat)
    df = distribution.df
    C = 1 / √(df * π) * gamma_function((df + 1) / 2) / gamma_function(df / 2)
    f = C * (1 + x^2 / df)^((-df - 1) / 2)
end


function t_cumulative_distribution_function(distribution::tDistributions, x::AbstractFloat)
    X = distribution
    F, error = quadgk(t -> t_probability_density_function(X, t), -Inf, x)
    return F
end


struct CauchyDistribution <: tDistributions 
    df::Integer
    function CauchyDistribution()
        df = 1
        new(df)
    end
end


function quantile_finder(distribution::ChiSquaredDistribution, α::AbstractFloat)
    f(x) = chi_squared_cumulative_distribution_function(distribution, x) - 1 + α
    quantile = find_zero(f, (0, 200))
    return quantile
end


function quantile_finder(distribution::NormalDistribution, α::AbstractFloat)
    f(x) = normal_cumulative_distribution_function(distribution, x) - 1 + α
    quantile = find_zero(f, (-5, 5))
    return quantile
end


function quantile_finder(distribution::tDistributions, α::AbstractFloat)
    f(x) = t_cumulative_distribution_function(distribution, x) - 1 + α
    quantile = find_zero(f, (1, 1000))
    return quantile
end
