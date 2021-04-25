using QuadGK


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


struct GammaDistribution <: ContiniousDistribution
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


struct ExponentialDistribution <: ContiniousDistribution
    β::AbstractFloat
    α::Integer
    function ExponentialDistribution(β::AbstractFloat)
        α = 1
        return new(β, α)
    end
end


function exponential_probability_density_function(distribution::ExponentialDistribution, x::AbstractFloat)
    f = gamma_probability_density_function(distribution, x)
end


function exponential_cumulative_distribution_function(distribution::ExponentialDistribution, x::AbstractFloat)
    F = gamma_cumulative_distribution_function(distribution, x)
end


struct ChiSquaredDistribution <: ContiniousDistribution
    f::AbstractFloat
end


function chi_squared_probability_density_function(distribution::ChiSquaredDistribution, t::AbstractFloat)
    f = distribution.f
    C = 2^(-f / 2) / gamma_function(f / 2)
end
