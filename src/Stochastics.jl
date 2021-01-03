module Stochastics


export BinomialDistribution, binomial_mean, binomial_variance, binomial_probability_mass_funcion


abstract type Distribution 
end


abstract type DiscreteDistribution <: Distribution 
end


"""
BinomialDistribution creates a struct for the binomial distribution. 

See also: [`binomial_mean`](@ref), [`binomial_variance`](@ref), [`binomial_probability_mass_funcion`](@ref)
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


abstract type ContiniousDistribution <: Distribution 
end


end
