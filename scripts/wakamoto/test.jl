# Import the package.
using AdvancedMH
using Distributions
using MCMCChains

using LinearAlgebra

# Generate a set of data from the posterior we want to estimate.
data = rand(Normal(0, 1), 30)

# Define the components of a basic model.
insupport(θ) = θ[2] >= 0
dist(θ) = Normal(θ[1], θ[2])
density(θ) = insupport(θ) ? sum(logpdf.(dist(θ), data)) : -Inf

# Construct a DensityModel.
model = DensityModel(density)

# Set up our sampler with a joint multivariate Normal proposal.
spl = RWMH(MvNormal(zeros(2), I))

# Sample from the posterior.
chain = sample(model, spl, 100000; param_names=["μ", "σ"], chain_type=Chains)
