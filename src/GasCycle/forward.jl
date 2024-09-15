include("../Inputs/Inputs.jl")
include("gas_cycle.jl")
using .Inputs


function EtoC(gcm::GasCycleModel, Eₜ::AbstractVector{<:Real}, pool_partition::AbstractMatrix{<:Real}, αₜ::AbstractVector{<:Real}, E₀::AbstractVector{<:Real}, Δt::Real)
    δ = Δt ./ (αₜ .* gcm.τ)
    e⁻ᵟ = exp.(-δ)
    pool_partition = gcm.a .* (Eₜ .- E₀) .* (1 ./ δ) .* Δt * (1 .- e⁻ᵟ) .+ pool_partition .* e⁻ᵟ
    airborneₜ₊₁ = sum(pool_partition, dims=2)
    Cₜ₊₁ = gcm.C₀ .+ gcm.EtoC .* airborneₜ₊₁
    return Cₜ₊₁, pool_partition, airborneₜ₊₁
end


function α(gcm::GasCycleModel, airborneₜ::AbstractVector{<:Real}, cumulativeₜ::AbstractVector{<:Real}, Tₜ::Real, iirfmax::Real)
    iirf = gcm.r0 .+ gcm.ru .* (cumulativeₜ .- airborneₜ) .+ gcm.rT .* Tₜ .+ gcm.ra .* airborneₜ
    iirf = min.(iirf, iirfmax)
    αₜ = gcm.g₀ .* exp.(iirf ./ gcm.g₁) # may be nans here need attention 
    return αₜ
end

# E = Emissions("src/defaults/ssp245-emissions.csv", ["CO2", "CH4", "SO2", "BC"])
# gascycle = GasCycleModel("src/defaults/gas_parameters.csv", ["CO2", "CH4", "SO2", "BC"])


# foo = GasCycleModel(["CO2", "CH4", "SO2", "BC"],
#                     ones(Real, 4, 4),
#                     ones(Real, 4, 4),
#                     [1., 2., 3., 4.],
#                     [1., 2., 3., 4.],
#                     [1., 2., 3., 4.],
#                     [1., 2., 3., 4.],
#                     ones(Real, 4),
#                     ones(Real, 4),
#                     ones(Real, 4, 4),
#                     [1., 2., 3., 4.],
#                     [1., 2., 3., 4.])


# x = convert(Vector{Real}, [1., 2., 3., 4.])
# y = [20.08553692,  109.19630007,  445.23947731, 1613.71517397]
# y ≈ α(foo, x, x, 1., 100.)


# n_pool = 4
# t = 1
# Eₜ = E.values[:, t]
# pool_partition = zeros(Real, 4, n_pool)
# airborneₜ = zeros(Real, 4)
# cumulativeₜ = E.cumulative[:, t]
# Tₜ = 0.
# iirfmax = 100.
# αₜ = α(gascycle, airborneₜ, cumulativeₜ, Tₜ, iirfmax)



# E₀ = E.values[:, 1]
# Δt = 1.
# Cₜ, Burdenₜ = EtoC(gascycle, Eₜ, pool_partition, αₜ, E₀, Δt)
