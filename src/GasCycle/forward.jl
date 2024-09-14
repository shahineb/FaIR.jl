include("../Inputs/Inputs.jl")
include("gas_cycle.jl")
using .Inputs


function EtoC(gcm::GasCycleModel, Eₜ::Vector{Real}, Burdenₜ₋₁::Matrix{Real}, αₜ::Vector{Real}, E₀::Vector{Real}, Δt::Real)
    δ = Δt ./ (αₜ .* gcm.τ)
    e⁻ᵟ = exp(-δ)

    Burdenₜ = gcm.a .* (Eₜ .- E₀) .* (1 ./ δ) .* Δt * (1 .- e⁻ᵟ) + Burdenₜ₋₁ * e⁻ᵟ
    Airborneₜ = sum(Burdenₜ, dims=2)
    Cₜ = gcm.C₀ .+ gcm.EtoC .* Airborneₜ
    return Cₜ, Burdenₜ
end


function α(gcm::GasCycleModel, Airborneₜ::Matrix{Real}, Cumulativeₜ::Matrix{Real}, T::Vector{Real}, iirfmax::Real)
    iirf = gcm.r0 + gcm.ru * (Cumulativeₜ - Airborneₜ) + gcm.rT * T + gcm.ra * Airborneₜ
    iirf = max(irrf, iirfmax)
    αₜ = gcm.g₀ * exp(iirf / gcm.g₁) # may be nans here need attention 
    return αₜ
end

E = Emissions("src/defaults/historical-emissions.csv", ["CO2", "CH4", "SO2", "BC"])
gascycle = GasCycleModel("src/defaults/gas_parameters.csv", ["CO2", "CH4", "SO2", "BC"])

n_pool = 4
Eₜ = E.data[:, 1]
Burdenₜ₋₁ = zeros(Real, 4, n_pool)
αₜ = ones(Real, 4)
E₀ = E.data[:, 1]
Δt = 1.

Cₜ, Burdenₜ = EtoC(gascycle, Eₜ, Burdenₜ₋₁, αₜ, E₀, Δt)
