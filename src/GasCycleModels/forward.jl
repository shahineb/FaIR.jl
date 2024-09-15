include("../Inputs/Inputs.jl")
include("reservoir_model.jl")
using .Inputs


function EtoC(gcm::ReservoirModel, Eₜ::AbstractVector{<:Real}, pool_partition::AbstractMatrix{<:Real}, αₜ::AbstractVector{<:Real}, E₀::AbstractVector{<:Real}, Δt::Real)
    δ = Δt ./ (αₜ .* gcm.τ)
    e⁻ᵟ = exp.(-δ)
    pool_partition = gcm.a .* (Eₜ .- E₀) .* (1 ./ δ) .* Δt .* (1 .- e⁻ᵟ) .+ pool_partition .* e⁻ᵟ
    airborneₜ₊₁ = sum(pool_partition, dims=2)
    Cₜ₊₁ = gcm.C₀ .+ gcm.EtoC .* airborneₜ₊₁
    return Cₜ₊₁, pool_partition, airborneₜ₊₁
end


function α(gcm::ReservoirModel, airborneₜ::AbstractVector{<:Real}, cumulativeₜ::AbstractVector{<:Real}, Tₜ::Real, iirfmax::Real)
    iirf = gcm.r0 .+ gcm.ru .* (cumulativeₜ .- airborneₜ) .+ gcm.rT .* Tₜ .+ gcm.ra .* airborneₜ
    iirf = min.(iirf, iirfmax)
    αₜ = gcm.g₀ .* exp.(iirf ./ gcm.g₁)
    return αₜ
end