
struct ReservoirModel
    species::Vector{String}
    a::AbstractMatrix{<:Real}
    τ::AbstractMatrix{<:Real}
    r0::AbstractVector{<:Real}
    ru::AbstractVector{<:Real}
    rT::AbstractVector{<:Real}
    ra::AbstractVector{<:Real}
    C₀::AbstractVector{<:Real}
    EtoC::AbstractVector{<:Real}
    g₀::AbstractVector{<:Real}
    g₁::AbstractVector{<:Real}
end


function ReservoirModel(species::Vector{String}, a::AbstractMatrix{<:Real}, τ::AbstractMatrix{<:Real}, r0::AbstractVector{<:Real}, ru::AbstractVector{<:Real}, rT::AbstractVector{<:Real}, ra::AbstractVector{<:Real}, C₀::AbstractVector{<:Real}, EtoC::AbstractVector{<:Real})
    g₀, g₁ = compute_g(a, τ)
    return ReservoirModel(species, a, τ, r0, ru, rT, ra, C₀, EtoC, g₀, g₁)
end


function ReservoirModel(path::String, species::Vector{String})
    params = load_gascycle_params(path, species)
    ReservoirModel(species, params.a, params.τ, params.r0, params.ru, params.rT, params.ra, params.C₀, params.EtoC)
end


function compute_g(a, τ)
    δ = 100 ./ τ
    e⁻ᵟ = exp.(-δ)
    g₁ = sum(a .* τ .* (1 .- (1 .+ δ) .* e⁻ᵟ), dims=2)
    g₀ = exp.(-sum(a .* τ .* (1 .- e⁻ᵟ), dims=2) ./ g₁)
    return vec(g₀), vec(g₁)
end

function EtoC(gcm::ReservoirModel, Eₜ, pool_partition, αₜ, E₀, Δt)
    δ = Δt ./ (αₜ .* gcm.τ)
    e⁻ᵟ = exp.(-δ)
    pool_partition = gcm.a .* (Eₜ .- E₀) .* (1 ./ δ) .* Δt .* (1 .- e⁻ᵟ) .+ pool_partition .* e⁻ᵟ
    airborneₜ₊₁ = sum(pool_partition, dims=2)
    Cₜ₊₁ = gcm.C₀ .+ gcm.EtoC .* airborneₜ₊₁
    return Cₜ₊₁, pool_partition, airborneₜ₊₁
end