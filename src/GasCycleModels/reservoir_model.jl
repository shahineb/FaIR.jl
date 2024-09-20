struct ReservoirModel <: AbstractGasCycleModel
    species::Vector{String}
    a::AbstractMatrix{<:Real}
    τ::AbstractMatrix{<:Real}
    r0::AbstractVector{<:Real}
    ru::AbstractVector{<:Real}
    rT::AbstractVector{<:Real}
    ra::AbstractVector{<:Real}
    C₀::AbstractVector{<:Real}
    CperE::AbstractVector{<:Real}
    g₀::AbstractVector{<:Real}
    g₁::AbstractVector{<:Real}
    molecular_weight::AbstractVector{<:Real}
    χ_sensitivity_τᶜᴴ⁴::AbstractVector{<:Real}
    T_sensitivity_τᶜᴴ⁴::AbstractVector{<:Real}
    idx_E::AbstractVector{<:Real}
    idx_C::AbstractVector{<:Real}
end


function ReservoirModel(species::Vector{String},
                        a::AbstractMatrix{<:Real},
                        τ::AbstractMatrix{<:Real},
                        r0::AbstractVector{<:Real},
                        ru::AbstractVector{<:Real},
                        rT::AbstractVector{<:Real},
                        ra::AbstractVector{<:Real},
                        C₀::AbstractVector{<:Real},
                        molecular_weight::AbstractVector{<:Real},
                        χ_sensitivity_τᶜᴴ⁴::AbstractVector{<:Real},
                        T_sensitivity_τᶜᴴ⁴::AbstractVector{<:Real},
                        idx_E::AbstractVector{<:Real},
                        idx_C::AbstractVector{<:Real})
    g₀, g₁ = compute_g(a, τ)
    CperE = compute_CperE(molecular_weight)
    return ReservoirModel(species, a, τ, r0, ru, rT, ra, C₀, CperE, g₀, g₁, molecular_weight, χ_sensitivity_τᶜᴴ⁴, T_sensitivity_τᶜᴴ⁴, idx_E, idx_C)
end


function ReservoirModel(path::String, species::Vector{String})
    params = load_gascycle_params(path, species)
    return ReservoirModel(species,
                          params.a,
                          params.τ,
                          params.r0,
                          params.ru,
                          params.rT,
                          params.ra,
                          params.C₀,
                          params.molecular_weight,
                          params.χ_sensitivity_τᶜᴴ⁴,
                          params.T_sensitivity_τᶜᴴ⁴,
                          params.idx_E,
                          params.idx_C)
end



function compute_g(a, τ)
    δ = 100 ./ τ
    e⁻ᵟ = exp.(-δ)
    g₁ = sum(a .* τ .* (1 .- (1 .+ δ) .* e⁻ᵟ), dims=2)
    g₀ = exp.(-sum(a .* τ .* (1 .- e⁻ᵟ), dims=2) ./ g₁)
    return vec(g₀), vec(g₁)
end


function compute_CperE(molecular_weight)
    CperE = 1 ./ (MASS_ATMOSPHERE ./ 1e18 .* molecular_weight ./ MOLECULAR_WEIGHT_AIR)
    return CperE
end


function EtoC(gcm::ReservoirModel, Eₜ, pool_partition, αₜ, E₀, Δt)
    δ = Δt ./ (αₜ .* gcm.τ)
    e⁻ᵟ = exp.(-δ)
    pool_partition = gcm.a .* (Eₜ .- E₀) .* (1 ./ δ) .* Δt .* (1 .- e⁻ᵟ) .+ pool_partition .* e⁻ᵟ
    airborneₜ₊₁ = sum(pool_partition, dims=2)
    Cₜ₊₁ = gcm.C₀ .+ gcm.CperE .* airborneₜ₊₁
    return Cₜ₊₁, pool_partition, airborneₜ₊₁
end