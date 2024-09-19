struct ARIForcing <: OptionalForcing
    species::Vector{String}
    scaling::AbstractVector{<:Real}
    radiative_efficiency::AbstractVector{<:Real}
    E₀::AbstractVector{<:Real}
    C₀::AbstractVector{<:Real}
    idx_E::AbstractVector{<:Real}
    idx_C::AbstractVector{<:Real}
end



function ARIForcing(path::String, species::Vector{String})
    params = load_ari_forcing_params(path, species)
    return ARIForcing(species, params.scaling, params.radiative_efficiency, params.E₀, params.C₀, params.idx_E, params.idx_C)
end


function computeF(fm::ARIForcing, E, C)
    ΔE = E .- fm.E₀
    ΔC = C .- fm.C₀
    ΔE = ifelse.(isnan.(ΔE), 0.0, ΔE)
    ΔC = ifelse.(isnan.(ΔC), 0.0, ΔC)
    F = (ΔE .+ ΔC) .* fm.radiative_efficiency .* fm.scaling
    return F
end