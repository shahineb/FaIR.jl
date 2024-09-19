struct ACIForcing <: OptionalForcing
    species::Vector{String}
    scaling::AbstractVector{<:Real}
    s::AbstractVector{<:Real}  # sensitivity = shape
    β::AbstractVector{<:Real}  # scale
    E₀::AbstractVector{<:Real}
    C₀::AbstractVector{<:Real}
    idx_E::AbstractVector{<:Real}
    idx_C::AbstractVector{<:Real}
end



function ACIForcing(path::String, species::Vector{String})
    params = load_aci_forcing_params(path, species)
    return ACIForcing(species, params.scaling, params.s, params.β, params.E₀, params.C₀, params.idx_E, params.idx_C)
end


function computeF(fm::ACIForcing, E, C)
    sᴱ = fm.s .* fm.idx_E
    sᶜ = fm.s .* fm.idx_C

    E₀ = ifelse.(isnan.(fm.E₀), 0.0, fm.E₀)
    E = ifelse.(isnan.(E), 0.0, E)
    C₀ = ifelse.(isnan.(fm.C₀), 0.0, fm.C₀)
    C = ifelse.(isnan.(C), 0.0, C)

    R₀ = log.(1 .+ E₀ .* sᴱ .+ C₀ .* sᶜ)
    R = log.(1 .+ E .* sᴱ .+ C .* sᶜ)
    F = fm.scaling .* fm.β .* (R - R₀)
    return F
end