struct ACIForcing <: OptionalForcing
    species::Vector{String}
    scaling::AbstractVector{<:Real}
    s::AbstractVector{<:Real}
    β::AbstractVector{<:Real}
    idx_E::AbstractVector{<:Real}
    idx_C::AbstractVector{<:Real}
end



function ACIForcing(path::String, species::Vector{String})
    params = load_aci_forcing_params(path, species)
    return ACIForcing(species, params.scaling, params.s, params.β, params.idx_E, params.idx_C)
end


function computeF(fm::ACIForcing, E, C, E₀, C₀)
    sᴱ = fm.s .* fm.idx_E
    sᶜ = fm.s .* fm.idx_E

    R₀ = log(1 + sum(E₀ .* sᴱ) + sum(C₀ .* sᶜ))
    R = log(1 + sum(E .* sᴱ) + sum(C .* sᶜ))
    F = fm.scaling .* fm.β .* (R - R₀)
    return F
end