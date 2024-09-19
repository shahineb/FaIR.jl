struct ARIForcing
    species::Vector{String}
    scaling::AbstractVector{<:Real}
    radiative_efficiency::AbstractVector{<:Real}
    idx_E::AbstractVector{<:Real}
    idx_C::AbstractVector{<:Real}
end



function ARIForcing(path::String, species::Vector{String})
    params = load_ari_forcing_params(path, species)
    return ARIForcing(species, params.scaling, params.radiative_efficiency, params.idx_E, params.idx_C)
end


function compute_ari_forcing(fm::ARIForcing, E, C, E₀, C₀)
    ΔE = E .- E₀
    ΔC = C .- C₀
    F = (ΔE .+ ΔC) .* fm.radiative_efficiency .* fm.scaling
    return F
end


# Could even load the E0 and C0 into the struct, that would make a lot of sense
# Also would ideally want to feed in directly the whole emissions/concentration vector, need to find a swifty way to do that