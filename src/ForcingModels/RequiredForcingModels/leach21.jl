struct Leach21 <: RequiredForcing
    f₁ᶜᴼ²::Real
    f₃ᶜᴼ²::Real
    f₃ᶜᴴ⁴::Real
    f₃ᴺ²ᴼ::Real
    C₀::AbstractVector{<:Real}
    scaling::AbstractVector{<:Real}
end


function Leach21(path::String, species::Vector{String})
    params = load_leach21_params(path, species)
    index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
    scaling = ones(length(species))
    return Leach21(params.f[index.CO2, 1], params.f[index.CO2, 3], params.f[index.CH4, 3], params.f[index.N2O, 3], params.C₀, scaling)
end

function computeF(fm::Leach21, C, idxCO₂, idxCH₄, idxN₂O)
    # TODO : remove indexing from arguments
    C_CO₂ = C[idxCO₂]
    C_CH₄ = C[idxCH₄]
    C_N₂O = C[idxN₂O]
    C₀_CO₂ = fm.C₀[idxCO₂]
    C₀_CH₄ = fm.C₀[idxCH₄]
    C₀_N₂O = fm.C₀[idxN₂O]

    F = zeros(Real, 3)
    F[idxCO₂] = fm.scaling[idxCO₂] .* (fm.f₁ᶜᴼ² .* log.(C_CO₂ ./ C₀_CO₂) .+ fm.f₃ᶜᴼ² .* (.√C_CO₂ .- .√C₀_CO₂))
    F[idxCH₄] = fm.scaling[idxCH₄] .* fm.f₃ᶜᴴ⁴ .* (.√C_CH₄ .- .√C₀_CH₄)
    F[idxN₂O] = fm.scaling[idxN₂O] .* fm.f₃ᴺ²ᴼ .* (.√C_N₂O .- .√C₀_N₂O)
    return F
end