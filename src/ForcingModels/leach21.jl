include("load_csv.jl")

struct Leach21
    f₁ᶜᴼ²::Real
    f₃ᶜᴼ²::Real
    f₃ᶜᴴ⁴::Real
    f₃ᴺ²ᴼ::Real
    C₀::AbstractVector{<:Real}

    function Leach21(f₁ᶜᴼ²::Real, f₃ᶜᴼ²::Real, f₃ᶜᴴ⁴::Real, f₃ᴺ²ᴼ::Real, C₀::AbstractVector{<:Real})
        new(f₁ᶜᴼ², f₃ᶜᴼ², f₃ᶜᴴ⁴, f₃ᴺ²ᴼ, C₀)
    end

    function Leach21(path::String, species::Vector{String})
        params = load_leach21_params(path, species)
        index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
        new(params.f[index.CO2, 1], params.f[index.CO2, 3], params.f[index.CH4, 3], params.f[index.N2O, 3], params.C₀)
    end
end


function CtoF(fm::Leach21, C::AbstractVector{<:Real}, scaling::Real, idxCO₂::Int, idxCH₄::Int, idxN₂O::Int)
    C_CO₂ = C[idxCO₂]
    C_CH₄ = C[idxCH₄]
    C_N₂O = C[idxN₂O]
    C₀_CO₂ = fm.C₀[idxCO₂]
    C₀_CH₄ = fm.C₀[idxCH₄]
    C₀_N₂O = fm.C₀[idxN₂O]

    F = zeros(Real, 3)
    F[idxCO₂] = scaling .* (fm.f₁ᶜᴼ² .* log.(C_CO₂ ./ C₀_CO₂) .+ fm.f₃ᶜᴼ² .* (.√C_CO₂ .- .√C₀_CO₂))
    F[idxCH₄] = scaling .* fm.f₃ᶜᴴ⁴ .* (.√C_CH₄ .- .√C₀_CH₄)
    F[idxN₂O] = scaling .* fm.f₃ᴺ²ᴼ .* (.√C_N₂O .- .√C₀_N₂O)
    return F
end