struct Meinshausen2020 <: RequiredForcing
    a₁::Real
    b₁::Real
    c₁::Real
    d₁::Real
    a₂::Real
    b₂::Real
    c₂::Real
    d₂::Real
    a₃::Real
    b₃::Real
    d₃::Real
    C₀::AbstractVector{<:Real}
    scaling::AbstractVector{<:Real}
    radiative_efficiency::AbstractVector{<:Real}
    active_dims::Vector{Int}
end


function Meinshausen2020(C₀::AbstractVector{<:Real}, scaling::AbstractVector{<:Real}, radiative_efficiency::AbstractVector{<:Real}, active_dims::Vector{Int})
    return Meinshausen2020(-2.4785e-07, 0.00075906, -0.0021492, 5.2488, -0.00034197, 0.00025455, -0.00024357, 0.12173, -8.9603e-05, -0.00012462, 0.045194,
                           C₀, scaling, radiative_efficiency, active_dims)
end 


function Meinshausen2020(path::String, active_dims::Vector{Int})
    species = ["CO2", "CH4", "N2O"]
    params = load_meinshausen20_params(path, species)
    return Meinshausen2020(params.C₀, params.scaling, params.radiative_efficiency, active_dims)
end


function computeF(fm::Meinshausen2020, C)
    C_CO₂, C_CH₄, C_N₂O = C[fm.active_dims]
    C₀_CO₂, C₀_CH₄, C₀_N₂O = fm.C₀
    F = zeros(Real, 3)

    # CO₂
    Cmax_CO₂ = C₀_CO₂ - fm.b₁ / (2 * fm.a₁)
    if C_CO₂ > Cmax_CO₂
        αₚ = fm.d₁ - fm.b₁^2 / (4 * fm.a₁)
    elseif C_CO₂ ≤ C₀_CH₄
        αₚ = fm.d₁
    else
        αₚ = fm.d₁ + fm.a₁ * (C_CO₂ - C₀_CO₂)^2 + fm.b₁ * (C_CO₂ - C₀_CO₂)
    end
    αᴺ²ᴼ = fm.c₁ * √C_N₂O
    F_CO₂ = fm.scaling[1] * (αₚ + αᴺ²ᴼ) * log(C_CO₂ / C₀_CO₂)

    # CH₄
    F_CH₄ = fm.scaling[2] * (fm.a₃ * √C_CH₄ + fm.b₃ * √C_N₂O + fm.d₃) * (√C_CH₄ - √C₀_CH₄)

    # N₂O
    F_N₂O = fm.scaling[3] * (fm.a₂ * √C_CO₂ + fm.b₂ * √C_N₂O + fm.c₂ * √C_CH₄ + fm.d₂) * (√C_N₂O - √C₀_N₂O)

    F = [F_CO₂; F_CH₄; F_N₂O]
    return F
end