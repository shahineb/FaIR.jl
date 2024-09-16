include("load_csv.jl")

struct Meinshausen2020
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

    function Meinshausen2020(a₁::Real,
                             b₁::Real,
                             c₁::Real,
                             d₁::Real,
                             a₂::Real,
                             b₂::Real,
                             c₂::Real,
                             d₂::Real,
                             a₃::Real,
                             b₃::Real,
                             d₃::Real,
                             C₀::AbstractVector{<:Real})
        new(a₁, b₁, c₁, d₁, a₂, b₂, c₂, d₂, a₃, b₃, d₃, C₀)
    end

    function Meinshausen2020(C₀::AbstractVector{<:Real})
        new(-2.4785e-07, 0.00075906, -0.0021492, 5.2488, -0.00034197, 0.00025455, -0.00024357, 0.12173, -8.9603e-05, -0.00012462, 0.045194, C₀)
    end
end


function CtoF(fm::Meinshausen2020, C::AbstractVector{<:Real}, scaling::Real, idxCO₂::Int, idxCH₄::Int, idxN₂O::Int)
    C_CO₂ = C[idxCO₂]
    C_CH₄ = C[idxCH₄]
    C_N₂O = C[idxN₂O]
    C₀_CO₂ = fm.C₀[idxCO₂]
    C₀_CH₄ = fm.C₀[idxCH₄]
    C₀_N₂O = fm.C₀[idxN₂O]
    F = zeros(Real, 3)

    # CO₂
    Cmax_CO₂ = C₀_CO₂ - fm.b₁ / (2 * fm.a₁)
    if C_CO₂ > Cmax_CO₂
        αₚ = fm.d₁ - fm.b₁^2 / (4 * fm.a₁)
    elseif C_CO₂ ≤ C₀_CH₄
        αₚ = fm.d₁
    else
        αₚ = fm.d₁ + fm.a₁ * (C_CO₂ - C₀_CO₂)^2 + f.b₁ * (C_CO₂ - C₀_CO₂)
    end
    αᴺ²ᴼ = fm.c₁ * √C_N₂O
    F[idxCO₂] = scaling * (αₚ + αᴺ²ᴼ) * log(C_CO₂ / C₀_CO₂)

    # CH₄
    F[idxCH₄] = scaling * (fm.a₃ * √C_CH₄ + fm.b₃ * √C_N₂O + fm.d₃) * (√C_CH₄ - √C₀_CH₄)

    # N₂O
    F[idxN₂O] = scaling * (fm.a₂ * √C_CO₂ + fm.b₂ * √C_N₂O + fm.c₂ * √C_CH₄ + fm.d₂) * (√C_N₂O - √C₀_N₂O)
    return F
end