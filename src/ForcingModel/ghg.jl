function leach2021ghg(C, C₀, scaling, idxCO₂, idxCH₄, idxN₂O, f₁ᶜᴼ², f₃ᶜᴼ², f₃ᶜᴴ⁴, f₃ᴺ²ᴼ)
    F = fill(NaN, size(C))

    C_CO₂ = C[:, :, :, idxCO₂]
    C_CH₄ = C[:, :, :, idxCH₄]
    C_N₂O = C[:, :, :, idxN₂O]
    C₀_CO₂ = C₀[:, :, :, idxCO₂]
    C₀_CH₄ = C₀[:, :, :, idxCH₄]
    C₀_N₂O = C₀[:, :, :, idxN₂O]

    F[:, :, :, idxCO₂] = scaling .* (f₁ᶜᴼ² .* log.(C_CO₂ ./ C₀_CO₂) .+ f₃ᶜᴼ² .* (.√C_CO₂ .- .√C₀_CO₂))
    F[:, :, :, idxCH₄] = scaling .* f₃ᶜᴴ⁴ .* (.√C_CH₄ .- .√C₀_CH₄)
    F[:, :, :, idxN₂O] = scaling .* f₃ᴺ²ᴼ .* (.√C_N₂O .- .√C₀_N₂O)

    return F
end