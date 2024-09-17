function compute_aci_forcing(E, C, E₀, C₀, scaling, β, s, idx_E, idx_C)
    sᴱ = s[idx_from_E]
    sᶜ = s[idx_from_C]

    R₀ = log(1 + sum(E₀[idx_E] .* sᴱ) + sum(C₀[idx_C] .* sᶜ))
    R = log(1 + sum(E[idx_E] .* sᴱ) + sum(C[idx_C] .* sᶜ))
    F = scaling * β * (R - R₀)
    return F
end

