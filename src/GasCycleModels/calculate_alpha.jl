
function α(gcm::ReservoirModel, airborneₜ, cumulativeₜ, Tₜ, iirfmax=100.)
    iirf = gcm.r0 .+ gcm.ru .* (cumulativeₜ .- airborneₜ) .+ gcm.rT .* Tₜ .+ gcm.ra .* airborneₜ
    iirf = min.(iirf, iirfmax)
    αₜ = gcm.g₀ .* exp.(iirf ./ gcm.g₁)
    return αₜ
end

function αᶜᴴ⁴(gcm::ReservoirModel, Eₜᶜᴴ⁴, Cₜᶜᴴ⁴, Tₜ, E₀ᶜᴴ⁴, C₀ᶜᴴ⁴)
    logα = log(1 + (Eₜᶜᴴ⁴ - E₀ᶜᴴ⁴) * gcm.χ_sensitivity_τᶜᴴ⁴) +
           log(1 + (Cₜᶜᴴ⁴ - C₀ᶜᴴ⁴) * gcm.χ_sensitivity_τᶜᴴ⁴) +
           log(1 + Tₜ *  gcm.T_sensitivity_τᶜᴴ⁴)
    return exp(logα)
end

