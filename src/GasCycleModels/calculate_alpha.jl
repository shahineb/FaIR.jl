
function α(gcm::ReservoirModel, airborneₜ::AbstractVector{<:Real}, cumulativeₜ::AbstractVector{<:Real}, Tₜ::Real, iirfmax::Real)
    iirf = gcm.r0 .+ gcm.ru .* (cumulativeₜ .- airborneₜ) .+ gcm.rT .* Tₜ .+ gcm.ra .* airborneₜ
    iirf = min.(iirf, iirfmax)
    αₜ = gcm.g₀ .* exp.(iirf ./ gcm.g₁)
    return αₜ
end

function αᶜᴴ⁴(gcm::ReservoirModel, Eₜᶜᴴ⁴::Real, Cₜᶜᴴ⁴::Real, Tₜ::Real, E₀ᶜᴴ⁴::Real, C₀ᶜᴴ⁴::Real,)
    logα = log(1 + (Eₜᶜᴴ⁴ - E₀ᶜᴴ⁴) * gcm.χ_sensitivity_τᶜᴴ⁴) +
           log(1 + (Cₜᶜᴴ⁴ - C₀ᶜᴴ⁴) * gcm.χ_sensitivity_τᶜᴴ⁴) +
           log(1 + Tₜ *  gcm.T_sensitivity_τᶜᴴ⁴)
    return exp(logα)
end

