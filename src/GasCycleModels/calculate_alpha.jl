
function α(gcm::ReservoirModel, airborneₜ, cumulativeₜ, Tₜ, iirfmax=100.)
    iirf = gcm.r0 .+ gcm.ru .* (cumulativeₜ .- airborneₜ) .+ gcm.rT .* Tₜ .+ gcm.ra .* airborneₜ
    iirf = min.(iirf, iirfmax)
    αₜ = gcm.g₀ .* exp.(iirf ./ gcm.g₁)
    return αₜ
end

function αᶜᴴ⁴(gcm::ReservoirModel, Eₜ, Cₜ, E₀, C₀, Tₜ)
    ΔE = (Eₜ .- E₀) .* gcm.idx_E
    ΔC = (Cₜ .- C₀) .* gcm.idx_C
    ΔE = ifelse.(isnan.(ΔE), 0., ΔE)
    ΔC = ifelse.(isnan.(ΔC), 0., ΔC)
    logα = sum(log.(1 .+ ΔE .* gcm.χ_sensitivity_τᶜᴴ⁴)) +
           sum(log.(1 .+ ΔC .* gcm.χ_sensitivity_τᶜᴴ⁴)) +
           sum(log.(1 .+ Tₜ .*  gcm.T_sensitivity_τᶜᴴ⁴))
    return exp(logα)
end

