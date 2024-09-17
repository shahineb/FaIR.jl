function compute_linear_forcing(driver, driver₀, scaling, radiative_efficiency)
    return ((driver .- driver₀) .* radiative_efficiency) .* scaling
end