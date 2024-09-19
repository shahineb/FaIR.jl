abstract type AbstractForcingModel end
struct RequiredForcing <: AbstractForcingModel end
struct OptionalForcing <: AbstractForcingModel end

struct Forcing
    required_forcing::RequiredForcing
    optional_forcing::AbstractVector{<:OptionalForcing}
end


function computeF(fm::Forcing; E=E, C=C, Eᶜᵘᵐ=Eᶜᵘᵐ)
    # Dispatching over all the forcing models for different inputs
    # would require types for all kinds of inputs
    # Let's leave that for another time
end


