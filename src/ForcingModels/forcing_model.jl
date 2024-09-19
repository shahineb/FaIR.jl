abstract type AbstractForcingModel end
abstract type RequiredForcing <: AbstractForcingModel end
abstract type OptionalForcing <: AbstractForcingModel end

struct ForcingModel
    required_forcing::RequiredForcing
    optional_forcing::AbstractVector{<:OptionalForcing}
end


function computeF(fm::ForcingModel; E=E, C=C, Eᶜᵘᵐ=Eᶜᵘᵐ)
    # Dispatching over all the forcing models for different inputs
    # would require types for all kinds of inputs
    # Let's leave that for another time
end


