struct LinearForcing{T} <: OptionalForcing
    species::Vector{String}
    scaling::AbstractVector{<:Real}
    radiative_efficiency::AbstractVector{<:Real}
    E₀::AbstractVector{<:Real}
    C₀::AbstractVector{<:Real}
end


const MinorGHGForcing = LinearForcing{:GHG}
const ContrailsForcing = LinearForcing{:Contrails}
const LAPSIForcing = LinearForcing{:LAPSI}
const StratosphericVapourForcing = LinearForcing{:StratosphericVapour}
const LandUseForcing = LinearForcing{:LandUse}



function loading_method(T)
    if T == :GHG
        return load_minor_ghg_forcing_params
    elseif T == :Contrails
        return load_contrails_forcing_params
    elseif T == :LAPSI
        return load_lapsi_forcing_params
    elseif T == :StratosphericVapour
        return load_stratospheric_vapour_forcing_params
    elseif T == :LandUse
        return load_landuse_forcing_params
    else
        throw(ArgumentError("Unknown forcing type."))
    end
end


function LinearForcing{T}(path::String, species::Vector{String}) where T
    params = loading_method(T)(path, species)
    return LinearForcing{T}(species, params.scaling, params.radiative_efficiency, params.E₀, params.C₀)
end



function computeF(fm::LinearForcing, driver, driver₀)
    return ((driver .- driver₀) .* fm.radiative_efficiency) .* fm.scaling
end
computeF(fm::MinorGHGForcing, C) = computeF(fm, C, fm.C₀)
computeF(fm::ContrailsForcing, E) = computeF(fm, E, zeros(size(E)))
computeF(fm::LAPSIForcing, E) = computeF(fm, E, fm.E₀)
computeF(fm::StratosphericVapourForcing, C) = computeF(fm, C, fm.C₀)
computeF(fm::LandUseForcing, Eᶜᵘᵐ) = computeF(fm, Eᶜᵘᵐ, zeros(size(Eᶜᵘᵐ)))