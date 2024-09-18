
struct LinearForcing{T}
    species::Vector{String}
    scaling::AbstractVector{<:Real}
    radiative_efficiency::AbstractVector{<:Real}
end

const GHGForcing = LinearForcing{:GHG}
const ContrailsForcing = LinearForcing{:Contrails} 
const LAPSIForcing = LinearForcing{:LAPSI}
const StratosphericVapourForcing = LinearForcing{:StratosphericVapour}
const LandUseForcing = LinearForcing{:LandUse}



function GHGForcing(path::String, species::Vector{String})
    params = load_ghg_linear_forcing_params(path, species)
    return GHGForcing(species, params.scaling, params.radiative_efficiency)
end


# path = "src/defaults/species_configs_properties.csv"
# species = ["SO2", "CFC-11", "CO", "CH3Cl"]


function compute_linear_forcing(fm::LinearForcing, driver, driver₀)
    return ((driver .- driver₀) .* fm.radiative_efficiency) .* fm.scaling
end