using CSV, DataFrames

function load_leach21_params(path, species)
    speciesdf = CSV.read(path, DataFrame)
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    C₀ = row.PI_conc
    f = Matrix(row[:, [:f1, :f2, :f3]])
    output = (species=species, C₀=C₀, f=f)
    return output
end


function load_meinshausen20_params(path, species)
    speciesdf = CSV.read(path, DataFrame)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    scaling = row.forcing_scale
    C₀ = row.forcing_reference_concentration
    ghg_radiative_efficiency = row.greenhouse_gas_radiative_efficiency
    output = (scaling=scaling, C₀=C₀, ghg_radiative_efficiency=ghg_radiative_efficiency)
    return output
end

# path = "src/defaults/species_configs_properties.csv"
# species = lowercase.(["CO2", "CH4", "N2O"])
# load_meinshausen20_params(path, species)
# speciesdf = CSV.read(path, DataFrame)
# row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
# row


# a = Matrix(row[:, [:partition_fraction0, :partition_fraction1, :partition_fraction2, :partition_fraction3]])
# τ = Matrix(row[:, [:unperturbed_lifetime0, :unperturbed_lifetime1, :unperturbed_lifetime2, :unperturbed_lifetime3]])
# χ_sensitivity_τᶜᴴ⁴ = row.ch4_lifetime_chemical_sensitivity
# T_sensitivity_τᶜᴴ⁴ = row.lifetime_temperature_sensitivity