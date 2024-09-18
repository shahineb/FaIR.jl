using CSV, DataFrames


function load_and_filter_species(path, species)
    speciesdf = CSV.read(path, DataFrame)
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    return row
end


function load_ghg_linear_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.greenhouse_gas_radiative_efficiency
    return (scaling=scaling, radiative_efficiency=radiative_efficiency)
end




# load_ghg_linear_forcing_params("src/defaults/species_configs_properties.csv", ["SO2", "CFC-11", "CO", "CH3Cl"])
# speciesdf = CSV.read("src/defaults/species_configs_properties.csv", DataFrame)
# species = lowercase.(["SO2", "CFC-11", "CO", "CH3Cl"])
# row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)


