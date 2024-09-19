using CSV, DataFrames


function load_and_filter_species(path, species)
    speciesdf = CSV.read(path, DataFrame)
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    return row
end


function load_minor_ghg_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = 0.01 * row.forcing_scale  # 0.01 is a unit handling coefficient, need to be moved elsewhere
    radiative_efficiency = row.greenhouse_gas_radiative_efficiency
    return (scaling=scaling, radiative_efficiency=radiative_efficiency)
end


function load_contrails_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.contrails_radiative_efficiency
    return (scaling=scaling, radiative_efficiency=radiative_efficiency)
end


function load_lapsi_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.lapsi_radiative_efficiency
    return (scaling=scaling, radiative_efficiency=radiative_efficiency)
end


function load_stratospheric_vapour_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.h2o_stratospheric_factor
    return (scaling=scaling, radiative_efficiency=radiative_efficiency)
end


function load_landuse_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.land_use_cumulative_emissions_to_forcing
    return (scaling=scaling, radiative_efficiency=radiative_efficiency)
end


function load_aci_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    β = row.aci_scale
    s = row.aci_shape
    idx_E = row.aerosol_chemistry_from_emissions
    idx_C = row.aerosol_chemistry_from_concentration
    return (scaling=scaling, β=β, s=s, idx_E=idx_E, idx_C=idx_C)
end


function load_ari_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.erfari_radiative_efficiency
    idx_E = row.aerosol_chemistry_from_emissions
    idx_C = row.aerosol_chemistry_from_concentration
    return (scaling=scaling, radiative_efficiency=radiative_efficiency, idx_E=idx_E, idx_C=idx_C)
end


# load_ghg_linear_forcing_params("src/defaults/species_configs_properties.csv", ["SO2", "CFC-11", "CO", "CH3Cl"])
# speciesdf = CSV.read("src/defaults/species_configs_properties.csv", DataFrame)

# row = filter(row -> row.greenhouse_gas == 1, speciesdf)

# println(row.name)

# species = lowercase.(["SO2", "CFC-11", "CO", "CH3Cl"])
# row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)


