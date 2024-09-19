using CSV, DataFrames

function load_and_filter_species(path, species)
    speciesdf = CSV.read(path, DataFrame)
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    return row
end

function load_leach21_params(path, species)
    row = load_and_filter_species(path, species)
    C₀ = row.PI_conc
    f = Matrix(row[:, [:f1, :f2, :f3]])
    output = (species=species, C₀=C₀, f=f)
    return output
end


function load_meinshausen20_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    C₀ = row.forcing_reference_concentration
    radiative_efficiency = row.greenhouse_gas_radiative_efficiency
    output = (scaling=scaling, C₀=C₀, radiative_efficiency=radiative_efficiency)
    return output
end


function load_minor_ghg_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = 0.01 * row.forcing_scale  # 0.01 is a unit handling coefficient, need to be moved elsewhere
    radiative_efficiency = row.greenhouse_gas_radiative_efficiency
    C₀ = row.forcing_reference_concentration
    E₀ = row.forcing_reference_emissions
    return (scaling=scaling, radiative_efficiency=radiative_efficiency, C₀=C₀, E₀=E₀)
end


function load_contrails_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.contrails_radiative_efficiency
    C₀ = row.forcing_reference_concentration
    E₀ = row.forcing_reference_emissions
    return (scaling=scaling, radiative_efficiency=radiative_efficiency, C₀=C₀, E₀=E₀)
end


function load_lapsi_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.lapsi_radiative_efficiency
    C₀ = row.forcing_reference_concentration
    E₀ = row.forcing_reference_emissions
    return (scaling=scaling, radiative_efficiency=radiative_efficiency, C₀=C₀, E₀=E₀)
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
    C₀ = row.forcing_reference_concentration
    E₀ = row.forcing_reference_emissions
    return (scaling=scaling, radiative_efficiency=radiative_efficiency, C₀=C₀, E₀=E₀)
end


function load_aci_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    β = row.aci_scale
    s = row.aci_shape
    C₀ = row.forcing_reference_concentration
    E₀ = row.forcing_reference_emissions
    idx_E = row.aerosol_chemistry_from_emissions
    idx_C = row.aerosol_chemistry_from_concentration
    return (scaling=scaling, β=β, s=s, C₀=C₀, E₀=E₀, idx_E=idx_E, idx_C=idx_C)
end


function load_ari_forcing_params(path, species)
    row = load_and_filter_species(path, species)
    scaling = row.forcing_scale
    radiative_efficiency = row.erfari_radiative_efficiency
    idx_E = row.aerosol_chemistry_from_emissions
    idx_C = row.aerosol_chemistry_from_concentration
    return (scaling=scaling, radiative_efficiency=radiative_efficiency,  C₀=C₀, E₀=E₀, idx_E=idx_E, idx_C=idx_C)
end