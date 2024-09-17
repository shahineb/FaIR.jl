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
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    scaling = row.forcing_scale
    C₀ = row.forcing_reference_concentration
    ghg_radiative_efficiency = row.greenhouse_gas_radiative_efficiency
    output = (scaling=scaling, C₀=C₀, ghg_radiative_efficiency=ghg_radiative_efficiency)
    return output
end