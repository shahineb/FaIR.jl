using CSV, DataFrames

function load_leach21_params(path::String, species::Vector{String})
    speciesdf = CSV.read(path, DataFrame)
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    C₀ = row.PI_conc
    f = Matrix(row[:, [:f1, :f2, :f3]])
    output = (species=species, C₀=C₀, f=f)
    return output
end