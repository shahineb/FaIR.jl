using CSV, DataFrames

function load_csv_emission_data(path::String, species::Vector{String})
    df = CSV.read(path, DataFrame)
    species = lowercase.(species)
    values = Matrix(filter(row -> lowercase(row.name) ∈ species, df)[:, 3:end])
    units = df.units
    year = map(x -> parse(Int, x), names(df)[3:end])
    return (year=year, units=units, values=values)
end