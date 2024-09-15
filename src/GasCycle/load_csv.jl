using CSV, DataFrames

function load_gascycle_params(path::String, species::Vector{String})
    speciesdf = CSV.read(path, DataFrame)
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    a = Matrix(row[:, [:a1, :a2, :a3, :a4]])
    τ = Matrix(row[:, [:tau1, :tau2, :tau3, :tau4]])
    r0 = row.r0
    ru = row.rC
    rT = row.rT
    ra = row.rA
    C₀ = row.PI_conc
    EtoC = row.emis2conc
    output = (species=species, a=a, τ=τ, r0=r0, ru=ru, rT=rT, ra=ra, C₀=C₀, EtoC=EtoC, f=f)
    return output
end