using CSV, DataFrames

function load_gascycle_params(path, species)
    speciesdf = CSV.read(path, DataFrame)
    species = lowercase.(species)
    row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    a = Matrix(row[:, [:partition_fraction0, :partition_fraction1, :partition_fraction2, :partition_fraction3]])
    τ = Matrix(row[:, [:unperturbed_lifetime0, :unperturbed_lifetime1, :unperturbed_lifetime2, :unperturbed_lifetime3]])
    r0 = row.iirf_0
    ra = row.iirf_airborne
    ru = row.iirf_uptake
    rT = row.iirf_temperature
    g₀ = row.g0
    g₁ = row.g1
    C₀ = row.baseline_concentration
    molecular_weight = row.molecular_weight
    χ_sensitivity_τᶜᴴ⁴ = row.ch4_lifetime_chemical_sensitivity
    T_sensitivity_τᶜᴴ⁴ = row.lifetime_temperature_sensitivity
    output = (a=a,
              τ=τ,
              r0=r0,
              ra=ra,
              ru=ru,
              rT=rT,
              g₀=g₀,
              g₁=g₁,
              C₀=C₀,
              molecular_weight=molecular_weight,
              χ_sensitivity_τᶜᴴ⁴=χ_sensitivity_τᶜᴴ⁴,
              T_sensitivity_τᶜᴴ⁴=T_sensitivity_τᶜᴴ⁴)
    return output
end