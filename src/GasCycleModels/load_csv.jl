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
    return (species=species, a=a, τ=τ, r0=r0, ru=ru, rT=rT, ra=ra, C₀=C₀, EtoC=EtoC)
end


function load_gascycle_params(path, species)
    # speciesdf = CSV.read(path, DataFrame)
    # species = lowercase.(species)
    # row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
    # a = Matrix(row[:, [:partition_fraction0, :partition_fraction1, :partition_fraction2, :partition_fraction3]])
    # τ = Matrix(row[:, [:unperturbed_lifetime0, :unperturbed_lifetime1, :unperturbed_lifetime2, :unperturbed_lifetime3]])
    # χ_sensitivity_τᶜᴴ⁴ = row.ch4_lifetime_chemical_sensitivity
    # T_sensitivity_τᶜᴴ⁴ = row.lifetime_temperature_sensitivity
end


path = "src/defaults/species_configs_properties.csv"
species = ["CO2", "CH4", "N2O"]

speciesdf = CSV.read(path, DataFrame)
species = lowercase.(species)
row = filter(row -> lowercase.(row.name) ∈ species, speciesdf)
a = Matrix(row[:, [:partition_fraction0, :partition_fraction1, :partition_fraction2, :partition_fraction3]])
τ = Matrix(row[:, [:unperturbed_lifetime0, :unperturbed_lifetime1, :unperturbed_lifetime2, :unperturbed_lifetime3]])
χ_sensitivity_τᶜᴴ⁴ = row.ch4_lifetime_chemical_sensitivity
T_sensitivity_τᶜᴴ⁴ = row.lifetime_temperature_sensitivity
r0 = row.iirf_0
ra = row.iirf_airborne
ru = row.iirf_uptake
rT = row.iirf_temperature
g0 = row.g0
g1 = row.g1
C₀ = row.baseline_concentration

row.baseline_concentration
row.forcing_reference_concentration
