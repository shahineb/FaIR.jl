using FaIR
using CSV, DataFrames
using GLMakie

# Path to parameter files
emission_csv = "src/defaults/ssp245-emissions.csv"
species_csv = "src/defaults/species_configs_properties.csv"
ebm_csv = "src/defaults/4xCO2_cummins_ebm3.csv"

# Define emissions
species = ["CO2", "CH4", "N2O", "SO2", "BC"]
E = Emissions(emission_csv, species)

# Define gas cycle model
gas_model = ReservoirModel(species_csv, species)

# Define forcing models
CO₂_CH₄_N₂O_forcing_model = Meinshausen2020(species_csv, ["CO2", "CH4", "N2O"])
ari_forcing = ARIForcing(species_csv, ["CH4", "N2O", "SO2", "BC"])
aci_forcing = ACIForcing(species_csv, ["SO2", "BC"])

volcanic_df = CSV.read("src/defaults/volcanic_ERF_monthly_175001-201912.csv", DataFrame)
v = volcanic_df.erf
reshaped_v = reshape(v, 12, :)
block_means = sum(reshaped_v, dims=1) ./ 12
volcanic_forcing = VolcanicForcing([vec(block_means); zeros(351 - 270)])


# Define energy balance model
seed = 2
Δt = 1.
Nₜ = length(E.year)
ebm = BoxModel(ebm_csv, seed, Δt, Nₜ)
eᴬ, bd, wd = ebm_dynamics(ebm)

# Run model
n_species = length(species)
n_pool = 4
pool_partition = zeros(Real, n_species, n_pool)
E₀ = zeros(length(species))
iirfmax = 100.
airborneₜ = zeros(Real, n_species)
αs = zeros(Real, size(E.values))
C = zeros(Real, size(E.values))
C[:, 1] = gas_model.C₀
F = zeros(Real, size(E.values))
T = zeros(Real, ebm.Nₜ, ebm.Nbox + 1)


for t in 2:Nₜ
    αs[:, t - 1] = α(gas_model, airborneₜ, E.cumulative[:, t - 1], T[t - 1, 2], iirfmax)
    C[:, t], pool_partition = EtoC(gas_model, E.values[:, t - 1], pool_partition, αs[:, t - 1], E₀, Δt)
    FCO₂_CH₄_N₂O = computeF(CO₂_CH₄_N₂O_forcing_model, C[:, t], E.index.CO2, E.index.CH4, E.index.N2O)
    Fari = computeF(ari_forcing, E.values[2:end, t], C[2:end, t])
    Faci = computeF(aci_forcing, E.values[4:end, t], C[4:end, t])
    F[:, t] = [FCO₂_CH₄_N₂O; 0; 0] .+ [0.; Fari] .+ [0.; 0.; 0; Faci]
    T[t, :] = FtoT(T[t - 1, :], eᴬ, bd, wd[t - 1, :], sum(F[:, t]) + volcanic_forcing.values[t])
end


# Emissions plot
fig = Figure()
units = ["GtC/yr", "Mt CH₄/yr", "Mt N₂O/yr", "Mt SO₂/yr", "Mt BC/yr"]
for i in 1:n_species
    ax = Axis(fig[1, i]; xlabel = "Year", ylabel = units[i], title = species[i])
    
    # Plot the emissions over time
    lines!(ax, E.year, E.values[i, :], label=species[i])
    
    # Enable legend
    axislegend(ax, position=:rb)
end
display(fig)


# Concentrations plot
fig = Figure()
units = ["ppm", "ppb", "ppb"]
for i in 1:n_species-2
    ax = Axis(fig[1, i]; xlabel = "Year", ylabel = units[i], title = species[i])
    
    # Plot the concentration over time
    lines!(ax, E.year, C[i, :], label=species[i])
    
    # Enable legend
    axislegend(ax, position=:rb)
end
display(fig)


# Forcing plot
fig = Figure()
for i in 1:n_species
    ax = Axis(fig[1, i]; xlabel = "Year", ylabel = "Wm⁻²", title = species[i])
    
    # Plot the forcing over time
    lines!(ax, E.year, F[i, :], label=species[i])
    
    # Enable legend
    axislegend(ax, position=:lt)
end
ax = Axis(fig[1, n_species + 1]; xlabel = "Year", ylabel = "Wm⁻²", title = "Total")
lines!(ax, E.year, vec(sum(F, dims=1)) + volcanic_forcing.values, label="Total")
axislegend(ax, position=:rb)
display(fig)


# Temperature plot
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Year", ylabel = "ΔT (K)", title = "GMST")
lines!(ax, E.year, T[:, 2])
display(fig)

GLMakie.closeall()