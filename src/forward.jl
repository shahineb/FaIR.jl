include("GasCycle/GasCycle.jl")
include("ForcingModel/ForcingModel.jl")
include("EnergyBalanceModel/EnergyBalanceModel.jl")
include("Inputs/Inputs.jl")

using .GasCycle, .ForcingModel, .EnergyBalanceModel, .Inputs

species = ["CO2", "CH4", "N2O"]
E = Emissions("src/defaults/ssp245-emissions.csv", species)

Nₜ = length(E.year)
Δt = 1.
seed = 2

gascycle = GasCycleModel("src/defaults/gas_parameters.csv", species)
fm = Leach21("src/defaults/gas_parameters.csv", species)

ebm = EBM("src/defaults/4xCO2_cummins_ebm3.csv", seed, Δt, Nₜ)
A = computeA(ebm)
eᴬ = exp(A)
bd = compute_bd(ebm)
wd = samplevariability(ebm)



n_species = length(species)
n_pool = 4
pool_partition = zeros(Real, n_species, n_pool)
# E₀ = [0., 0., 2.44004844, 2.09777075]
E₀ = [0., 0., 0.]
iirfmax = 100.
airborneₜ = zeros(Real, n_species)
αs = zeros(Real, size(E.values))
C = zeros(Real, size(E.values))
F = zeros(Real, size(E.values))
T = zeros(Real, ebm.Nₜ, ebm.Nbox + 1)
for t in 2:Nₜ
    αs[:, t - 1] = α(gascycle, airborneₜ, E.cumulative[:, t - 1], T[t - 1, 2], iirfmax)
    C[:, t], pool_partition = EtoC(gascycle, E.values[:, t - 1], pool_partition, αs[:, t - 1], E₀, Δt)
    F[:, t] = CtoF(fm, C[:, t], 1., E.index.CO2, E.index.CH4, E.index.N2O)
    T[t, :] = FtoT(T[t - 1, :], eᴬ, bd, wd[t - 1, :], sum(F[:, t]))
end


using GLMakie

# Concentrations plot
fig = Figure()
units = ["ppm", "ppb", "ppb"]
for i in 1:n_species
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
    axislegend(ax, position=:rb)
end
ax = Axis(fig[1, n_species + 1]; xlabel = "Year", ylabel = "Wm⁻²", title = "Total")
lines!(ax, E.year, vec(sum(F, dims=1)), label="Total")
axislegend(ax, position=:rb)
display(fig)


# Temperature plot
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Year", ylabel = "ΔT (K)", title = "GMST")
lines!(ax, E.year, T[:, 2])
display(fig)

GLMakie.closeall()