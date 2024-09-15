include("GasCycle/GasCycle.jl")
include("ForcingModel/ForcingModel.jl")
include("EnergyBalanceModel/EnergyBalanceModel.jl")
include("Inputs/Inputs.jl")

using .GasCycle, .ForcingModel, .EnergyBalanceModel, .Inputs

species = ["CO2", "CH4"]
E = Emissions("src/defaults/ssp245-emissions.csv", species)
gascycle = GasCycleModel("src/defaults/gas_parameters.csv", species)

n_species = length(species)
n_pool = 4
pool_partition = zeros(Real, n_species, n_pool)
# E₀ = [0., 0., 2.44004844, 2.09777075]
E₀ = [0., 0.]
Δt = 1.
iirfmax = 100.
airborneₜ = zeros(Real, n_species)
Tₜ = 1.
αs = zeros(Real, size(E.values))
C = zeros(Real, size(E.values))
F = zeros(Real, size(E.values))
for t in 1:length(E.year)
    αs[:, t] = α(gascycle, airborneₜ, E.cumulative[:, t], Tₜ, iirfmax)
    C[:, t], pool_partition = EtoC(gascycle, E.values[:, t], pool_partition, αs[:, t], E₀, Δt)
    # F[:, t] = NaN
end


using GLMakie
fig = Figure()
units = ["ppm", "ppb"]
for i in 1:n_species
    ax = Axis(fig[1, i]; xlabel = "Year", ylabel = units[i],
              title = species[i])
    
    # Plot the concentration over time
    lines!(ax, E.year, C[i, :], label=species[i])
    
    # Enable legend
    axislegend(ax, position=:rb)
end
display(fig)
GLMakie.closeall()