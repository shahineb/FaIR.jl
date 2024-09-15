include("GasCycle/GasCycle.jl")
include("ForcingModel/ForcingModel.jl")
include("EnergyBalanceModel/EnergyBalanceModel.jl")
include("Inputs/Inputs.jl")

using .GasCycle, .ForcingModel, .EnergyBalanceModel, .Inputs

species = ["CO2", "CH4", "SO2", "BC"]
E = Emissions("src/defaults/ssp245-emissions.csv", species)
gascycle = GasCycleModel("src/defaults/gas_parameters.csv", species)


n_pool = 4
pool_partition = zeros(Real, 4, n_pool)
E₀ = [0., 0., 2.44004844, 2.09777075]
Δt = 1.
iirfmax = 100.
airborneₜ = zeros(Real, 4)
Tₜ = 0
αs = zeros(Real, size(E.values))
C = zeros(Real, size(E.values))
for t in 1:length(E.year)
    αs[:, t] = α(gascycle, airborneₜ, E.cumulative[:, t], Tₜ, iirfmax)
    C[:, t], pool_partition = EtoC(gascycle, E.values[:, t], pool_partition, αs[:, t], E₀, Δt)
end

# Cₜ, Burdenₜ = EtoC(gascycle, Eₜ, Burdenₜ₋₁, αₜ, E₀, Δt)
C
αₜ