include("GasCycle/GasCycle.jl")
include("ForcingModel/ForcingModel.jl")
include("EnergyBalanceModel/EnergyBalanceModel.jl")
include("Inputs/Inputs.jl")

using CSV, DataFrames
using .GasCycle, .ForcingModel, .EnergyBalanceModel, .Inputs

E = Emissions("src/defaults/historical-emissions.csv", ["CO2", "CH4", "SO2", "BC"])
gascycle = GasCycleModel("src/defaults/gas_parameters.csv", ["CO2", "CH4", "SO2", "BC"])

n_pool = 4
Eₜ = E.data[:, 1]
Burdenₜ₋₁ = zeros(4, n_pool)
αₜ = ones(4)
E₀ = E.data[:, 1]
Δt = 1.

Cₜ, Burdenₜ = EtoC(gascycle, Eₜ, Burdenₜ₋₁, αₜ, E₀, Δt)
