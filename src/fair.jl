module FaIR

export Emissions
export ReservoirModel, EtoC, α
export Meinshausen2020, Leach21, CtoF, compute_ghg_forcing
export BoxModel, ebm_dynamics, FtoT, run


include("Inputs/Inputs.jl")
include("GasCycleModels/GasCycleModels.jl")
include("ForcingModels/ForcingModels.jl")
include("EnergyBalanceModels/EnergyBalanceModels.jl")
using .Inputs
using .GasCycleModels
using .ForcingModels
using .EnergyBalanceModels

end
