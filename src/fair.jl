module FaIR

export Emissions
export ReservoirModel, EtoC, α
export Leach21, CtoF
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
