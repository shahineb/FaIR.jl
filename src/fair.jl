module fair

export GasCycleModels
export ReservoirModel

export ForcingModels
export Leach21

export EnergyBalanceModels
export BoxModel

export Inputs
export Emissions

include("EnergyBalanceModels/EnergyBalanceModels.jl")
include("EnergyBalanceModels/box_model.jl")
include("ForcingModels/ForcingModel.jl")
include("ForcingModels/leach21.jl")
include("GasCycleModels/GasCycleModels.jl")
include("GasCycleModels/reservoir_model.jl")
include("Inputs/Inputs.jl")
include("Inputs/emissions.jl")

end
