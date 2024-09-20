module GasCycleModels

export ReservoirModel, EtoC, α, αᶜᴴ⁴

include("gas_cycle_model.jl")
include("constants.jl")
include("load_csv.jl")
include("reservoir_model.jl")
include("calculate_alpha.jl")

end