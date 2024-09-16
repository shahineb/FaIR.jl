module GasCycleModels

export ReservoirModel, EtoC, α, αᶜᴴ⁴

include("load_csv.jl")
include("reservoir_model.jl")
include("calculate_alpha.jl")

end