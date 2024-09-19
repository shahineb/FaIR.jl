module ForcingModels

# Required forcing models
export Meinshausen2020,
       Leach21

# Optional forcing models
export MinorGHGForcing,
       ContrailsForcing,
       LAPSIForcing,
       StratosphericVapourForcing,
       LandUseForcing,
       ACIForcing,
       ARIForcing,
       VolcanicForcing

# Dispath of forcing generic computation function
export computeF


include("forcing_model.jl")
include("load_csv.jl")
include("RequiredForcingModels/meinshausen2020.jl")
include("RequiredForcingModels/leach21.jl")
include("OptionalForcingModels/linear_forcing.jl")
include("OptionalForcingModels/ari_forcing.jl")
include("OptionalForcingModels/aci_forcing.jl")
include("PrescribedForcings/volcanic_forcing.jl")

end