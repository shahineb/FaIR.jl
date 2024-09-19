module ForcingModels

export Meinshausen2020,
       Leach21

export MinorGHGForcing,
       ContrailsForcing,
       LAPSIForcing,
       StratosphericVapourForcing,
       LandUseForcing,
       ACIForcing,
       ARIForcing

export compute_forcing


include("RequiredForcingModels/RequiredForcingModels.jl")
using .RequiredForcingModels

include("OptionalForcingModels/OptionalForcingModels.jl")
using .OptionalForcingModels


function compute_forcing(fm::RequiredForcing, args...)
    return RequiredForcingModels.compute_forcing(fm, args...)
end

function compute_forcing(fm::OptionalForcing, args...)
    return OptionalForcingModels.compute_forcing(fm, args...)
end

end