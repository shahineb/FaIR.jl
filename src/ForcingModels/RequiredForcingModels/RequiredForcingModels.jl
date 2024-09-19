module RequiredForcingModels

export RequiredForcing, Meinshausen2020, Leach21, compute_forcing

include("load_csv.jl")
include("required_forcing.jl")
include("leach21.jl")
include("meinshausen2020.jl")

end