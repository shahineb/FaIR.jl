module Inputs

export Emissions

include("load_csv.jl")

struct Emissions
    year::Vector{Real}
    species::Vector{String}
    units::Vector{String}
    values::Matrix{Real}
    cumulative::Matrix{Real}
    index::NamedTuple

    function Emissions(year::Vector{Real}, species::Vector{String}, units::Vector{String}, values::Matrix{Real})
        index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
        cumulative = cumsum(values, dims=2)
        new(year, species, units, values, cumulative, index)
    end

    function Emissions(path::String, species::Vector{String})
        emissions = load_csv_emission_data(path, species)
        cumulative = cumsum(emissions.values, dims=2)
        index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
        new(emissions.year, species, emissions.units, emissions.values, cumulative, index)
    end
end

end