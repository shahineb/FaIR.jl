module Inputs

export Emissions

include("load_csv.jl")

struct Emissions
    year::Vector{Real}
    species::Vector{String}
    units::Vector{String}
    data::Matrix{Real}
    index::NamedTuple

    function Emissions(year::Vector{Real}, species::Vector{String}, units::Vector{String}, data::Matrix{Real})
        index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
        new(year, species, units, data, index)
    end

    function Emissions(path::String, species::Vector{String})
        emissions = load_csv_emission_data(path, species)
        index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
        new(emissions.year, species, emissions.units, emissions.data, index)
    end
end

end