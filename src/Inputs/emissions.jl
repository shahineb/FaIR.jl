struct Emissions
    year::AbstractVector{<:Real}
    species::Vector{String}
    units::Vector{String}
    values::AbstractMatrix{<:Real}
    cumulative::AbstractMatrix{<:Real}
    index::NamedTuple
end


function Emissions(path::String, species::Vector{String})
    emissions = load_csv_emission_data(path, species)
    cumulative = cumsum(emissions.values, dims=2)
    index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
    Emissions(emissions.year, species, emissions.units, emissions.values, cumulative, index)
end