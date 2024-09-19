struct Emissions
    year::AbstractVector{<:Real}
    species::Vector{String}
    units::Vector{String}
    values::AbstractMatrix{<:Real}
    cumulative::AbstractMatrix{<:Real}
    index::NamedTuple
end


function Emissions(path, species)
    emissions = load_csv_emission_data(path, species)
    index = NamedTuple{Tuple(Symbol.(species))}(1:length(species))
    cumulative = cumsum(emissions.values, dims=2)
    return Emissions(emissions.year, species, emissions.units, emissions.values, cumulative, index)
end