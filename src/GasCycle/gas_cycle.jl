include("load_csv.jl")


struct GasCycleModel
    species::Vector{String}
    a::Matrix{Real}
    τ::Matrix{Real}
    r0::Vector{Real}
    ru::Vector{Real}
    rT::Vector{Real}
    ra::Vector{Real}
    C₀::Vector{Real}
    EtoC::Vector{Real}
    f::Matrix{Real}
    g₀::Vector{Real}
    g₁::Vector{Real}

    function GasCycleModel(species::Vector{String},
                           a::Matrix{Real},
                           τ::Matrix{Real},
                           r0::Vector{Real},
                           ru::Vector{Real},
                           rT::Vector{Real},
                           ra::Vector{Real},
                           C₀::Vector{Real},
                           EtoC::Vector{Real},
                           f::Matrix{Real})
            g₀, g₁ = compute_g(a, τ)
            new(species, a, τ, r0, ru, rT, ra, C₀, EtoC, f, g₀, g₁)
    end

    function GasCycleModel(path::String, species::Vector{String})
        params = load_gascycle_params(path, species)
        g₀, g₁ = compute_g(params.a, params.τ)
        println(size(g₀))
        new(species, params.a, params.τ, params.r0, params.ru, params.rT, params.ra, params.C₀, params.EtoC, params.f, g₀, g₁)
    end
end


function compute_g(a, τ)
    δ = 100 ./ τ
    e⁻ᵟ = exp(-δ)
    g₁ = sum(a .* τ .* (1 .- (1 .+ δ) .* e⁻ᵟ), dims=2)
    g₀ = exp.(-sum(a .* τ .* (1 .- e⁻ᵟ), dims=2) ./ g₁)
    return vec(g₀), vec(g₁)
end

# gascycle = GasCycleModel("src/defaults/gas_parameters.csv", ["CO2", "CH4", "SO2", "BC"])
