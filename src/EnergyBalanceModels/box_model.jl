using LinearAlgebra, Random, Distributions

include("load_csv.jl")

struct BoxModel
    Nbox::Int                  # Number of boxes
    C::AbstractVector{<:Real}  # Heat capacity of each box W m⁻² yr K⁻¹
    κ::AbstractVector{<:Real}  # Heat exchange coefficients W m⁻² K⁻¹
    ε::Real                    # Deep ocean heat uptake efficacy
    F₄ₓ::Real                  # Effective RF for 4xCO₂
    ση::Real                   # Stddev of the white noise in the radiative forcing
    σξ::Real                   # Stddev of the white noise in the temperature
    γ::Real                    # Autocorrelation parameter of stochastic forcing
    seed::Int                  # Random seed for stochastic modelling
    Δt::Real                   # Time step yr
    Nₜ::Int                    # Number of time steps

    function BoxModel(C::AbstractVector{<:Real},
                      κ::AbstractVector{<:Real},
                      ε::Real,
                      F₄ₓ::Real,
                      ση::Real,
                      σξ::Real,
                      γ::Real,
                      seed::Int,
                      Δt::Real,
                      Nₜ::Int)
        Nbox = length(C)
        new(Nbox, C ./ Δt, κ, ε, F₄ₓ, ση, σξ, γ, seed, Δt, Nₜ)
    end

    function BoxModel(path::String,
                      seed::Int,
                      Δt::Real,
                      Nₜ::Int)
        params = load_box_model_params(path)
        Nbox = length(params.C)
        new(Nbox, params.C ./ Δt, params.κ, params.ε, params.F₄ₓ, params.ση, params.σξ, params.γ, seed, Δt, Nₜ)
    end
end


function computeA(ebm::BoxModel)
    # Instantiate matrix
    A = zeros(ebm.Nbox + 1, ebm.Nbox + 1)
    ε⁺ = ones(ebm.Nbox)
    ε⁺[ebm.Nbox - 1] = ebm.ε

    # Include stochastic terms
    A₁₁ = -ebm.γ
    A₂₁ = 1 / ebm.C[1]
    A[1:2, 1] = [A₁₁, A₂₁]

    # First box
    A₂₂ = -(ebm.κ[1] + ε⁺[1] * ebm.κ[2]) / ebm.C[1]
    A₂₃ = ε⁺[1] * ebm.κ[2] / ebm.C[1]
    A[2, 2:3] = [A₂₂, A₂₃]

    # Last box
    Aₙₙ₋₁ = ebm.κ[end] / ebm.C[end]
    Aₙₙ = -Aₙₙ₋₁
    A[end, end-1:end] = [Aₙₙ₋₁, Aₙₙ]

    # Intermediate boxes if Nbox > 2
    for i in 2:ebm.Nbox-1
        Aᵢᵢ₋₁ = ebm.κ[i] / ebm.C[i]
        Aᵢᵢ = -(ebm.κ[i] + ε⁺[i] * ebm.κ[i + 1]) / ebm.C[i]
        Aᵢᵢ₊ᵢ = ε⁺[i] * ebm.κ[i + 1] / ebm.C[i]
        A[i + 1, i:i+2] = [Aᵢᵢ₋₁, Aᵢᵢ, Aᵢᵢ₊ᵢ]
    end
    return A
end


function compute_bd(ebm::BoxModel)
    # Compute exp(A)
    A = computeA(ebm)
    eᴬ = exp(A)

    # Compute vector forcing update bd = A⁻¹(eᴬ - I)b
    b = zeros(ebm.Nbox + 1)
    b[1] = ebm.γ
    bd = A \ ((eᴬ - I) * b)
    return bd
end


function samplevariability(ebm::BoxModel)
    A = computeA(ebm)
    Nₐ = size(A, 1)
    Q = zeros(Nₐ, Nₐ)
    Q[1, 1] = ebm.ση^2
    Q[2, 2] = (ebm.σξ / ebm.C[1])^2
    H = [-A Q ; zeros(Nₐ, Nₐ) A']
    eᴴ = exp(H)
    Qd = eᴴ[Nₐ + 1:end, Nₐ + 1:end]' * eᴴ[1:Nₐ, Nₐ + 1:end]
    Qd = 0.5 .* (Qd .+ Qd')
    Random.seed!(ebm.seed)
    𝒩 = MvNormal(zeros(Nₐ), Qd)
    wd = rand(𝒩, ebm.Nₜ)'
    return wd
end


function impulseresponse(ebm::BoxModel)
    # Compute eigendecomposition of A
    A = computeA(ebm)
    Λ, Φ = eigen(A[2:end, 2:end])

    # Compute timescales d and equilibrium response q (notations from Millar et al. 2017)
    d = -ebm.Δt ./ real(Λ)
    q = d .* (Φ[1, :] .* inv(Φ)[:, 1]) ./ (ebm.C[1] * ebm.Δt)
    return d, q
end


function emergentparameters(ebm::BoxModel, ratio₂ₓ₄ₓ=0.5)
    d, q = impulseresponse(ebm)
    τ₂ₓ = log(2) / log(1.01)
    ecs = ebm.F₄ₓ * ratio₂ₓ₄ₓ * sum(q)
    tcr = ebm.F₄ₓ * ratio₂ₓ₄ₓ * sum(q .* (1 .- (d ./ τ₂ₓ) .* (1 .- exp.(-τ₂ₓ ./ d))))
    return ecs, tcr
end


function ebm_dynamics(ebm::BoxModel)
    # Compute exp(A)
    A = computeA(ebm)
    eᴬ = exp(A)

    # Compute vector forcing update bd = A⁻¹(eᴬ - I)b
    b = zeros(ebm.Nbox + 1)
    b[1] = ebm.γ
    bd = A \ ((eᴬ - I) * b)

    # Sample internal variability updates
    wd = samplevariability(ebm)
    return eᴬ, bd, wd
end


function FtoT(Tₜ₋₁, eᴬ, bd, wdₜ₋₁, RFₜ₋₁)
    Tₜ = eᴬ * Tₜ₋₁ .+ bd .* RFₜ₋₁ .+ wdₜ₋₁
    return Tₜ
end


function run(ebm::BoxModel, RF::AbstractVector{<:Real})
    # Compute exp(A)
    A = computeA(ebm)
    eᴬ = exp(A)

    # Compute vector forcing update bd = A⁻¹(eᴬ - I)b
    b = zeros(ebm.Nbox + 1)
    b[1] = ebm.γ
    bd = A \ ((eᴬ - I) * b)

    # Sample internal variability updates
    wd = samplevariability(ebm)

    # Compute boxes temperatures
    T = zeros(ebm.Nₜ, ebm.Nbox + 1)
    for i in 2:ebm.Nₜ
        T[i, :] = eᴬ * T[i - 1, :] .+ bd .* RF[i - 1] .+ wd[i - 1, :]
    end

    # Return surface temperature and stochastic forcing time serie
    Tₛ = T[:, 2]
    𝓕 = T[:, 1]
    return Tₛ, 𝓕
end