module EnergyBalanceModel

export EBM, FtoT, run, computeA, compute_bd, samplevariability

include("EBM.jl")

# C = @SVector [1.0, 2.0, 3.0]  # Static vector of size 3
# κ = @SVector [0.5, 1.0, 1.5]  # Same size static vector
# ebm = EBM(3, C, κ, 0.7, 1.0, 0.1, 0.05, 0.9, 123, 0.01, 1000)
# RF = collect(1:ebm.Nₜ) ./ 300.
# Tₛ, 𝓕 = run(ebm, RF)
# plot(collect(1:length(Tₛ)), 𝓕)

end

