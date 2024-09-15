using CSV, DataFrames

function load_ebm_params(path::String)
    df = CSV.read(path, DataFrame)
    row = filter(row -> row.model == "CAMS-CSM1-0", df)
    γ = row.gamma[1]
    C = vcat(row.C1, row.C2, row.C3)
    κ = vcat(row.kappa1, row.kappa2, row.kappa3)
    ε = row.epsilon[1]
    ση = row.sigma_eta[1]
    σξ = row.sigma_xi[1]
    F₄ₓ = row.F_4xCO2[1]
    output = (γ=γ, C=C, κ=κ, ε=ε, ση=ση, σξ=σξ, F₄ₓ=F₄ₓ)
    return output
end