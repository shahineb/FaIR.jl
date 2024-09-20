using FaIR
using Statistics
using GLMakie


# ## Path to files to load
#
# See files in src/defaults/ for the format these files must have.
#
# At the moment this code support emissions as an input type. The input data 
# for FaIR is yearly global emission data for different species.

emission_csv = "src/defaults/ssp245-emissions.csv"                                                 # Input emission data
constrained_ensemble_csv = "src/defaults/calibrated_constrained_parameters_calibration1.4.1.csv"   # Set of all parameters from constrained ensemble 
species_csv = "src/defaults/species_configs_properties.csv"                                        # Default species configs (for the parameters not sampled in the ensemble)
ebm_csv = "src/defaults/MPI-ESM1-2-LR_4xCO2_cummins_ebm3.csv"                                      # Default EBM params (for the parameters not sampled in the ensemble) TODO : this shouldn't be ESM specific


# ### Load emissions
# 
# Emissions must include at least CO₂, CH₄ and N₂O. If data is missing for one, you need to
# set it to zero in the emission file. (TODO : take care of this internally)
#
# The order in which the species are specified specifies the rows on which they act
# in the emission, concentration and forcing arrays.
species = ["CO2", "CH4", "N2O", "SO2", "BC"]
E = Emissions(emission_csv, species)




# ### Load constrained ensemble dataframe
using CSV, DataFrames
df = CSV.read(constrained_ensemble_csv, DataFrame)

γ = df.gamma_autocorrelation
C = Matrix(df[:, ["ocean_heat_capacity[0]", "ocean_heat_capacity[1]", "ocean_heat_capacity[2]"]])
κ = Matrix(df[:, ["ocean_heat_transfer[0]", "ocean_heat_transfer[1]", "ocean_heat_transfer[2]"]])
ε = df.deep_ocean_efficacy
ση = df.sigma_eta
σξ = df.sigma_xi
F₄ₓ = df.forcing_4co2

r0CO2 = df[:, "iirf_0[CO2]"]
ruCO2 = df[:, "iirf_uptake[CO2]"]
rTCO2 = df[:, "iirf_temperature[CO2]"]
raCO2 = df[:, "iirf_airborne[CO2]"]

radiative_efficiency_ari = Matrix(df[:, ["erfari_radiative_efficiency[CH4]", "erfari_radiative_efficiency[N2O]", "erfari_radiative_efficiency[Sulfur]", "erfari_radiative_efficiency[BC]"]])
s = Matrix(df[:, ["aci_shape[Sulfur]", "aci_shape[BC]"]])
β = df.aci_scale
meinshausen_scaling = Matrix(df[:, ["forcing_scale[CO2]", "forcing_scale[CH4]", "forcing_scale[N2O]"]])
forcing_scale_SO2 = ones(size(meinshausen_scaling)[1])
forcing_scale_BC = ones(size(meinshausen_scaling)[1])
ari_scaling = hcat(meinshausen_scaling[:, 2:end], forcing_scale_SO2, forcing_scale_BC)
aci_scaling = hcat(forcing_scale_SO2, forcing_scale_BC)


C₀_CO₂ = df[:, "baseline_concentration[CO2]"]

N_members = size(df)[1]




# ### Initialise default gas cycle model
gm = ReservoirModel(species_csv, species)




# ### Initialise main greenhouse gas forcing model
# 
# This is the forcing model for CO₂, CH₄, N₂O. It must be included in any run.
#
# We need to provide it with row index on which CO₂, CH₄ and  N₂O are active
# in the emission, concentration and forcing arrays (in this particular order).
# In the example below, CO₂ corresponds to the 1ˢᵗ row, CH₄ on the 2ⁿᵈ and N₂O the 3ʳᵈ
active_dims_CO₂_CH₄_N₂O = [1, 2, 3]
CO₂_CH₄_N₂O_forcing = Meinshausen2020(species_csv, active_dims_CO₂_CH₄_N₂O)





# ### Initialise other optional forcing models
#
# These are additional models that allow to account for additional sources of
# forcing such as minor greenhouse gas, aerosols, contrails, land use etc.
#
# Each forcing model needs to be instantiated with a list of species which effect
# we want to account for, and the dimensions (row index) on which they are active
# in the emission, concentration and forcing arrays

# #### Forcing model for aerosol direct effect
species_ari = ["CH4", "N2O", "SO2", "BC"]
active_dims_ari = [2, 3, 4, 5]
ari_forcing = ARIForcing(species_csv, species_ari, active_dims_ari)

# #### Forcing model for aerosol-cloud interactions
species_aci = ["SO2", "BC"]
active_dims_aci = [4, 5]
aci_forcing = ACIForcing(species_csv, species_aci, active_dims_aci)





# ### Initialise energy balance model
#
# The energy balance model is here setup as a 3-box model and includes
# a stochastic internal variability components following Cummins et al. 2020.
#
# We define variables for the energy balance model but only initialise it
# within the loop below since we need one different model per set of parameters
Nbox = 3                                # number of temperature boxes in the box model
use_internal_variability = false        # whether or not sample from internal variability
Δt = 1.                                 # yearly time step
Nₜ = length(E.year)                      # number of years 
ebm = BoxModel(ebm_csv, Δt, Nₜ)





# ### Allocate arrays for the run
#
function allocate_gas_cycle(n_species)
    n_pool = 4
    pool_partition = zeros(Real, n_species, n_pool)
    airborneₜ = zeros(Real, n_species)
    E₀ = zeros(length(species))
    Cₜ₋₁ = zeros(Real, length(species))
    C₀ = gas_model.C₀
    return pool_partition, airborneₜ, E₀, Cₜ₋₁, C₀
end


# For the energy balance model, we intialise an array of each box temperature at
# each time step, with an additional box dedicated to the TOA forcing fluctuations
# following the procedure in Cummins et al. 2020. 
#
# Here we add an extra dimension for ensemble member we're running.
T = zeros(Real, N_members, Nₜ, Nbox + 1)
active_dims_ghg = [E.index.CO2, E.index.CH4, E.index.N2O] # TODO : not needing this anymore





# ### Run FaIR once for every seed
#
n_species = length(species)
for ω in collect(1:N_members)
    # Initialise gas cycle model
    r0 = vcat(r0CO2[ω], gm.r0[2:end])
    ra = vcat(raCO2[ω], gm.ra[2:end])
    ru = vcat(ruCO2[ω], gm.ru[2:end])
    rT = vcat(rTCO2[ω], gm.rT[2:end])
    C₀ = vcat(C₀_CO₂[ω], gm.C₀[2:end])
    gm = ReservoirModel(species, gm.a, gm.τ, r0, ru, rT, ra, C₀, gm.molecular_weight, gm.χ_sensitivity_τᶜᴴ⁴, gm.T_sensitivity_τᶜᴴ⁴, gm.idx_E, gm.idx_C)

    # Initialise forcing models
    C₀ = vcat(C₀_CO₂[ω], CO₂_CH₄_N₂O_forcing.C₀[2:end])
    CO₂_CH₄_N₂O_forcing = Meinshausen2020(C₀, meinshausen_scaling[ω, :], CO₂_CH₄_N₂O_forcing.radiative_efficiency, active_dims_CO₂_CH₄_N₂O)
    ari_forcing = ARIForcing(species_ari, ari_scaling[ω, :], ari_forcing.radiative_efficiency, ari_forcing.E₀, ari_forcing.C₀, ari_forcing.idx_E, ari_forcing.idx_C, active_dims_ari)
    aci_forcing = ACIForcing(species_aci, aci_scaling[ω, :], s[ω, :], β[ω] .* ones(2), aci_forcing.E₀, aci_forcing.C₀, aci_forcing.idx_E, aci_forcing.idx_C, active_dims_aci)


    # Initialise a box-model with the seed
    ebm = BoxModel(C[ω, :], κ[ω, :], ε[ω], F₄ₓ[ω], ση[ω], σξ[ω], γ[ω], Δt, Nₜ)
    eᴬ, bd, wd = ebm_dynamics(ebm, use_internal_variability)

    # Allocate arrays
    pool_partition, airborneₜ, E₀, Cₜ₋₁ , C₀ = allocate_gas_cycle(n_species)
    for t in 2:Nₜ
        # Compute feedback gas lifetime scaling factor
        αₜ = α(gas_model, airborneₜ, E.cumulative[:, t - 1], T[ω, t - 1, 2])
        αₜ[E.index.CH4] = αᶜᴴ⁴(gas_model, E.values[:, t - 1], Cₜ₋₁ , E₀, C₀, T[ω, t - 1, 2])

        # Compute concentrations and reservoirs partition
        Cₜ, pool_partition = EtoC(gas_model, E.values[:, t - 1], pool_partition, αₜ, E₀, Δt)

        # Compute forcings from the 3 major greenhouse gases
        FCO₂_CH₄_N₂O = computeF(CO₂_CH₄_N₂O_forcing, Cₜ)

        # Compute aerosol direct effect forcing
        Fari = computeF(ari_forcing, E.values[:, t], Cₜ)

        # Compute aerosol-cloud interaction forcing
        Faci = computeF(aci_forcing, E.values[:, t], Cₜ)

        # Compute total effective radiative forcing
        Fₜ = sum(map(sum, [FCO₂_CH₄_N₂O, Fari, Faci]))
        
        # Compute boxes temperatures responses to forcing
        T[ω, t, :] = FtoT(T[ω, t - 1, :], eᴬ, bd, wd[t - 1, :], Fₜ)
    end
end





# Temperature plot
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Year", ylabel = "ΔT (K)")
for ω in 1:N_members
    lines!(ax, E.year, T[ω, :, 2], color=(:blue, 0.05))
end
lines!(ax, E.year, vec(mean(T[:, :, 2], dims=1)), color=(:purple, 0.8), label="Mean")
axislegend(ax, position=:lt)
display(fig)
save("examples/constrained-ensemble-projections.png", fig)



GLMakie.closeall()