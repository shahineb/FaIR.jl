using FaIR
using Statistics
using GLMakie


# ## Path to files to load
#
# See files in src/defaults/ for the format these files must have.
# The species config file follows the format from https://github.com/OMS-NetZero/FAIR/blob/master/src/fair/defaults/data/ar6/species_configs_properties.csv
# The energy balance parameters file follows the format from https://github.com/OMS-NetZero/FAIR/blob/master/tests/test_data/4xCO2_cummins_ebm3.csv
#
# At the moment this code support emissions as an input type. The input data 
# for FaIR is yearly global emission data for different species.

emission_csv = "src/defaults/ssp245-emissions.csv"             # Input emission data
species_csv = "src/defaults/species_configs_properties.csv"    # Species configs
ebm_csv = "src/defaults/MPI-ESM1-2-LR_4xCO2_cummins_ebm3.csv"  # EBM params



# ## FaIR components instantiation
#
# ### Load emissions
# 
# Emissions must include at least CO₂, CH₄ and N₂O. If data is missing for one, you need to
# set it to zero in the emission file. (TODO : take care of this internally)
#
# The order in which the species are specified specifies the rows on which they act
# in the emission, concentration and forcing arrays.
species = ["CO2", "CH4", "N2O", "SO2", "BC"]
E = Emissions(emission_csv, species)




# ### Initialise gas cycle model
gas_model = ReservoirModel(species_csv, species)




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
Nseeds = 50                             # number of seeds to sample internal variability
Δt = 1.                                 # yearly time step
Nₜ = length(E.year)                      # number of years 
ebm = BoxModel(ebm_csv, Δt, Nₜ)
eᴬ, bd, _ = ebm_dynamics(ebm)




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
# Here we add an extra dimension for each seed we will sample internal variability from.
T = zeros(Real, Nseeds, ebm.Nₜ, ebm.Nbox + 1)
active_dims_ghg = [E.index.CO2, E.index.CH4, E.index.N2O] # TODO : not needing this anymore






# ### Run FaIR once for every seed
#
n_species = length(species)
for ω in collect(1:Nseeds)
    # Sample from internal variability
    wd = samplevariability(ebm, ω)

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
for ω in 1:Nseeds
    lines!(ax, E.year, T[ω, :, 2], color=(:blue, 0.05))
end
lines!(ax, E.year, vec(mean(T[:, :, 2], dims=1)), color=(:blue, 0.8), label="Mean")
axislegend(ax, position=:lt)
display(fig)
save("examples/demo-internal-variability.png", fig)



GLMakie.closeall()