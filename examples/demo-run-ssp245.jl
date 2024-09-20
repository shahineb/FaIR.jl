using FaIR
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
#
# We precompute the temperature feedback matrix eᴬ, forcing feedback vector bd
# and the sampled internal variability trajectory wd TODO : find a way to make wd optional

use_internal_variability = false        # whether or not sample from internal variability
Δt = 1.                                 # yearly time step
Nₜ = length(E.year)                      # number of years 
ebm = BoxModel(ebm_csv, Δt, Nₜ)
eᴬ, bd, wd = ebm_dynamics(ebm, use_internal_variability)






# ### Allocate arrays for the run
#
# For the gas cycle model, we initialise four arrays :
#   - `pool_partition` : running per specie burden in each atmospheric lifetime pool
#   - `airborneₜ` : per specie total cumulative emission remaining in the atmosphere
#   - `E₀` : baseline emissions for each specie
#   - `C` : per specie total atmospheric concentration at each time step
n_pool = 4
n_species = length(species)
pool_partition = zeros(Real, n_species, n_pool)
airborneₜ = zeros(Real, n_species)
E₀ = zeros(length(species))
C = zeros(Real, size(E.values))
C[:, 1] = gas_model.C₀


# For the forcing we initialise an array of forcing per specie at each time step
F = zeros(Real, size(E.values))

# For the energy balance model, we intialise an array of each box temperature at
# each time step, with an additional box dedicated to the TOA forcing fluctuations
# following the procedure in Cummins et al. 2020.
T = zeros(Real, ebm.Nₜ, ebm.Nbox + 1)
active_dims_ghg = [E.index.CO2, E.index.CH4, E.index.N2O] # TODO : not needing this anymore






# ### Run FaIR
#
for t in 2:Nₜ
    # Compute feedback gas lifetime scaling factor
    αₜ = α(gas_model, airborneₜ, E.cumulative[:, t - 1], T[t - 1, 2])
    αₜ[E.index.CH4] = αᶜᴴ⁴(gas_model, E.values[:, t - 1], C[:, t - 1], E₀, gas_model.C₀, T[t - 1, 2])

    # Compute concentrations and reservoirs partition
    C[:, t], pool_partition = EtoC(gas_model, E.values[:, t - 1], pool_partition, αₜ, E₀, Δt)

    # Compute forcings from the 3 major greenhouse gases
    FCO₂_CH₄_N₂O = computeF(CO₂_CH₄_N₂O_forcing, C[:, t])
    F[active_dims_ghg, t] .+= FCO₂_CH₄_N₂O

    # Compute aerosol direct effect forcing
    Fari = computeF(ari_forcing, E.values[:, t], C[:, t])
    F[ari_forcing.active_dims, t] .+= Fari

    # Compute aerosol-cloud interaction forcing
    Faci = computeF(aci_forcing, E.values[:, t], C[:, t])
    F[aci_forcing.active_dims, t] .+= Faci
    
    # Compute boxes temperatures responses to forcing
    T[t, :] = FtoT(T[t - 1, :], eᴬ, bd, wd[t - 1, :], sum(F[:, t]))
end





# Emissions plot
fig = Figure()
units = ["GtC/yr", "Mt CH₄/yr", "Mt N₂O/yr", "Mt SO₂/yr", "Mt BC/yr"]
for i in 1:n_species
    ax = Axis(fig[1, i]; xlabel = "Year", ylabel = units[i], title = species[i])
    lines!(ax, E.year, E.values[i, :], label=species[i])
    axislegend(ax, position=:rb)
end
display(fig)


# Concentrations plot
fig = Figure()
units = ["ppm", "ppb", "ppb"]
for i in 1:n_species-2
    ax = Axis(fig[1, i]; xlabel = "Year", ylabel = units[i], title = species[i])
    lines!(ax, E.year, C[i, :], label=species[i])
    axislegend(ax, position=:rb)
end
display(fig)


# Forcing plot
fig = Figure()
for i in 1:n_species
    ax = Axis(fig[1, i]; xlabel = "Year", ylabel = "Wm⁻²", title = species[i])
    lines!(ax, E.year, F[i, :], label=species[i])
    axislegend(ax, position=:lt)
end
ax = Axis(fig[1, n_species + 1]; xlabel = "Year", ylabel = "Wm⁻²", title = "Total")
lines!(ax, E.year, vec(sum(F, dims=1)), label="Total")
axislegend(ax, position=:lt)
display(fig)


# Temperature plot
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Year", ylabel = "ΔT (K)", title = "GMST")
lines!(ax, E.year, T[:, 2])
display(fig)



GLMakie.closeall()