# fair.jl


[x] Clean up ebm (make the forward methods compact like for the other modules)
[x] Change dirty constructors to clean constructing functions
[x] Implement more sophisticated forcing function that allows to load the actual file

[ ] Implement linear forcing and aerosol forcing modules that load all the parameters needed from csv (scaling, efficiency, concentration_idx, emission_idx) so that I can then only pass that to the calculation function
[ ] Figure out whether implementation of functions makes sense in terms of vectors dimensions (related to above task)
[ ] Implement everything dirtily with SO2 being handled on the side as linear + ACI


[ ] Reconvene to find a sensible way to refactor forcing module. Ideally, I would want to have a single forcing model initialisation which takes as arguments the other forcing models. The problem is that the split isn't clear because linear also applies to some ghg, not all ghg need to be used at the same time.
    - Category 1 : necessary agents (CO2, CH4, N2O) -> need to be used so might as well set emissions to zero if they're not there. Can manipulate their indices explicitely since they will always be there. Maybe I can define a function that maps any concentration vector to a 3-vector with zeros in case some are missing

    - Category 2 : optional agents (SO2, BC, minor GHG) -> do not need to be used, so need to find a way to distinguish them in the emission/concentration vector. I think the Input instance should have a split between necessary and unnecessary forcings.

    - Category 3 : non-species specific forcing (aci, ari, volcanic) -> a bit odd, stands on the side, it's not tied with one specific species, might still using emission/concentration data from multiple species or none if prescribed. Yet, in the end it's an additional entry in the radiative forcing vector.

Ooooh but maybe I can define an addition between forcings? Bo, that would be so sweet (but also unnecessary I know)

Maybe I can have a general forcing struct which takes all of these into it, but doesn't necessarily need all of them (except Category 1)

[ ] Reimplement the forcing module in terms of categories of forcing, along with the Input module which should have indices for categories of forcing

filter(kv -> kv ∈ [:CO2, :CH4], keys(E.index)) -> maybe a way to filter based on species but looks cumbersome


[ ] Methane lifetime update -> might mean I need to think of a way to distinguish species in alpha
[ ] Allow to use prescribed forcing for Volcanic activity
[ ] Make it possible to load multiple configs (not too hard, only applies to ebm)
[ ] Start documenting code
[x] Harmonise modules, struct, function naming
[ ] Add an option to sample multiple trajectories from internal variability
[ ] Sampling from parameters?
[ ] Tuning against parameters?