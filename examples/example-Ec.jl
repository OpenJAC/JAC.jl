
#
println("Ec) Three-step cascade computations and simulations for the decay of the neon 1s^-1 3p hole states: AverageSCA model.")
using JLD2
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
#
if false 

elseif  false    ## Stepwise decay cascade
setDefaults("print summary: open", "zzz-Cascade-computation-following-decay.sum")

name = "Cascade after neon 1s --> 3p excitation"
grid = Radial.Grid(false)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                           scheme=Cascade.StepwiseDecayScheme([Auger, Radiative], 5, 0, Shell[], Shell[]),
                           initialConfigs=[Configuration("1s^1 2s^2 2p^6 3p")] )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

else
setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

data = [wb]
name = "Simulation of the neon 1s^-1 3p decay"

wc   = Cascade.Simulation(Cascade.Simulation(), name=name, properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution(), Cascade.FinalLevelDistribution()], 
                          settings=Cascade.SimulationSettings(0., 0., 0., 0., 0., [(1, 2.0), (2, 1.0), (3, 0.5)]), computationData=data )
println(wc)
wd = perform(wc; output=true)
setDefaults("print summary: close", "")

end
