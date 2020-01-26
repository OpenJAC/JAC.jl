#
println("Ea) Two-step cascade computations and simulations following the 1s-3p photo-excitation of Si^-: AverageSCA model.")
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
#
if  false    ## Stepwise decay cascade
setDefaults("print summary: open", "zzz-Cascade-computation.sum")

name = "Cascade after Si- 1s-3p excitation"
grid = Radial.Grid(false)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                           scheme=Cascade.StepwiseDecayScheme([Auger, Radiative], 2, 0, Shell[], Shell[]),
                           initialConfigs=[Configuration("1s^1 2s^2 2p^6 3s^2 3p^4")] )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

else
setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

using JLD
## data = [JLD.load("zzz-cascade-decay-computations-2020-01-25T11.jld"), JLD.load("zzz-cascade-decay-computations-2020-01-25T11.jld")]
data = [JLD.load("zzz-cascade-decay-computations-2020-01-25T21.jld")]
name = "Simulation after Si- 1s-3p excitation"

wc   = Cascade.Simulation(Cascade.Simulation(), name=name, properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution()], ## , Cascade.FinalLevelDistribution()], 
                          settings=Cascade.SimulationSettings(0., 0., 0., 0., [(1, 2.0), (2, 1.0), (3, 0.5)]), computationData=data )
println(wc)
wd = perform(wc; output=true)
setDefaults("print summary: close", "")

end
