#
println("Ea) Three-step cascade computation and simulation after 1s-3p photo-excitation of Si^-: AverageSCA model.")
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
#
if  false    ## Stepwise decay cascade
setDefaults("print summary: open", "zzz-Cascade-computation.sum")

name = "Si- (1s^-1) cascade"
grid = Radial.Grid(false)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                           scheme=Cascade.StepwiseDecayScheme([Auger, Radiative], 1, 0, Shell[], Shell[]),
                           initialConfigs=[Configuration("1s^1 2s^2 2p^6 3s^2 3p^4")] )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

elseif  true    ## Photoionization cascade
setDefaults("print summary: open", "zzz-Cascade-computation.sum")

name = "neutral Si- photoionization cascade"
grid = Radial.Grid(false)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                           scheme=Cascade.PhotonIonizationScheme([Photo], 1, [30.0, 80.0]),
                           initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p^3"), Configuration("1s^2 2s^2 2p^6 3s^1 3p^4")] )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

else
setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

using JLD
data = JLD.load("zzz-Cascade-computation-2020-01-20T20.jld")["cascade decay data:"]
name = "Si- (1s^-1) simulation"

wc   = Cascade.Simulation(Cascade.Simulation(), name=name, properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution(), Cascade.FinalLevelDistribution()], 
                          settings=Cascade.SimulationSettings(0., 0., 0., 0., [(1, 1.0), (2, 1.0), (3, 0.5)]), cascadeData=data )
println(wc)
wd = perform(wc; output=true)
setDefaults("print summary: close", "")

end
