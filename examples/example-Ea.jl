#
println("Ea) Three-step cascade computation and simulation after 1s-3p photo-excitation of Si^-: AverageSCA model.")
#
if  false
setDefaults("print summary: open", "zzz-Cascade-computation.sum")
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

name = "Si- (1s^-1) cascade"
grid = Radial.Grid(false)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                           processes=[Auger, Radiative], initialConfigs=[Configuration("1s^1 2s^2 2p^6 3s^2 3p^4")],
                           maxElectronLoss=4 )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

else
setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

using JLD
data = JLD.load("zzz-Cascade-computation-2020-01-17T21.jld")["cascade data:"]
name = "Si- (1s^-1) simulation"

wc   = Cascade.Simulation(Cascade.Simulation(), name=name, properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution(), Cascade.FinalLevelDistribution()], 
                          settings=Cascade.SimulationSettings(0., 0., 0., 0., [(1, 1.0), (2, 1.0), (3, 0.5), (60, 0.25)]), cascadeData=data )
println(wc)
wd = perform(wc; output=true)
setDefaults("print summary: close", "")

end
