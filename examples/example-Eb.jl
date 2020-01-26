#
println("Eb) Three-step cascade computation and simulation for the photo-ionization of Si^- and its subsequent decay: AverageSCA model.")
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
#
if  false    ## Photoionization cascade
setDefaults("print summary: open", "zzz-Cascade-computation-photoionization.sum")

name = "Photoionization of Si- "
grid = Radial.Grid(false)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                           scheme=Cascade.PhotonIonizationScheme([Photo], 1, [30.0, 80.0]),
                           initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")] )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

elseif  false    ## Stepwise decay cascade
setDefaults("print summary: open", "zzz-Cascade-computation-following-decay.sum")

using JLD
data = JLD.load("zzz-cascade-ionizing-computations-2020-01-26T10.jld")
iniMultiplets = data["generated multiplets:"]
## lineData      = data["photo-ionizing line data:"]

name = "Si- (1s^-1) decay cascade"
grid = Radial.Grid(false)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                           scheme=Cascade.StepwiseDecayScheme([Auger, Radiative], 1, 0, Shell[], Shell[]),
                           initialMultiplets=iniMultiplets )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

else
setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

using JLD
data = [JLD.load("zzz-cascade-ionizing-computations-2020-01-26T10.jld"), JLD.load("zzz-cascade-decay-computations-2020-01-26T10.jld")]
name = "Simulation after Si- 1s and 2s ionization and subsequent decay"

wc   = Cascade.Simulation(Cascade.Simulation(), name=name, properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution()], ## , Cascade.FinalLevelDistribution()], 
                          settings=Cascade.SimulationSettings(0., 0., 0., 0., [(1, 2.0), (2, 1.0), (3, 0.5)]), computationData=data )
println(wc)
wd = perform(wc; output=true)
setDefaults("print summary: close", "")

end
