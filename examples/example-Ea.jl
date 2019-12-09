#
println("Ea) Three-step cascade after 1s-3p photo-excitation of Si^+: Configuration model.")
#
if  false
setDefaults("print summary: open", "zzz-Cascade-computation.sum")
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

name = "Si- (1s^-1) cascade"
grid = Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600)
wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                           processes=[Auger, Radiative], initialConfigs=[Configuration("1s^1 2s^2 2p^6 3s^2 3p^4")], initialLevels=[(1, 1.)],
                           maxElectronLoss=2 )
println(wa)
wb = perform(wa; output=true)
setDefaults("print summary: close", "")

else
setDefaults("print summary: open", "zzz-Cascade-simulatiaon.sum")

using JLD
data = JLD.load("zzz-Cascade-2019-12-03T20.jld")["cascade data:"]
name = "Si- (1s^-1) simulation"

wc   = Cascade.Simulation(Cascade.Simulation(), name=name )
## wc   = Cascade.Simulation(Cascade.Simulation(); name=name, cascadeData=data )
println(wc)
## wd = perform(wc; output=true)
setDefaults("print summary: close", "")

end
