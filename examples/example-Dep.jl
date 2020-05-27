#
println("Dep) Test of the AutoIonization and PlasmaShift modules with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-AutoIonization-Plasma.sum")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

grid = JAC.Radial.Grid(false)
grid = JAC.Radial.Grid(grid, rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, NoPoints = 900)
wa   = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                          initialConfigs=[Configuration("1s 2s^2 2p^6")],
                          finalConfigs  =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], 
                          process = AugerInPlasma(),
                          processSettings = PlasmaShift.AugerSettings(PlasmaShift.DebyeHueckel(), 0.25, 0., 0, true, true, Tuple{Int64,Int64}[(1,0)]) )

wb = perform(wa)
setDefaults("print summary: close", "")


