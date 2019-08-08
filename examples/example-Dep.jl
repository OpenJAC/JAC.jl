#
println("Dep) Test of the AutoIonization and PlasmaShift modules with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-AutoIonization-Plasma.sum")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(10.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, NoPoints = 900), 
                        initialConfigs=[Configuration("1s 2s^2 2p^6")],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], 
                        process = JAC.AugerInPlasma,
                        processSettings = PlasmaShift.AugerSettings(PlasmaShift.DebyeHueckel(), 0.25, 0., 0, true, true, Tuple{Int64,Int64}[(1,0)]) )

wb = perform(wa)
setDefaults("print summary: close", "")


