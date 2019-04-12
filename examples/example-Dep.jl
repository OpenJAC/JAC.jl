#
println("Dep) Test of the AutoIonization and PlasmaShift modules with ASF from an internally generated initial- and final-state multiplet.")
#
JAC.define("print summary: open", "zzz-AutoIonization-Plasma.sum")
JAC.define("method: continuum, Galerkin")
JAC.define("method: normalization, pure sine")   ## JAC.define("method: normalization, pure Coulomb")    JAC.define("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(10.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, NoPoints = 600), 
                        initialConfigs=[Configuration("1s 2s^2 2p^6")],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], 
                        process = JAC.AugerInPlasma,
                        processSettings = PlasmaShift.AutoIonizationSettings(PlasmaShift.DebyeHueckel, 0.25, 0., 0, true, true, Tuple{Int64,Int64}[(1,0)]) )

wb = perform(wa)
JAC.define("print summary: close", "")


