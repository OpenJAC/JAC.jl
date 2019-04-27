#
println("Dcp) Test of the Photoionization and PlasmaShift module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoIonization-Plasma.sum")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(36.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, NoPoints = 600),
                        initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                        process = JAC.PhotoInPlasma, 
                        processSettings=PlasmaShift.PhotoSettings(PlasmaShift.DebyeHueckel, 0.25, 0., 0, [E1], [JAC.UseCoulomb], 
                                                                  [3000., 3200.], true, false, Tuple{Int64,Int64}[]) )

wb = perform(wa)
setDefaults("print summary: close", "")

