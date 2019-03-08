#
println("Dc) Test of the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
JAC.define("print summary: open", "zzz-PhotoIonization.sum")
JAC.define("method: continuum, Galerkin")
JAC.define("method: normalization, pure sine")   ## JAC.define("method: normalization, pure Coulomb")    JAC.define("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(36.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, NoPoints = 600),
                        initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                        process = JAC.Photo, 
                        processSettings=PhotoIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [3000., 3200.], false, true, true, true, 
                                        false, Tuple{Int64,Int64}[]) )

wb = perform(wa)
JAC.define("print summary: close", "")

