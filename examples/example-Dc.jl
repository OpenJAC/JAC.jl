#
println("Dc) Test of the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
JAC.define("print summary: open", "zzz-PhotoIonization.sum")
JAC.define("method: continuum, asymptotic Coulomb")  ## JAC.define("method: continuum, Galerkin")
JAC.define("method: normalization, pure Coulomb")   ## JAC.define("method: normalization, pure Coulomb")    JAC.define("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(20.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 800),
                        initialConfigs=[Configuration("[Ne] 3s^2 3p^6")],
                        finalConfigs  =[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6") ], 
                        process = JAC.Photo, 
                        processSettings=PhotoIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [50., 60., 70., 80., 90., 100.], false, true, true, true, 
                                        false, Tuple{Int64,Int64}[], ExpStokes(1., 0., 0.)) )

wb = perform(wa)
JAC.define("print summary: close", "")


#==
                        processSettings=PhotoIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [3000., 3200.], false, true, true, true, 
                                        true, Tuple{Int64,Int64}[(1,1)]) )
==#
