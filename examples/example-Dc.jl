#
println("Dc) Test of the PhotoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoIonization.sum")
setDefaults("method: continuum, asymptotic Coulomb")  ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure Coulomb")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(false), nuclearModel=Nuclear.Model(20.), 
                        initialConfigs=[Configuration("[Ne] 3s^2 3p^6")],
                        finalConfigs  =[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6") ], 
                        process = JAC.Photo, 
                        processSettings=PhotoIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [50., 60., 70., 80., 90., 100.], false, true, true, true, 
                                        false, Tuple{Int64,Int64}[], JAC.ExpStokes(1., 0., 0.)) )

wb = perform(wa)
setDefaults("print summary: close", "")

