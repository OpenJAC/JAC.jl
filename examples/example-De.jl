#
println("De) Test of the AutoIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-AutoIonization.sum")
setDefaults("method: continuum, asymptotic Coulomb")  ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")       ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(12.); grid=JAC.Radial.Grid(false), 
                        initialConfigs  =[Configuration("1s 2s^2 2p^6")],
                        finalConfigs    =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], process = JAC.Auger,
                        processSettings = AutoIonization.Settings(true, true, true, Tuple{Int64,Int64}[(1,0)], 0., 1.0e6, 2, "Coulomb") )

wb = perform(wa)
setDefaults("print summary: close", "")


