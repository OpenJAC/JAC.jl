
#
println("Dd) Test of the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoRecombination.sum")
setDefaults("method: continuum, asymptotic Coulomb")  ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure Coulomb")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(12.); grid=JAC.Radial.Grid(false),
                        initialConfigs=[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^6")], 
                        process = JAC.Rec, 
                        processSettings=PhotoRecombination.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [10., 30., 50.], [0.], 
                                        false, true, true, true, true, Tuple{Int64,Int64}[(1,1)]) )

wb = perform(wa)
setDefaults("print summary: close", "")

