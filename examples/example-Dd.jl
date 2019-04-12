
#
println("Dd) Test of the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.")
#
JAC.define("print summary: open", "zzz-PhotoRecombination.sum")
JAC.define("method: continuum, asymptotic Coulomb")  ## JAC.define("method: continuum, Galerkin")
JAC.define("method: normalization, pure Coulomb")   ## JAC.define("method: normalization, pure Coulomb")    JAC.define("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(12.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600),
                        initialConfigs=[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^6")], 
                        process = JAC.Rec, 
                        processSettings=PhotoRecombination.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [10., 30., 50.], [0.], 
                                        false, true, true, true, true, Tuple{Int64,Int64}[(1,1)]) )

wb = perform(wa)
JAC.define("print summary: close", "")

