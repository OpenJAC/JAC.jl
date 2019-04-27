
#
#=
Example: Test of the Dielectronic module      
=#
println("Df) Test of the Dielectronic module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
#
setDefaults("print summary: open", "zzz-Dielectronic.sum")
setDefaults("method: continuum, asymptotic Coulomb")  ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")   ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600),
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                        finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                        process = JAC.Dierec, 
                        processSettings=Dielectronic.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                                              true, Tuple{Int64,Int64,Int64}[(1,1,0), (1,2,0), (1,3,0)], 0., 0., 0., "Coulomb")  )

wb = perform(wa)
setDefaults("print summary: close", "")


