
#
#=
Example: Test of the Dielectronic module      
=#
println("Df) Test of the Dielectronic module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
#
JAC.define("print summary: open", "zzz-Dielectronic.sum")
JAC.define("method: continuum, Galerkin")
JAC.define("method: normalization, pure sine")   ## JAC.define("method: normalization, pure Coulomb")    JAC.define("method: normalization, pure sine")

wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 2.0e-2, NoPoints = 600),
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                        finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                        process = JAC.Dierec, 
                        processSettings=Dielectronic.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, 
                                                              true, Tuple{Int64,Int64,Int64}[(1,1,0), (1,2,0), (1,3,0)], 0., 0., 0., "Coulomb")  )

wb = perform(wa)
JAC.define("print summary: close", "")


