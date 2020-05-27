
println("Dq) Test of the RadiativeAuger module with ASF from an internally generated initial- and final-state multiplet.")
#
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                        process = RAuger(), 
                        processSettings=RadiativeAuger.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], 6, true, false, 
                                                                Tuple{Int64,Int64}[], 0., 1.0e6, 2)  )

wb = perform(wa)

