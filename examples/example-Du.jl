#
println("Du) Test of the CoulombIonization module with ASF from an internally generated initial- and final-state multiplet.")

@warn("\n\n !!! This example does not work properly at present !!! \n\n")
#
setDefaults("print summary: open", "zzz-CoulombIonization.sum")
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(36.), 
                        initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                        process = JAC.Coulion, 
                        processSettings=CoulombIonization.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], [3000., 4000.], false, true, 
                                        false, Tuple{Int64,Int64}[]) )

wb = perform(wa)
setDefaults("print summary: close", "")

