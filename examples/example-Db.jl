
println("Db) Test of the PhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoExcitation.sum")
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(3.), 
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        finalConfigs  =[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                        process = PhotoExc(), 
                        processSettings=PhotoExcitation.Settings([E1, M1], [UseCoulomb, UseBabushkin], true, true, true, true, 
                                                                 LineSelection(), 0., 0., 1.0e6, JAC.ExpStokes(0., 0., 0.) ) )

wb = perform(wa)
setDefaults("print summary: close", "")


#     [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")], 

