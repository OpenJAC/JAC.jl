
println("Di) Test of the RayleighCompton module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
#
setDefaults("print summary: open", "zzz-RayleighCompton.sum")
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                        finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                        process = JAC.Compton, 
                        processSettings=RayleighCompton.Settings()  )

wb = perform(wa)
setDefaults("print summary: close", "")


## [E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], false, false, false, true, true, Tuple{Int64,Int64,Int64}[(1,1,0), (1,2,0), (1,2,0)], 0., 0., 0., "Coulomb"
