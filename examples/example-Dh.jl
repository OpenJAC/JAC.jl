

println("Dh) Test of the PhotoExcitationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
#
setDefaults("print summary: open", "zzz-PhotoExcitationAutoIon.sum")
grid = JAC.Radial.Grid(true)
grid = JAC.Radial.Grid(grid, rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, NoPoints = 600)
wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(26.), 
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        intermediateConfigs=[Configuration("1s 2s^2"), Configuration("1s 2p^2")],
                        finalConfigs=[Configuration("1s^2")], 
                        process = PhotoExcAuto(), 
                        processSettings=PhotoExcitationAutoion.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, PathwaySelection(), 2)  )

wb = perform(wa)
setDefaults("print summary: close", "")
