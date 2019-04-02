

println("Dh) Test of the PhotoExcitationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
#
JAC.define("print summary: open", "zzz-PhotoExcitationAutoIon.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.); grid=JAC.Radial.Grid("grid: by given parameters"; rnt = 2.0e-5, h = 5.0e-2, hp = 1.5e-2, NoPoints = 600), 
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        intermediateConfigs=[Configuration("1s 2s^2"), Configuration("1s 2p^2")],
                        finalConfigs=[Configuration("1s^2")], 
                        process = JAC.PhotoExcAuto, 
                        processSettings=PhotoExcitationAutoion.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, false, Tuple{Int64,Int64}[], 2)  )

wb = perform(wa)
JAC.define("print summary: close", "")

##                      processSettings=PhotoExcitationAutoion.Settings([E1, M1], [JAC.UseCoulomb, JAC.UseBabushkin], true, false, Tuple{Int64,Int64}[], 2)  )
