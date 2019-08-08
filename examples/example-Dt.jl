#
println("Dt) Test of the InternalConversion module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-InternalConversion.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(36.);
                        initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                        process = JAC.InternalConv, 
                        processSettings=InternalConversion.Settings() )

wb = perform(wa)
setDefaults("print summary: close", "")

