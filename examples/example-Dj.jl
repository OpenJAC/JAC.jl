println("Dj) Test of the MultiPhotonDeExcitation module with ASF from an internally generated initial- and final-state multiplet.")

@warn("\n\n !!! This example does not work properly at present !!! \n\n")
#
setDefaults("print summary: open", "zzz-MultiPhotonDeExcitation.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(26.), 
                        initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")],
                        intermediateConfigs=[Configuration("1s 2s^2 2p"), Configuration("1s 2s 2p^2") ],
                        finalConfigs  =[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p") ], 
                        process = JAC.MultiPhotonDE, 
                        processSettings=MultiPhotonDeExcitation.Settings()  )
#

wb = perform(wa)
setDefaults("print summary: close", "")


