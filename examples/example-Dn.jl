println("Dn) Test of the ImpactExcitation module with ASF from an internally generated initial- and final-state multiplet.")
#
JAC.define("print summary: open", "zzz-ImpactExcitation.sum")
wa = Atomic.Computation("xx",  Nuclear.Model(36.);
                        initialConfigs=[Configuration("1s^2 2s^2 2p^6")],
                        finalConfigs  =[Configuration("1s^2 2s^2 2p^5"), Configuration("1s^2 2s 2p^6") ], 
                        process = Eimex, 
                        processSettings=ImpactExcitation.Settings([100., 200.], false, true, false, Tuple{Int64,Int64}[], 2, 0.) )

wb = perform(wa)


