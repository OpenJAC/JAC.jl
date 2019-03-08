#
println("Da) Test of the Radiative module with ASF from an internally generated initial- and final-state multiplet.")
#
# JAC.define("unit: rate", "a.u.")
JAC.define("unit: rate", "1/s")
JAC.define("print summary: open", "zzz-radiative.sum")
#== wa = Atomic.Computation("xx",  Nuclear.Model(36.);
                        initialConfigs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                        finalConfigs  =[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")], 
                        process = JAC.RadiativeX, 
                        processSettings=Radiative.Settings([E1, M1, E2, M2], [JAC.UseCoulomb, JAC.UseBabushkin], true, true, 
                           true, Tuple{Int64,Int64}[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)], 0., 0., 10000. ) )  ==#
                           
wa = Atomic.Computation("Energies and Einstein coefficients for the spectrum Fe X",  Nuclear.Model(26.);
                        initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
                        finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")], 
                        process = RadiativeX, 
                        processSettings = Radiative.Settings([E1, M1, E2, M2], [JAC.UseCoulomb, JAC.UseBabushkin], 
                                          false, true, false, Tuple{Int64,Int64}[], 0., 0., 10000. ) );
wb = @time( perform(wa) )
JAC.define("print summary: close", "")


## , Configuration("[Ne] 3s^2 3p^4 3d")  Configuration("[Ne] 3s 3p^6")
