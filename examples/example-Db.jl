
println("Db) Apply & test the PhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoExcitation.sum")

if  true
    # Last successful:  12May2024
    # Compute the photoexcitation cross sections for neutral lithium
    pSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb, UseBabushkin], printBefore=true) 
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(3.), 
                            initialConfigs  = [Configuration("1s^2 2s")],
                            finalConfigs    = [Configuration("1s^2 2p")], 
                            processSettings = pSettings )

    wb = perform(wa)
    #
elseif true
    # Last successful:  12May2024
    # Test of absorption f-values for 1s^2 --> 1s2p ^1P_1 resonance  in the helium isoelectronic sequence; cf. Table 8.1 in Section 8.1.a
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    pSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb, UseBabushkin], printBefore=true) 
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(20.), 
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s 2p")], 
                            processSettings = pSettings )
    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")
