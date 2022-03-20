#
println("Da) Test of the PhotoEmission module with ASF from an internally generated initial- and final-state multiplet.")
#
# setDefaults("unit: rate", "a.u.")
setDefaults("unit: energy", "Kayser")
setDefaults("unit: rate", "1/s")
setDefaults("print summary: open", "zzz-radiative.sum")

if false

    grid = Radial.Grid(true)
    ## grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5, h = 1.0e-2, hp = 0., rbox = 5.0)
    setDefaults("standard grid", grid)
    defaultsSettings = PhotoEmission.Settings()
    photoSettings = PhotoEmission.Settings(defaultsSettings, multipoles=[E1], gauges=[UseCoulomb,UseBabushkin], printBefore=true)
    
    comp = Atomic.Computation(Atomic.Computation(), name="Energies and Einstein coefficients for the spectrum Fe X",  
              grid=grid, nuclearModel=Nuclear.Model(26.);
              initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
              finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")], 
              processSettings = photoSettings ); 
    @show comp
    perform(comp)          
              
elseif false
    # LineSelection(true, indexPairs=[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)])
    wa = Atomic.Computation(Atomic.Computation(), name="xx",  grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(3.),
                            initialConfigs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), Configuration("1s^2 3s"),  
                                            Configuration("1s^2 3p"), Configuration("1s^2 3d"), Configuration("1s^2 4s"),  
                                            Configuration("1s^2 4p"), Configuration("1s^2 4d")],
                            finalConfigs  =[Configuration("1s^2 2s"), Configuration("1s^2 2p"), Configuration("1s^2 3s"),  
                                            Configuration("1s^2 3p"), Configuration("1s^2 3d"), Configuration("1s^2 4s"),  
                                            Configuration("1s^2 4p"), Configuration("1s^2 4d")], 
                            processSettings= PhotoEmission.Settings([E1], [UseCoulomb,JAC.UseBabushkin], false, true, LineSelection(), 0., 0., 10000. ) )
    wb = @time( perform(wa) )
elseif false
    # LineSelection(true, indexPairs=[(5,0), (7,0), (10,0), (11,0), (12,0), (13,0), (14,0), (15,0), (16,0)])
    wa = Atomic.Computation(Atomic.Computation(), name="xx",  grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(90.),
                            initialConfigs=[Configuration("1s^2 2s 2p")],
                            finalConfigs  =[Configuration("1s^2 2s^2")], 
                            processSettings= PhotoEmission.Settings([E1,M1,E2,M2], [UseCoulomb,JAC.UseBabushkin], false, true, LineSelection(), 0., 0., 10000. ) )
    wb = @time( perform(wa) )
elseif false
    # Test: Lifetimes of various helium- and lithium-like ions; cf. Figure 8.1 in section 8.1.a
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    pSettings = PhotoEmission.Settings([E1,M1], [UseCoulomb,JAC.UseBabushkin], false, true, CorePolarization(), LineSelection(), 0., 0., 10000. )
    wa = Atomic.Computation(Atomic.Computation(), name="xx",  grid=grid, nuclearModel=Nuclear.Model(14.),
                            initialConfigs=[Configuration("1s 2s"), Configuration("1s 2p")],
                            finalConfigs  =[Configuration("1s^2")], 
                            processSettings=pSettings)
    wb = @time( perform(wa) )
elseif true 
    # Test: Lifetimes for lithium-like ions; cf. Figure 8.1 in section 8.1.a
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    pSettings = PhotoEmission.Settings([E1,M1], [UseCoulomb,JAC.UseBabushkin], false, true, CorePolarization(), LineSelection(), 0., 0., 10000. )
    wa = Atomic.Computation(Atomic.Computation(), name="Lifetimes for lithium-like ions", 
                            grid=grid, nuclearModel=Nuclear.Model(12.), 
                            initialConfigs = [Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs   = [Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                            processSettings = pSettings );
    wb = @time( perform(wa) )
    #
elseif true 
    # Test for Christophe Hoffmeister
    defaultsSettings = PhotoEmission.Settings()
    photoSettings    = PhotoEmission.Settings(defaultsSettings, multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], printBefore=true)
    wa = Atomic.Computation(Atomic.Computation(), name="Energies and Einstein coefficients for neon 2s^-1 3p --> 2p^6 ^1S_0", 
                            grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(10.1), 
                            initialConfigs = [Configuration("1s^2 2s 2p^6 3p")],
                            finalConfigs   = [Configuration("1s^2 2s^2 2p^6")], 
                            process = Radiative(), processSettings = photoSettings );
    wb = @time( perform(wa) )
    #
end
setDefaults("print summary: close", "")
