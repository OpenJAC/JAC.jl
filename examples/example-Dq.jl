
#
println("Dq)  Test of the TwoElectronOnePhoton module with ASF from an internally generated initial and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-TwoElectronOnePhoton.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")


grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

if  false
    # Two-electron-one-photon transition of Li-like Ar: 1s 2p 3d  -->  1s^2 2p: Generation of MeanFieldMultiplet
    name        = "Li-like Xe"
    refConfigs  = [Configuration("1s^2 2s 3s"), Configuration("1s^2 2s 3p")]
    mfSettings  = MeanFieldSettings()
    #
    wa          = Representation(name, Nuclear.Model(54.01), Radial.Grid(true), refConfigs, MeanFieldMultiplet(mfSettings) )
    wb          = generate(wa, output=true)
    #
elseif  true
    # Two-electron-one-photon transition of Li-like Ar: 1s 2p 3d  -->  1s^2 2p: Computation of rates
    nMultiplet    = wb["mean-field multiplet"]
    rydbergShells = Basics.generateShellList(6,  6, "p") 
    teopSettings  = TwoElectronOnePhoton.Settings([E1], [UseCoulomb,JAC.UseBabushkin], true, LineSelection(), -200., CoulombInteraction(), nMultiplet)
    wc = Atomic.Computation(Atomic.Computation(), name="TEOP", grid=grid, nuclearModel=Nuclear.Model(54.01), 
                            initialConfigs = [Configuration("1s^2 3s 3p")],
                            finalConfigs   = [Configuration("1s^2 2s^2")],  
                            processSettings= teopSettings )
    ## @show wc
    wd = perform(wc)
elseif  true
    ## setDefaults("standard grid", grid)
    defaultsSettings = PhotoEmission.Settings()
    photoSettings = PhotoEmission.Settings(defaultsSettings, multipoles=[E1], gauges=[UseCoulomb,UseBabushkin], printBefore=true)
    
    we = Atomic.Computation(Atomic.Computation(), name="OEOP",  grid=grid, nuclearModel=Nuclear.Model(54.01),
                            initialConfigs = [Configuration("1s^2 3s 3p")],
                            finalConfigs   = [Configuration("1s^2 2s 3s"), Configuration("1s^2 2p 3p")], 
                            processSettings = photoSettings ); 
    perform(we)          
elseif  false
    # Example Marek Pajek: TEOP of Xe^34+ [Ar] + (10s + 10p + 10d ...) -->  [Ar] 3d
    # branching fraction:  TEOP / OEOP = 3.8 * 10^-4 +- 9.6 * 10^-5
    # ... should be reasonable easy if the upper example works reasonable; perhaps first for the helium-like analogue
end
    
