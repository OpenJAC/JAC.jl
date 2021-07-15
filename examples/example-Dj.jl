println("Dj) Test of the MultiPhotonDeExcitation module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-MultiPhotonDeExcitation.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

## grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
name = "two-photon absorption or emission of lithium-like xenon"

if  false
    # Green function for two-photon absorption or emission of lithium-like xenon
    refConfigs       = [Configuration("[He] 2p")]
    levelSymmetries  = [LevelSymmetry(1//2, Basics.minus), LevelSymmetry(3//2, Basics.minus), 
                        LevelSymmetry(1//2, Basics.plus),  LevelSymmetry(3//2, Basics.plus)]
    greenSettings    = GreenSettings(3, [0, 1], 0.01, true, LevelSelection())

    greenRep         = Representation(name, Nuclear.Model(54.), grid, refConfigs, 
                                      GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 3, greenSettings) )  
                                   ## GreenExpansion( AtomicState.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                   ## GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron()
    greenOut         = generate(greenRep, output=true)
    mpGreen          = greenOut["Green channels"]
    
elseif  false
    #
    println("Calculate two-photon absorption cross sections by monochromatic and equally-polarized photons from the same beam")
    mpSettings       = MultiPhotonDeExcitation.Settings(MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic(), 
                                                        [E1,M1], [UseCoulomb], mpGreen, 0, true, LineSelection())
    wa = Atomic.Computation(Atomic.Computation(), name=name, grid=grid, nuclearModel=Nuclear.Model(54.), 
                            initialConfigs= [Configuration("1s^2 2s")],
                            finalConfigs  = [Configuration("1s^2 3s")], 
                            processSettings= mpSettings  )
    wb = perform(wa)
    
elseif  true
    #
    println("Calculate two-photon emission cross sections: energy-differential + total cross sections")
    mpSettings       = MultiPhotonDeExcitation.Settings(MultiPhotonDeExcitation.TwoPhotonEmission(), 
                                                        [E1,M1], [UseCoulomb], mpGreen, 4, true, LineSelection())
    
    wc = Atomic.Computation(Atomic.Computation(), name=name, grid=grid, nuclearModel=Nuclear.Model(54.), 
                            initialConfigs= [Configuration("1s^2 3s")],
                            finalConfigs  = [Configuration("1s^2 2s")], 
                            processSettings= mpSettings    )
    wd = perform(wc)
    
end

setDefaults("print summary: close", "")


