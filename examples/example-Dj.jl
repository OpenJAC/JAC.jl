println("Dj) Test of the MultiPhotonDeExcitation module with ASF from an internally generated initial- and final-state multiplet.")

#
setDefaults("print summary: open", "zzz-MultiPhotonDeExcitation.sum")

if  false
    #
    println("Calculate two-photon absorption cross sections by monochromatic and equally-polarized photons from the same beam")
    name          = "Lithium 1s^2 2s ground configuration"
    refConfigs    = [Configuration("[He] 2p")]
    greenSettings = GreenSettings(5, [0, 1, 2], 0.01, true, false, Int64[])

    wa          = Representation(name, Nuclear.Model(54.), Radial.Grid(true), refConfigs, 
                                 GreenExpansion( AtomicState.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                                 ## GreenExpansion( AtomicState.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                 ## GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                                                 [LevelSymmetry(1//2, Basics.minus), LevelSymmetry(3//2, Basics.minus)], 3, greenSettings) )
    println(wa)
    wb            = generate(wa, output=true)
    greenChannels = wb["Green channels"]
    
elseif  false
    #
    println("Calculate two-photon absorption cross sections by monochromatic and equally-polarized photons from the same beam ... continue")
    
    wc = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(54.), 
                            initialConfigs=[Configuration("1s^2 2s")],
                            finalConfigs  =[Configuration("1s^2 3s")], 
                            process = JAC.MultiPhotonDE, 
                            processSettings=MultiPhotonDeExcitation.Settings(MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic(), 
                                                                    [E1,M1], [UseCoulomb], greenChannels, true, false, Tuple{Int64,Int64}[])  )
    wd = perform(wc)
    
elseif  false
    #
    println("Calculate two-photon emission cross sections: energy-differential + total cross sections")
    name          = "Lithium 1s^2 2s ground configuration"
    refConfigs    = [Configuration("[He] 2p")]
    greenSettings = GreenSettings(5, [0, 1, 2], 0.01, true, false, Int64[])

    wa          = Representation(name, Nuclear.Model(54.), Radial.Grid(true), refConfigs, 
                                 GreenExpansion( AtomicState.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                                 ## GreenExpansion( AtomicState.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                 ## GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                                                 [LevelSymmetry(1//2, Basics.minus), LevelSymmetry(3//2, Basics.minus)], 3, greenSettings) )
    println(wa)
    wb            = generate(wa, output=true)
    greenChannels = wb["Green channels"]
    
elseif  true
    #
    println("Calculate two-photon emission cross sections: energy-differential + total cross sections ... continued")
    
    wc = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(54.), 
                            initialConfigs=[Configuration("1s^2 3s")],
                            finalConfigs  =[Configuration("1s^2 2s")], 
                            process = JAC.MultiPhotonDE, 
                            processSettings=MultiPhotonDeExcitation.Settings(MultiPhotonDeExcitation.TwoPhotonEmission(), 
                                                                    [E1,M1], [UseCoulomb], greenChannels, true, false, Tuple{Int64,Int64}[])  )
    wd = perform(wc)
    
end

setDefaults("print summary: close", "")


