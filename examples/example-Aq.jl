
#
println("Aq) Test of the Green(function) expansion.")

if  false
    name          = "Lithium 1s^2 2s ground configuration"
    refConfigs    = [Configuration("[He] 2s")]
    greenSettings = GreenSettings(5, [0, 1, 2], 0.01, true, LevelSelection())
    #
    wa          = Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, 
                                ## GreenExpansion( Atomic.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                                ## GreenExpansion( Atomic.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                    GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                                                    [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.plus)], 3, greenSettings) )
    println(wa)

    wb = generate(wa, output=true)
    #
elseif true
    # Test example for MS "Approximate atomic Green functions" in molecules (2021)
    name            = "Approximate Green function for neon-like ions."
    refConfigs      = [Configuration("1s^2 2s^2 2p^6")]
    levelSymmetries = [LevelSymmetry(1, Basics.minus), LevelSymmetry(2, Basics.minus)]
    greenSettings   = GreenSettings(15, [0, 1, 2, 3], 0.01, true, LevelSelection())
    greenExpansion  = GreenExpansion(AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                                 levelSymmetries, 10, greenSettings)
    #
    rep = Representation(name, Nuclear.Model(18.), Radial.Grid(true), refConfigs, greenExpansion)
    #
    generate(rep)
    #
end



