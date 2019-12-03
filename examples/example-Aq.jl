
#
println("Aq) Test of the Green(function) expansion.")

name          = "Lithium 1s^2 2s ground configuration"
refConfigs    = [Configuration("[He] 2s")]
greenSettings = GreenSettings(5, [0, 1, 2], 0.01, true, false, Int64[])
#
wa          = Representation(name, Nuclear.Model(8.), Radial.Grid("grid: exponential"), refConfigs, 
                             ## GreenExpansion( Atomic.SingleCSFwithoutCI(), Basics.DeExciteSingleElectron(), 
                             ## GreenExpansion( Atomic.CoreSpaceCI(), Basics.DeExciteSingleElectron(), 
                                GreenExpansion( Atomic.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                                             [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.plus)], 3, greenSettings) )
println(wa)

wb = generate(wa, output=true)



