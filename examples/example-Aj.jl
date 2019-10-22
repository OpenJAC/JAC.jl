#
println("Aj) Test of restricted-active-space (RAS) computations.")
#

name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state"
refConfigs  = [Configuration("[He] 2s^2")]
rasSettings = RasSettings([1], 24, 1.0e-6, false, true, [1,2,3] )
frozen      = [Shell("1s")]
from        = [Shell("2s")]
#
to          = [Shell("2s"), Shell("2p")]
step1       = RasStep(seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=frozen)
#
append!(to, [Shell("3s"), Shell("3p"), Shell("3d")])
step2       = RasStep(seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=frozen)
#
append!(to, [Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")])
step3       = RasStep(seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=frozen)
#
wa          = RasComputation(name, Nuclear.Model(4.), Radial.Grid("grid: exponential"), 
                             refConfigs, LevelSymmetry(0, Basics.plus), 4, 15, [step1, step2, step3], rasSettings)

wb = perform(wa)

