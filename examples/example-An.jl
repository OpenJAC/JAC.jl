#
println("An) Test of a mean-field basis and configuration-interaction (CI) expansion.")
#

name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
refConfigs  = [Configuration("[He] 2s^2 2p^4")]
mfSettings  = MeanFieldSettings()
#
wa          = Representation(name, Nuclear.Model(8.), Radial.Grid("grid: exponential"), refConfigs, 
                             MeanFieldBasis(mfSettings) )
println("wa = $wa")

wb = generate(wa, output=true)

orbitals    = wb["mean-field basis"].orbitals
ciSettings  = CiSettings(true, false, Int64[], false, LevelSymmetry[] )
from        = [Shell("2s")]
#
frozen      = [Shell("1s")]
to          = [Shell("2s"), Shell("2p")]
excitations = RasStep()
#             RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=deepcopy(frozen))
#
wc          = Representation(name, Nuclear.Model(8.), Radial.Grid("grid: exponential"), refConfigs, 
                             CiExpansion(orbitals, excitations, ciSettings) )
println("wc = $wc")

wd = generate(wc, output=true)

