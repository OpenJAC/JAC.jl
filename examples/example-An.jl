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

orbitals    = wa.orbitals
ciSettings  = CiSettings(true, false, Int64[], false, LevelSymmetry[] )
from        = [Shell("2s")]
#
frozen      = [Shell("1s")]
to          = [Shell("2s"), Shell("2p")]
excitations = RasStep()
#             RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=deepcopy(frozen))
#
wb          = Representation(name, Nuclear.Model(8.), Radial.Grid("grid: exponential"), refConfigs, 
                             CiExpansion(orbitals, excitations, ciSettings) )
println("wb = $wb")

wb = generate(wa, output=true)

