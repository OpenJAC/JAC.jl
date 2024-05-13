#
println("Ai) Apply & test for restricted-active-space (RAS) expansions.")

if  true
    # Last successful:  unknown ...
    # Compute 
    name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state"
    refConfigs  = [Configuration("[He] 2s^2")]
    rasSettings = RasSettings([1], 24, 1.0e-6, CoulombInteraction(), LevelSelection(true, indices=[1,2,3,4,5,6]) )
    from        = [Shell("2s")]
    #
    frozen      = [Shell("1s")]
    to          = [Shell("2s"), Shell("2p")]
    step1       = RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    append!(frozen, [Shell("2s"), Shell("2p")])
    append!(to,     [Shell("3s"), Shell("3p"), Shell("3d")])
    step2       = RasStep(step1; seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    append!(frozen, [Shell("3s"), Shell("3p"), Shell("3d")])
    append!(to,     [Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")])
    step3       = RasStep(step2, seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))
    #
    wa          = Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs, 
                                RasExpansion([LevelSymmetry(0, Basics.plus), LevelSymmetry(1, Basics.minus)], 4, [step1, step2, step3], rasSettings) )
    println("wa = $wa")

    wb = generate(wa, output=true)
    #
end

