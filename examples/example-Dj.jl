
println("Dj) Test of the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoDoubleIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  true
    # Last visit:  01-Nov-2026
    # Last successful:  unknown ...
    # Photo-double ionization of he-like Be: 1s^2   -->  epsilon1 kappa1 + epsilon2 kappa2:
    #   Generation of the MeanFieldMultiplet of IONIC CORES.
    #
    # THE GREEN SPACE IS THE (N-1)-ELECTRON IONIC CORES, AND IT WAS NOT, UNTIL 01-Nov-2026.  This branch built
    #   it as `AddElectrons(1, quasiShells)` applied to [1s], i.e. the TWO-electron configurations 1s4s and
    #   1s4p -- the "intermediate levels" that PhotoDoubleIonization.Settings' gMultiplet field then asked for
    #   in words.  But PhotoDoubleIonization.generateIntermediateLevels generates a continuum orbital IN THE
    #   FIELD OF each gMultiplet level and ATTACHES it, so those became THREE-electron intermediates against a
    #   TWO-electron initial state, and every junction of the second-order sum paired unequal electron numbers.
    #   Measured 01-Nov-2026: initial 2, Green 2, Green + one continuum electron 3.  That is why this branch has
    #   never carried a `Last successful` date, and the module now REFUSES the mismatch instead of returning a
    #   plausible number.  The 4l quasi-shells were a hand discretisation of the continuum and are superseded:
    #   generateIntermediateLevels carries its own energy grid (NoIntermediateEnergies).
    name        = "He-like Be"
    refConfigs  = [Configuration("1s"), Configuration("2s"), Configuration("2p")]
    mfSettings  = MeanFieldSettings()
    #
    wa          = Representation(name, Nuclear.Model(4.01), Radial.Grid(true), refConfigs, MeanFieldMultiplet(mfSettings) )
    wb          = generate(wa, output=true)
    #
elseif true
    # Last visit:  01-Nov-2026 (NEEDS AN EARLIER BRANCH OF THIS FILE)   ## with the Green space corrected to the
    #                            (N-1)-electron cores it now COMPLETES with no warning, at maxKappa = 1 and one
    #                            photon energy: 400 eV gives 2.998258e+02 barn Coulomb, 1.701679e+02 Babushkin.
    # Last successful:  unknown ...   ## and it stays unknown: a Coulomb/Babushkin ratio of 1.76 is not a
    #                            converged cross section.  Three ionic cores and maxKappa = 1 is a demonstration
    #                            of the pipeline, not a number to quote;  the convergence study is open work.
    # Photo-double ionization of he-like Be: 1s^2   -->  epsilon1 kappa1 + epsilon2 kappa2: cross sections
    # `wb` is produced by an earlier branch of this file (line 21).
    if  !@isdefined(wb)   error("Run the earlier branch of this file first; it defines `wb`.")   end
    nMultiplet     = wb[ResultKeys.MeanFieldMultiplet]
    # The KEYWORD form, so that a new Settings field cannot silently break this call -- which had already happened:
    # the positional version passed `quasiShells` third, where `photonEnergies` now sits, and a field of that name
    # no longer exists.  A second field was added on 30-Sep-2026 (NoIntermediateEnergies) and would have been the
    # eleventh of twelve.
    doubleSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                         multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[400., 500.],
                         NoEnergySharings=2, maxKappa=2, calcDifferentialCs=true, printBefore=true,
                         lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                         eeInteraction=CoulombInteraction(), gMultiplet=nMultiplet )
    
    wc = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(4.01), 
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s^0")],
                            processSettings = doubleSettings )

    wd = perform(wc)
    #
end
#
setDefaults("print summary: close", "")


