
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
    ## THE GREEN SPACE MUST BE BUILT ON THE GRID THE COMPUTATION USES, and until 01-Nov-2026 this line passed
    ## `Radial.Grid(true)` -- the global STANDARD grid, 392 points out to 614 a.u. -- while branch 2 computes on
    ## the grid defined at the top of this file, 1302 points out to 10 a.u.  The orbitals come back with
    ## `useStandardGrid = true` and an EMPTY stored grid, so their radial arrays were then read index by index
    ## against a mesh they do not belong to.  MEASURED: <a|a> on the computation's own grid came out 0.3581
    ## (1s_1/2), 0.1149 (2s_1/2), 0.1242 / 0.1241 (2p_1/2 / 2p_3/2) where it must be 1, i.e. norms wrong by
    ## factors of 3 to 9.  That is not two legitimate orbital sets differing -- it is ONE set mis-indexed, and no
    ## choice of transition operator can repair it.  Passing `grid` here restores <a|a> = 1.0000000000 exactly,
    ## with length(P) matching the mesh;  `setDefaults("standard grid", grid)` would do the same.
    wa          = Representation(name, Nuclear.Model(4.01), grid, refConfigs, MeanFieldMultiplet(mfSettings) )
    wb          = generate(wa, output=true)
    #
elseif true
    # Last visit:  01-Nov-2026 (NEEDS AN EARLIER BRANCH OF THIS FILE)
    #
    # THE maxKappa = 1 NUMBERS THIS COMMENT USED TO QUOTE -- 400 eV, 2.998258e+02 barn Coulomb and
    # 1.701679e+02 Babushkin -- ARE WITHDRAWN, 01-Nov-2026.  They were computed with the Green orbitals read on
    # the wrong radial mesh (see the note in branch 1), i.e. with norms wrong by factors of 3 to 9, so they were
    # never a result.  With the Green space on ONE grid, at this branch's own maxKappa = 2 and both energies:
    #
    #     omega [eV]      Coulomb [barn]    Babushkin [barn]    Cou/Bab
    #        400           8.024051e+02       1.268623e+03       0.632
    #        500           7.132792e+02       2.485340e+03       0.287
    #
    # and for the record, the same run BEFORE the grid was fixed gave 7.802171e+02 / 5.444154e+02 (1.433) and
    # 9.421547e+02 / 1.107223e+03 (0.851).  So the repair moves BABUSHKIN by a factor ~2.3 and Coulomb barely,
    # and the gauge ratio gets WORSE rather than better.  **That is reported rather than explained away:** the
    # fix is right on its own terms, since an orbital whose <a|a> reads 0.358 on the mesh it is integrated on is
    # wrong whatever it is multiplied by -- but the worsening says this calculation's LARGER errors lie
    # elsewhere (maxKappa = 2, three ionic cores, two sharings), and until those are addressed its gauge ratio
    # is not a figure of merit that can be read in either direction.
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


