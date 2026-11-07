
println("Dj) Test of the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-PhotoDoubleIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  true
    # Last visit:  07-Nov-2026
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
elseif true
    # Last visit:  05-Nov-2026 (RE-DATED -- the cross sections below are NOT those of 02-Nov, see WHY THE NUMBERS MOVED)
    # Last reproduced:  05-Nov-2026   ## threshold 371.4828 eV;  400 eV Coulomb 8.659709e+02 / Babushkin
    #                            2.874448e+03 barn;  500 eV Coulomb 5.411607e+02 / Babushkin 1.637261e+03
    #                            barn  (13 of 13 numbers below, 94 s)
    # Last successful:  unknown ...   ## the pipeline is complete here and the ENERGIES are right, but the cross
    #                            sections are NOT CONVERGED -- the PARTIAL-WAVE SUM is the leading error, see
    #                            WHAT IS NOT CONVERGED below.  The gauge split is NOT the reason;  see below.
    #
    # Photo-double ionization of he-like Be, 1s^2 --> epsilon1 kappa1 + epsilon2 kappa2, with ALL THREE multiplets
    # -- initial, final and Green -- carried on ONE orbital set.  This is the branch the file was missing:  branches
    # 1 and 2 above run the same physics with the Green space standing on its own mean field, and a second-order
    # amplitude is only defined on one orthonormal set.
    #
    # WHY ONE SET.  <r|H-E|s> = <r|h_D|s> + <r|V_ee|s> - E <r|s>, and ManyElectron.matrixElement_HminusE drops the
    # third term on an orthonormality assumption that holds ACROSS two independent SCF sets only by accident.  At a
    # junction of the second-order sum where the two CSFs differ in ONE orbital, the dropped term is not small:
    # measured 02-Nov-2026 on this very case, with separate sets <101s_1/2|1s_initial> = -0.0235 and the dropped
    # term came to 2.05x the matrix element it would have corrected.  On the shared set below that overlap falls to
    # -0.003355 and the dropped term to 0.294x -- a factor 7.  Both PhotoDoubleIonization precondition checks
    # (checkGreenSpaceOnGrid and checkSecondOrderConditions) are SILENT in this branch, which no earlier branch of
    # this file has managed.
    #
    # THE SPANNING REQUIREMENT, which is why the SCF runs over three configurations rather than one.  The Green
    # space occupies 1s, 2s and 2p, so the shared set must contain all three;  an SCF on [1s^2] alone yields
    # 1s_1/2 and nothing else, and no choice of transition operator repairs a Green space built on orbitals the
    # initial state has never seen.  The price is paid on the initial state and it is small but not zero:
    #
    #     dedicated SCF on [1s^2]                       -13.684605803997 Ha
    #     frozen CI on the [1s^2, 1s 2s, 1s 2p] set     -13.651740782683 Ha     +0.032865 Ha, i.e. 0.240 %
    #
    # ENERGY BOOKKEEPING -- the part of this branch that IS trustworthy.  The final state is bare Be (1s^0, E = 0),
    # so the double-ionization threshold is simply |E(1s^2)| = 371.4828 eV, and the sharing table must respect it:
    #
    #     omega = 400 eV      eps_1 + eps_2 =   6.02640 +  22.49080 =  28.5172 eV   -->  threshold 371.4828 eV
    #     omega = 500 eV      eps_1 + eps_2 =  27.15890 + 101.35800 = 128.5169 eV   -->  threshold 371.4831 eV
    #
    # consistent to 3e-4 eV, and the two excess energies differ by 99.9997 eV against a photon-energy step of 100.
    # Against the real ion:  the 1/Z expansion of the exact non-relativistic two-electron energy,
    # E = -Z^2 + 0.625 Z - 0.158 Ha, gives -13.732 Ha at Z = 4.01 and so a threshold near 373.7 eV.  The 2.2 eV
    # shortfall decomposes as 0.9 eV of shared-set constraint and about 1.3 eV of the correlation a
    # single-configuration CI cannot hold -- 0.6 % low, which is what this mean field is worth and no worse.
    # AND NOTE THAT Z = 4.01 IS NOT COSMETIC:  dE/dZ = -2Z + 0.625 = -7.4 Ha, so the 0.01 carried by this file's
    # nuclear model moves the threshold by 2.0 eV, more than the shared orbital set costs.
    #
    # EXCHANGE SYMMETRY, a check that is NOT built in.  Each energy sharing is computed separately, with its own
    # pair of continuum orbitals, so the symmetry of dsigma/d eps_1 under eps_1 <--> eps_2 tests the amplitude
    # rather than restating an identity of the code:
    #
    #     400 eV      1.299101e+03  against  1.298984e+03  barn/eV      9e-5 relative
    #     500 eV      2.211048e+02  against  2.210789e+02  barn/eV      1e-4 relative
    #
    # THE TWO GAUGES, WHICH ARE NOT EXPECTED TO AGREE.  Total cross sections on the shared set:
    #
    #     omega [eV]      Coulomb [barn]    Babushkin [barn]    Cou/Bab
    #        400           1.361381e+03       4.492780e+03       0.303
    #        500           1.044199e+03       3.157561e+03       0.331
    #
    # A factor of three, and it is NOT the reason to distrust this branch.  Length and velocity are tied by the
    # commutator identity <a|p|b> = i (E_a - E_b) <a|r|b>, i.e. by a LEVEL difference, whereas a second-order
    # amplitude's photon operator carries the PHYSICAL omega against an intermediate state that is off shell.  The
    # two forms then differ by the factor (E_i - E_n - omega), which vanishes only ON shell, and they are
    # equivalent only through CLOSURE over a COMPLETE intermediate set -- four ionic cores are not that.  See the
    # module's STATUS block and Fortun et al., arXiv:1610.03571:  extended gauge invariance holds only at exact
    # resonance.  So READ THE COULOMB (velocity) COLUMN, and do NOT repair the split by substituting a level
    # difference for omega -- that restores agreement matrix element by matrix element while making the amplitude
    # wrong.
    #
    # What the shared orbital set DID change is the SHAPE of the split:  on separate sets the ratio was 0.632 at
    # 400 eV and 0.287 at 500 eV, strongly energy-dependent, and it is now nearly constant at 0.303 / 0.331.  That
    # is recorded as an observation and not as a figure of merit:  Babushkin is the unstable member of the pair, so
    # the ratio cannot serve as a convergence indicator in either direction.
    #
    # WHAT IS NOT CONVERGED.  With 2 sharings, maxKappa = 2 and four ionic cores this branch demonstrates a pipeline;
    # it does not deliver a cross section.  Prefer velocity, and do NOT read the gauge ratio above as a convergence
    # indicator -- that part stands.  BUT THE TWO REASONS THIS BRANCH GAVE UNTIL 05-Nov-2026 WERE BOTH WRONG, and the
    # correction is worth more than the original claim was.  It read "the sharing integral drifts +29 % between 3 and
    # 21 Gauss-Legendre points and is still rising, so only the first digit carries information", and it blamed
    # Babushkin's instability on the intermediate set.  Measured on helium at 200 eV in branch (d) below, with the
    # continuum amplitude FITTED rather than read off a window:  the sharing quadrature moves +0.58 % over 5 -> 21
    # points, and the intermediate mesh -0.87 % then +0.24 % over 9 -> 17 -> 25 -- non-monotone, hence converged.
    # THE LEADING TRUNCATION IS THE PARTIAL-WAVE SUM, +21.6 % over maxKappa 2 -> 3 and +0.8 % over 3 -> 4.  The
    # +29 % was an artefact of the continuum normalization described next, and it disappeared with it.
    #
    # WHY THE NUMBERS MOVED -- the continuum normalization, found 03-Nov-2026 while chasing the same question in
    # helium and REPAIRED THE SAME DAY, in the very method this file selects.  Continuum.normalizeOrbitalPureSine
    # used to scale an orbital by max|P| over a FIXED 300 grid points, which at hp = 1e-2 is 2.98 a.u. and so half a
    # wavelength only ABOVE 14.9 eV;  below that the window need contain no crest and the scale was wrong by up to
    # 23 %, differing between radial boxes with no pattern.  THIS BRANCH SITS INSIDE THAT RANGE:  its 400 eV sharing
    # puts eps_1 at 6.03 eV, a factor 2.5 below the old floor.  The method now calls Continuum.asymptoticAmplitude,
    # which FITS the asymptotic form instead of reading a maximum off it, and that is why the 02-Nov cross sections
    # above were re-made on 05-Nov:  Coulomb fell 36 % at 400 eV and 48 % at 500 eV, Babushkin by the same factors,
    # so the gauge RATIO is preserved to three digits (0.303 -> 0.3013, 0.331 -> 0.3305) while the SCALE changed.
    # A pure rescale of both gauges is the signature of a normalization repair and not of a physics regression.
    #
    # AND READ THE WARNING THIS BRANCH NOW PRINTS, because the repair also made the defect visible:  the 400 eV
    # 101d_3/2 continuum orbital reports rms misfit / amplitude = 2.59e-01, above the 0.20 the module calls suspect.
    # At 0.2215 Hartree its wavelength is 9.4 a.u. and this file's box is 10.0, so barely one oscillation fits and
    # the fitted asymptotic form is being extracted from almost no asymptotic region.  That is a BOX statement in
    # the sense of Rule 12, and it is the reason the 400 eV row remains the less trustworthy of the two.  What is
    # still missing from src/ is a check of the FORM rather than the amplitude:  the misfit is reported, but nothing
    # refuses an orbital on it, so a catastrophically wrong orbital would still be scaled confidently.
    #
    # BOTH HALVES NOW GO THROUGH THE PUBLIC ROUTE (05-Nov-2026).  Until then the shared orbitals arrived through
    # Representation but the Green space did not -- no Representation carried a caller-supplied orbital set, so it
    # had to be built by a direct Hamiltonian.performCIwithFrozenOrbitals call, outside generate and perform alike.
    # AtomicState.Representation now carries orbitalSource, the same field Atomic.Computation has carried since
    # 02-Nov-2026, so the Green space below is an ordinary generate(Representation) call and a reader can reproduce
    # this branch without leaving the documented interface.  The numbers did not move, and BOTH branches were run
    # to establish it rather than one:  branch (c) reproduced 8.659709e+02 / 2.874448e+03 barn at 400 eV and
    # 5.411607e+02 / 1.637261e+03 at 500 eV, and branch (d) all 7 numbers on its own dated claim line, each
    # digit for digit against the pre-change run.  That is what licenses this note.
    name        = "He-like Be, one orbital set"
    nm          = Nuclear.Model(4.01)
    spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")]
    wa          = Representation(name, nm, grid, spanConfigs, MeanFieldMultiplet(MeanFieldSettings()) )
    wb          = generate(wa, output=true)
    orbitals    = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
    #
    ## the (N-1)-electron ionic cores, as a CI on that SAME set;  generateIntermediateLevels then attaches one
    ## continuum orbital to each, which is what makes the intermediates two-electron states again
    wg          = Representation(name * " -- Green space", nm, grid,
                                 [Configuration("1s"), Configuration("2s"), Configuration("2p")],
                                 MeanFieldMultiplet(MeanFieldSettings()); orbitalSource = Basics.GivenOrbitals(orbitals) )
    gMultiplet  = generate(wg, output=true)[ResultKeys.MeanFieldMultiplet]
    #
    doubleSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                         multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[400., 500.],
                         NoEnergySharings=2, maxKappa=2, calcDifferentialCs=true, printBefore=true,
                         lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                         eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )

    wc = Atomic.Computation(Atomic.Computation(), name=name, grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s^0")],
                            orbitalSource   = Basics.GivenOrbitals(orbitals),
                            processSettings = doubleSettings )

    wd = perform(wc)
    #
elseif true
    # Last visit:  05-Nov-2026
    # Last reproduced:  05-Nov-2026   ## 200 eV, 3 sharings, maxKappa = 1:  Coulomb 7.812346e+02 / Babushkin
    #                            1.876583e+03 barn, Babushkin/Coulomb 2.402;  the quadrature identity
    #                            sigma = 0.5 SUM w*d closes to 6 digits in both gauges;  the eps_1 <--> eps_2
    #                            symmetry of dsigma/d eps_1 is 5.0e-05 (Cou) and 3.2e-05 (Bab);  no continuum
    #                            misfit warning anywhere in the run  (7 of 7 numbers, 196 s)
    # Last successful:  unknown ...   ## and it will stay unknown until the overshoot below is understood:  this
    #                            branch's physics is SOUND in every truncation that has been measured and still
    #                            lands about a fifth above the measurement, which is a MODEL statement.
    #
    # Photo-double ionization of HELIUM at 200 eV on ONE orbital set -- the same construction as branch (c) above,
    # on the one system for which a measurement exists.  This branch is the written record of priority item 6
    # (Nov-2026 numbering), which was carried for four days and is closed onto the challenges list with the
    # question OPEN;  the numbers below are the result of about five hours of machine time that no single run of
    # this branch reproduces, so they are QUOTED here with the harness that produced each one named.
    #
    # WHAT THIS BRANCH ITSELF RUNS is a cheap demonstration -- 3 sharings, maxKappa = 1, 9 intermediate energies --
    # which exercises the whole pipeline in a couple of minutes and is NOT the converged point.  The converged
    # point costs 1.4 h per orbital set.  Raise maxKappa to 3 and the sharings to 9 to reach it.
    #
    # THE GRID IS BUILT HERE AND NOT TAKEN FROM THE TOP OF THE FILE, which matters in the sense of Rule 12:  this
    # file's shared grid has rbox = 10.0, right for he-like Be where the outer electron sees Z_eff = 3.01, and far
    # too small for neutral helium where it sees Z_eff = 1.  A 40 a.u. box was checked against 24 and 60:  the
    # total moves 0.05 % between 24 and 40, so 40 is converged for the TOTAL.  It is not converged for the n = 3
    # members of a widened span, which press against the wall at 60 a.u. as well -- see THE GREEN SPACE below.
    #
    # ===== THE CONVERGED RESULT, and it DOES NOT agree with the measurement =====================================
    #
    # At maxKappa = 3, 9 sharings, 9 intermediate energies, a 40 a.u. box and the n = 2 span (tools/
    # scanPdiHeliumConvergence.jl):
    #
    #     sigma(2+)                 2.8208e-03 Mb velocity        6.7563e-03 Mb length
    #     sigma(2+)/sigma(+)            4.741 % velocity             13.10 % length     (5.185 % on K&M's sigma(+))
    #
    # against a measured 3-4 % near 200 eV.  So velocity sits ABOVE the band by about a fifth.  An earlier reading
    # of 3.868 %, INSIDE the band, was written into the module docstring on 04-Nov-2026 and WITHDRAWN the same day:
    # it was taken at maxKappa = 2, where a 22 % truncation error happened to cancel part of a real overshoot.
    # THE LESSON OUTLASTS THE NUMBER -- A RESULT LANDING INSIDE A MEASURED BAND IS NOT EVIDENCE UNTIL THE
    # CONVERGENCE ERROR IS SMALLER THAN THE BAND.  Compare against the DATA and not against either column of
    # Kornberg & Miraglia, Phys. Rev. A 48, 3714 (1993):  their own C3 model gives 23.9 % here and their paper says
    # the ratio runs "up to five times larger than the data" at intermediate energies, so landing between C2 and
    # the data is the right place to land and agreeing with C3 would have been a failure.
    #
    # ===== EVERY NUMERICAL TRUNCATION IS MEASURED, AND EVERY ONE IS UNDER 1 % ====================================
    #
    #     radial box            0.05 %     over  24 -> 40 a.u.
    #     energy sharings      +0.58 %     over   5 -> 21 Gauss-Legendre points
    #     intermediate mesh    -0.87 %, then +0.24 %   over   9 -> 17 -> 25      non-monotone, hence converged
    #     partial waves        +0.8  %     over  maxKappa 3 -> 4  (the 2 -> 3 step is +21.6 %, so 2 is NOT enough)
    #
    # Against a 20-50 % discrepancy with the measurement.  THE REMAINING ERROR IS THEREFORE IN THE MODEL AND NOT IN
    # THE ARITHMETIC, and the direction is the informative part:  this module is knock-out only, and a model MISSING
    # a mechanism cannot overshoot on physics.  Each total also satisfies its own quadrature identity,
    # sigma = 0.5 SUM_k w_k dsigma/d eps_1, to seven digits -- a unit-blind check, since the weights are in Hartree
    # while the differential cross section is per eV.
    #
    # ===== THE GREEN SPACE IS NOT THE CAUSE (tools/scanPdiHeliumGreenSpace.jl) ===================================
    #
    # The four-level (N-1) Green space was the leading suspect, and the prediction was stated before the run.  It
    # was FALSIFIED.  Opened from 4 levels to 9 (He+ through n = 3) on an orbital span widened to n = 3 and held
    # COMMON to both runs, so that the Green space was the ONLY thing varied:
    #
    #     sigma(2+)        -0.014 % Coulomb,  -0.15 % Babushkin
    #     Babushkin/Coulomb      2.377 -> 2.374
    #     worst single energy-sharing node    -1.07 % Babushkin, at equal sharing
    #
    # against a 0.02 % left-right symmetry noise floor.  Widening the SPAN did move the total, +3.7 % Coulomb and
    # +5.1 % Babushkin, but that is a basis effect and not a Green-space one.  ONE CAVEAT ON THE NULL:  the n = 3
    # orbitals a configuration-averaged mean field generates are over-screened -- Z_eff = 0.67, 0.56, 0.44 for 3s,
    # 3p, 3d against the 2 a He+ state should see -- and all five press against the wall of a 60 a.u. box.  So what
    # is exonerated is the Green space AS A SHARED ORBITAL SET CAN REPRESENT IT.
    #
    # ===== NOR IS THE ORBITAL SET THE LEVER (tools/scanPdiHeliumGaugeVsOrbitals.jl) ==============================
    #
    # The remaining candidate was the shared orbital set:  one set must serve a bound 1s and two continuum
    # electrons, so perhaps the second-order gauge gap merely AMPLIFIES a first-order inconsistency.  Five sets
    # were built, spanning a factor 4.2 in first-order gauge quality, and each run to convergence at maxKappa = 3:
    #
    #     set          mean field / span        eps(1s) [Ha]   1st-order gauge   2nd-order gap   Bab/Cou   sigma(2+)
    #     V^(N-1)      AL on 1s,2s,2p,n=3         -2.0001          4.19 %          138.56 %       2.386    4582 barn
    #     AL n=3       AL on 1s^2 + n=2,n=3       -1.7144          6.99 %          121.29 %       2.213    3933
    #     DFS n=3      DFS on 1s^2 + n=2,n=3      -1.2103         13.13 %          137.70 %       2.377    2902
    #     DFS n=2      DFS on 1s^2 + n=2          -1.1008         13.33 %          132.70 %       2.346    2798
    #     HS n=2       HS  on 1s^2 + n=2          -1.0033         17.47 %          147.85 %       2.479    2268
    #
    # THE FIRST-ORDER GAUGE DEVIATION FALLS BY A FACTOR 4.2 AND THE SECOND-ORDER GAP DOES NOT FOLLOW IT.  The gap
    # spans 121-148 %, a factor 1.11, and it is not even monotone:  the BEST first-order set gives 138.56 % while
    # DFS n=2, three times worse at first order, gives 132.70 %.  A straight-line fit through the five points has
    # r = 0.551 -- on five points indistinguishable from noise -- and a slope of 0.999 points of second-order gap
    # per point of first-order.  TAKE THAT TREND ENTIRELY AT FACE VALUE, which is the most generous possible reading
    # for the orbital hypothesis, and extrapolate to a hypothetical set with ZERO first-order deviation:  the
    # second-order gap is still 124.6 %.  The conclusion is therefore robust to whether the trend is believed at
    # all, and it is that THE ORBITAL AXIS CANNOT CLOSE THIS GAP.
    #
    # WHAT THE ORBITAL SET DOES CONTROL is the MAGNITUDE, cleanly and monotonically:  sigma(2+) runs 2268 -> 4582
    # barn, a factor 2.02, strictly ordered by eps(1s).  So the two symptoms separate -- whatever is wrong
    # multiplies both gauges by an orbital-independent ~2.4 while the absolute scale rides on the 1s binding.  That
    # is a sharper constraint on the defect than the overshoot alone.
    #
    # AND THE TWO ANCHORS CANNOT BE SATISFIED TOGETHER, which is worth knowing before anyone retunes the orbitals.
    # Gauge consistency of <eps|r|1s> optimizes at the He+ limit, eps(1s) = -2.0;  the observable sigma(+) optimizes
    # near eps(1s) = -1.0, where HS reproduces K&M to -0.0 % and gives the closest ratio of the five, 4.172 %;  and
    # the exact He ground state, AL on 1s^2 alone at eps(1s) = -0.917942, is at the wrong end of BOTH.  No single
    # shared set satisfies both.  The proper treatment is TWO sets with a biorthogonal transformation, which
    # conflicts head-on with the one-shared-orbital-set precondition this construction rests on -- see the WHY ONE
    # SET paragraph of branch (c), where the dropped <r|s> term is measured at 2.05x the matrix element on
    # independent sets and 0.294x on a shared one.
    #
    # ===== WHAT IS LEFT, AND WHY THE ITEM CLOSES WITH THE QUESTION OPEN =========================================
    #
    # By elimination the suspect is THE SECOND-ORDER CONSTRUCTION ITSELF -- the assembly of the amplitude from the
    # Green-function sum, not the spaces fed to it.  The first-order and second-order gauge gaps are 13.13 % and
    # 139 % on the production set;  two factors of first-order quality compound to about 28 %, so the construction
    # contributes roughly 5x of its own, and that factor is what remains unexplained.  It was not pursued further
    # because nothing is waiting on it:  the item is parked with the condition that it returns if an explicit
    # application needs a helium double-ionization cross section, on the same footing as the items closed that way
    # on 23-Oct-2026.  Until then this branch is the starting point, and it starts from a clear situation rather
    # than from a disagreement of unknown origin.
    name        = "He, one orbital set, 200 eV"
    nmHe        = Nuclear.Model(2.0)
    gridHe      = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 40.0)
    spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")]
    wa          = Representation(name, nmHe, gridHe, spanConfigs, MeanFieldMultiplet(MeanFieldSettings()) )
    wb          = generate(wa, output=true)
    orbitals    = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
    #
    ## the (N-1)-electron ionic cores on that SAME set, exactly as in branch (c);  four levels, which the 4 -> 9
    ## null above licenses and which is what makes a run of this size affordable
    wg          = Representation(name * " -- Green space", nmHe, gridHe,
                                 [Configuration("1s"), Configuration("2s"), Configuration("2p")],
                                 MeanFieldMultiplet(MeanFieldSettings()); orbitalSource = Basics.GivenOrbitals(orbitals) )
    gMultiplet  = generate(wg, output=true)[ResultKeys.MeanFieldMultiplet]
    #
    doubleSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                         multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[200.],
                         NoEnergySharings=3, NoIntermediateEnergies=9, maxKappa=1, calcDifferentialCs=true,
                         printBefore=true, lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                         eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )

    wc = Atomic.Computation(Atomic.Computation(), name=name, grid=gridHe, nuclearModel=nmHe,
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s^0")],
                            orbitalSource   = Basics.GivenOrbitals(orbitals),
                            processSettings = doubleSettings )

    wd = perform(wc)
    #
end
#
setDefaults("print summary: close", "")


