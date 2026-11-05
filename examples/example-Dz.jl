#
println("Dz) Apply & test the PhotoIonizationPseudostates module: a diagonalize-once CI mixing bound")
println("    autoionizing CSFs with a box-discretized continuum pseudostate ladder -- the alternative to")
println("    Dy's amplitude-sum route, built to finish challenges-list S29 by removing the q-parameter")
println("    extraction (and its phase-convention problem) entirely.")

setDefaults("unit: energy", "eV")
setDefaults("print summary: open", "zzz-PhotoIonizationPseudostates.sum")


# THE THREE BRANCHES, IN RISING ORDER OF DIFFICULTY, exactly mirroring example-Dy.jl's own structure (that
# module's branches a-c use its SAME system, this one starts with a cheaper one for branch a, then moves to
# Dy's own already-validated Fe K-LL system for the real cross-check, then the textbook He case).
#
# A DENSE, EXPLICIT-hp GRID IS REQUIRED THROUGHOUT -- found while building this module (05-Nov-2026,
# tools/photoionization-pseudostates-toy): Radial.Grid(true)'s default (pure exponential, hp = 0) gives a
# B-spline basis so sparse that almost no pseudostates fall inside any reasonable energy window. Every branch
# below uses Radial.Grid(Radial.Grid(false); hp=.., rbox=..) instead.


if  true
    # Last successful:  05-Nov-2026   ## H 1s photoionization at E_photon=0.7 Ha=19.048 eV against the EXACT
    #                                     Stobbe (1930)/Bethe-Salpeter closed form: sigma(exact)=9.039e-02 a0^2;
    #                                     JAC PhotoIonization reproduces it to ratio=1.0000 (0.00% gauge gap);
    #                                     THIS module gives sigma=9.100e-02 a0^2, ratio=1.0068 (0.68% high),
    #                                     both gauges identical (exactly as the exact one-electron case must).
    # Branch a: THE LITERATURE CHECK. Hydrogen is exactly one-electron -- no SCF approximation, no frozen-core
    #   assumption, no correlation anywhere -- so a closed-form answer exists and isolates the NUMERICAL METHOD
    #   alone, with nothing else that could be wrong. The exact cross section (Stobbe 1930; Bethe & Salpeter,
    #   "Quantum Mechanics of One- and Two-Electron Atoms", 1957; confirmed against the well-known threshold
    #   value sigma(0.5 Ha) = 6.31 Mb = 2.252e-01 a0^2, itself independently confirmed via
    #   https://arxiv.org/pdf/1810.08071 citing Stobbe's theory, and against the general arctan/arccot closed
    #   form given in https://arxiv.org/pdf/2008.06287):
    #
    #       sigma_1s(E_photon) = (2^5 pi^2 alpha / (3 E_photon^4)) * exp(-4 arctan(k)/k) / (1 - exp(-2 pi/k)),
    #       k = sqrt(2*(E_photon - 0.5))                                             [atomic units throughout]
    #
    #   REPORT: at E_photon = 0.7 Ha (electron energy 0.2 Ha = 5.442 eV above the exact 13.606 eV threshold),
    #   sigma(exact) = 9.039e-02 a0^2. `JAC.PhotoIonization` -- the EXISTING, already-trusted amplitude route,
    #   run here as an independent anchor rather than assumed -- reproduces this to ratio = 1.0000 (both gauges
    #   agree to 0.00%, as they must for an exact one-electron problem with no competing approximation). THIS
    #   module's diagonalize-once, density-of-states route gives sigma = 9.100e-02 a0^2, ratio = 1.0068 -- 0.7%
    #   high, both gauges IDENTICAL to the printed digits (again as the exact case requires).
    #
    #   THE FIRST VERSION OF THIS CHECK GAVE RATIO = 1.0636 (6.4% high), AND THE DIAGNOSIS IS NOW UNDERSTOOD
    #   AND FIXED, NOT JUST BOUNDED. Printing the individual stick-state sigma_k around the evaluation point
    #   showed them SMOOTH and monotonically falling with energy -- no state-to-state noise to average away --
    #   so the error was never statistical. `crossSectionAt` picked the single NEAREST available pseudostate's
    #   own energy rather than the exact requested one, on a cross section that falls off steeply here; linear
    #   interpolation between the two bracketing sticks (now the module's own `crossSectionAt`, replacing the
    #   earlier nearest-neighbour form) closed virtually all of the gap. A hand-computed interpolation between
    #   the two sticks either side of the target predicted ratio = 1.0068 before the fix was even made --
    #   confirmed by the actual re-run above to the digit.
    #
    #   System: H 1s -> p (bare proton) + e-. The "ion" has ZERO electrons and an EMPTY orbitals dict, so
    #   checkFrozenOrbitalConsistency has nothing to compare and passes trivially -- the one case where this
    #   module's own frozen-core requirement raises no question at all.
    grid   = Radial.Grid( Radial.Grid(false); hp=0.02, rbox=40.0 )
    nModel = Nuclear.Model(1., Nuclear.PointNucleus())

    fMult = SelfConsistent.performSCF([Configuration("1s^0")], nModel, grid, AsfSettings(); printout=false)
    finalIonLevel = fMult.levels[1]
    iMult = SelfConsistent.performSCF([Configuration("1s")], nModel, grid, AsfSettings(); printout=false)
    initialLevel = iMult.levels[1]
    @printf("\n>>> H 1s level energy = %.6f Ha  (exact = -0.5)\n", initialLevel.energy)

    photonEnergyHa = 0.7   ## = 0.5 (exact threshold) + 0.2 Ha electron energy, comfortably inside the window
    photonEv       = Defaults.convertUnits("energy: from atomic to eV", photonEnergyHa)

    function stobbeCrossSection(Egamma::Float64)   ## Stobbe's exact closed form, a0^2; see the REPORT above
        k     = sqrt(2*(Egamma - 0.5))
        alpha = Defaults.getDefaults("alpha")
        return( (2^5 * pi^2 * alpha) / (3 * Egamma^4) * exp(-4*atan(k)/k) / (1 - exp(-2*pi/k)) )
    end
    sigmaExact = stobbeCrossSection(photonEnergyHa)
    @printf(">>> EXACT (Stobbe) sigma = %.6e a0^2\n", sigmaExact)

    piSet   = PhotoIonization.Settings(PhotoIonization.Settings(); multipoles=[E1], gauges=[UseCoulomb, UseBabushkin],
                  photonEnergies=[photonEv], calcAnisotropy=false, printBefore=false)
    piLines = perform( Atomic.Computation(Atomic.Computation(), name="H PI reference", grid=grid, nuclearModel=nModel,
                  initialConfigs=[Configuration("1s")], finalConfigs=[Configuration("1s^0")], processSettings=piSet);
                  output=true )[PhotoIonization.Settings]
    for  line  in  piLines
        @printf(">>> JAC PhotoIonization:        sigma = (C=%.6e, B=%.6e);  ratio to exact (C) = %.4f\n",
                line.crossSection.Coulomb, line.crossSection.Babushkin, line.crossSection.Coulomb/sigmaExact)
    end

    ppSet    = PhotoIonizationPseudostates.Settings(PhotoIonizationPseudostates.Settings(); multipoles=[E1],
                   gauges=[UseCoulomb, UseBabushkin], pseudostateWindow=(0.05, 1.0), printBefore=false)
    spectra  = PhotoIonizationPseudostates.computeSpectra(fMult, initialLevel, nModel, grid, ppSet; output=true)

    sigmaPPTotal = EmProperty(0., 0.)
    for  spectrum  in  spectra
        sigmaPP = PhotoIonizationPseudostates.crossSectionAt(spectrum, photonEnergyHa)
        @printf(">>> kappa = %d:  sigma(Pseudostates) = (C=%.6e, B=%.6e)\n", spectrum.kappa, sigmaPP.Coulomb, sigmaPP.Babushkin)
        global sigmaPPTotal = sigmaPPTotal + sigmaPP
    end
    @printf(">>> PhotoIonizationPseudostates, summed over kappa: sigma = (C=%.6e, B=%.6e);  ratio to exact (C) = %.4f, (B) = %.4f\n",
            sigmaPPTotal.Coulomb, sigmaPPTotal.Babushkin, sigmaPPTotal.Coulomb/sigmaExact, sigmaPPTotal.Babushkin/sigmaExact)
    #
elseif  false
    # Last visit:  05-Nov-2026   ## RAN, BUT UNINFORMATIVE AS WRITTEN -- see the box-size finding below. The
    #                                comparison itself could not yet be made; no number is quoted as a result.
    #                                RE-CHECKED after branch a's interpolation fix to crossSectionAt (same
    #                                day): UNCHANGED -- still 0 pseudostates in the window, since this is a
    #                                box-COVERAGE problem (no states exist there at all), not an interpolation
    #                                one. The two findings are independent; fixing one does not touch the other.
    # Branch b: CROSS-CHECK AGAINST THE ALREADY-VALIDATED Fe K-LL SYSTEM of example-Dy.jl branches a-c.
    #   Two numerically independent techniques (amplitude-sum-plus-Fano-formula, and diagonalize-once-plus-
    #   pseudostates) computing the SAME observable near the SAME resonance is the strongest test available
    #   to either one -- stronger than either module's own internal consistency checks. THAT COMPARISON IS
    #   NOT YET REACHED, for a reason worth recording precisely rather than retrying blindly.
    #
    #   initialLevel's own basis spans BOTH the Be-like ground (1s^2 2s^2) AND the K-hole autoionizing
    #   configuration (1s 2s^2 2p) together, in ONE CI -- exactly the design this module's own docstring
    #   describes, replacing PhotoIonizationInterference's separate intermediateMultiplet.
    #
    #   THE REAL FINDING, 05-Nov-2026: this resonance's electron energy (resEv - threshold) is ~171 Ha --
    #   Fe is highly charged, so even a K-hole-to-continuum transition carries a huge kinetic energy. The
    #   grid here (rbox=4.0, matched to Fe's own tiny BOUND orbitals, exactly as example-Dy.jl's own branches
    #   use it) gives a box whose positive-energy level spacing GROWS with energy (free-particle-like,
    #   E_n ~ (n*pi/rbox)^2/2) -- at 171 Ha that spacing is of order TEN Ha, not a fraction of one. A first
    #   attempt with pseudostateWindow=(0.05,300.0) ran and returned a cross section ~300-600x SMALLER than
    #   example-Dy.jl branch c's own peak (sigma_resonant=1.336798); narrowing the window to (169.0,174.0),
    #   bracketing the resonance tightly, found ZERO pseudostates inside it at all. This is Rule 12's own
    #   lesson (a continuum orbital needs a box matched to ITS energy, not to the bound orbitals) landing on
    #   this module specifically: `buildPseudostateLadder` reuses the SAME grid/primitives the bound SCF
    #   used, with no mechanism to give the pseudostate ladder its OWN, independently-sized box. Doing so
    #   properly is a real design question for a later round (the two bases would then need to be reconciled
    #   at the level this module's "ONE consistent one-electron basis" design currently assumes they already
    #   are) -- NOT attempted here. A cheaper resonance, closer to its own threshold, would sidestep the
    #   problem but none of example-Dy.jl's own validated K-hole resonances qualifies; all four sit near the
    #   same ~171 Ha electron energy.
    grid   = Radial.Grid( Radial.Grid(false); hp=0.02, rbox=4.0 )
    nModel = Nuclear.Model(26.)
    # The INITIAL multiplet is built FIRST, ordinarily (self-consistently), spanning BOTH the ground and the
    # K-hole configuration together -- it needs its own 2p orbital that a bare "1s^2 2s" ion never would.  The
    # ion is then built AS A FROZEN-ORBITAL CI reusing the initial's own 1s/2s, the opposite order from branch
    # a, where the ion came first -- whichever side needs the richer orbital set must be computed first.
    iMult = SelfConsistent.performSCF([Configuration("1s^2 2s^2"), Configuration("1s 2s^2 2p")], nModel, grid,
                                       AsfSettings(); printout=false)
    initialLevel = iMult.levels[1]
    fMult = Hamiltonian.performCIwithFrozenOrbitals([Configuration("1s^2 2s")], initialLevel.basis.orbitals,
                                                     nModel, grid, AsfSettings(); printout=false)
    finalIonLevel = fMult.levels[1]
    resEv  = 6605.767   ## the J=1- K-hole resonance example-Dy.jl already validated
    thrEv  = Defaults.convertUnits("energy: from atomic to eV", finalIonLevel.energy - initialLevel.energy)
    @printf("\n>>> ionization threshold = %.4f eV;  resonance at %.4f eV\n", thrEv, resEv)

    ppSet   = PhotoIonizationPseudostates.Settings(PhotoIonizationPseudostates.Settings(); multipoles=[E1],
                  gauges=[UseCoulomb, UseBabushkin], pseudostateWindow=(0.05, 300.0), printBefore=false)
    spectra = PhotoIonizationPseudostates.computeSpectra(fMult, initialLevel, nModel, grid, ppSet; output=true)

    resHa = Defaults.convertUnits("energy: to atomic", resEv)
    for  spectrum  in  spectra
        sigma = PhotoIonizationPseudostates.crossSectionAt(spectrum, resHa)
        @printf(">>> kappa = %d:  sigma(Pseudostates, at resonance) = (C=%.6e, B=%.6e)\n",
                spectrum.kappa, sigma.Coulomb, sigma.Babushkin)
    end
    println(">>> compare against example-Dy.jl branch c's own peak (sigma_resonant = 1.336798, Coulomb, at")
    println(">>> omega = 6605.767 eV).")
    #
elseif  false
    # Last visit:  05-Nov-2026   ## RAN, but sigma came out PERFECTLY FLAT (8.543818e-05 Coulomb, kappa=-2;
    #                                8.322839e-05, kappa=1) across the whole +-0.3 eV scan -- the SAME
    #                                underlying resolution problem as branch b, in a different guise; see
    #                                below. No lineshape to judge, so no q and no date.
    #                                RE-CHECKED after branch a's interpolation fix, same day: the scan is no
    #                                longer artificially FLAT (interpolation now gives a smooth ramp between
    #                                the same two widely-spaced sticks, e.g. 3.80e-3 -> 5.94e-3 across the
    #                                last ten points of the scan) -- but there are still FAR too few sticks
    #                                across the resonance's own Gamma ~ 38 meV to show a real dip/peak, only a
    #                                straight-line interpolation between two points that happen to bracket it.
    #                                Less misleading than a flat plateau; still not a lineshape. The real
    #                                blocker is unchanged: box coverage, not interpolation.
    # Branch c: THE TEXTBOOK CASE -- He 1s^2 -> 2s2p ^1P deg (Madden & Codling 1963: q ~ -2.75), S29's own
    #   named validation target, left unsolved by Dy's own amplitude-ratio route (branch d there: q = 0.568,
    #   not a reproduction). Reused setup: SAME iConfigs/mConfigs/fConfigs/grid/EOL settings as example-Dy.jl
    #   branch d, so any difference in the result is attributable to the METHOD and not to the system.
    #
    #   initialLevel's basis spans the 1s^2/2p^2 ground CI AND the 2s2p autoionizing configuration together.
    #
    #   THE SAME ROOT CAUSE AS BRANCH B, MEASURED HERE IN ITS OTHER FORM. The resonance itself (Gamma ~ 38 meV
    #   = 0.0014 Ha) is narrower than the LOCAL pseudostate spacing this grid delivers at the resonance's own
    #   electron energy (~31.9 eV = 1.17 Ha) -- so every point of the +-0.3 eV (~0.011 Ha, already ~8 Gamma
    #   wide) scan landed on the SAME nearest stick state via crossSectionAt's nearest-neighbour lookup,
    #   giving a perfectly flat, uninformative result rather than a visible peak. Branch b's problem was a box
    #   too SMALL for a far-above-threshold resonance; this one is a box that is the right general size for
    #   He but still far too COARSE, LOCALLY, to resolve a sub-eV feature -- the same underlying limitation
    #   (the pseudostate ladder shares the bound SCF's own grid, with no way to ask for finer LOCAL resolution
    #   just where a narrow resonance sits) producing two different symptoms. See the module's own docstring,
    #   which now names this as a stated limitation, and challenges-list S29 for where this is held.
    grid   = Radial.Grid( Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=1.0e-2, rbox=25.0 )
    nModel = Nuclear.Model(2.)

    # iMult FIRST, as in branch b -- it needs the richer orbital set (2p, 2s besides 1s); the ion then reuses
    # iMult's OWN 1s orbital frozen, exactly the order checkFrozenOrbitalConsistency requires.
    iSettings = AsfSettings(AsfSettings(); scField=Basics.EOLField(), scfRoute=Basics.RotationRoute(60),
                            levelSelectionCI=LevelSelection(true, symmetries=[LevelSymmetry(AngularJ64(0), Basics.plus),
                                                                                LevelSymmetry(AngularJ64(1), Basics.minus)]))
    iMult = SelfConsistent.performSCF([Configuration("1s^2"), Configuration("2p^2"), Configuration("2s 2p")],
                                      nModel, grid, iSettings; printout=false)
    initialLevel = iMult.levels[1]
    fMult = Hamiltonian.performCIwithFrozenOrbitals([Configuration("1s")], initialLevel.basis.orbitals, nModel,
                                                     grid, AsfSettings(); printout=false)
    finalIonLevel = fMult.levels[1]
    dIdx = findfirst( lev -> LevelSymmetry(lev.J, lev.parity) == LevelSymmetry(AngularJ64(1), Basics.minus), iMult.levels )
    resEv = Defaults.convertUnits("energy: from atomic to eV", iMult.levels[dIdx].energy - initialLevel.energy)
    thrEv = Defaults.convertUnits("energy: from atomic to eV", finalIonLevel.energy - initialLevel.energy)
    @printf("\n>>> resonance (2s2p ^1P) at %.4f eV above ground;  ionization threshold = %.4f eV\n", resEv, thrEv)

    ppSet   = PhotoIonizationPseudostates.Settings(PhotoIonizationPseudostates.Settings(); multipoles=[E1],
                  gauges=[UseCoulomb], pseudostateWindow=(0.01, 2.0), printBefore=false)
    spectra = PhotoIonizationPseudostates.computeSpectra(fMult, initialLevel, nModel, grid, ppSet; output=true)

    for  omega in (resEv-0.3):0.02:(resEv+0.3)
        omegaHa = Defaults.convertUnits("energy: to atomic", omega)
        for  spectrum  in  spectra
            sigma = PhotoIonizationPseudostates.crossSectionAt(spectrum, omegaHa)
            @printf("   omega = %.4f eV,  kappa = %d:  sigma = %.6e\n", omega, spectrum.kappa, sigma.Coulomb)
        end
    end
    println(">>> report the LINESHAPE directly -- no q-parameter extraction, no phase convention to assume --")
    println(">>> against Madden & Codling's q ~ -2.75, Gamma ~ 38 meV.  `# Last successful` deliberately left")
    println(">>> blank: this needs a human physics judgement of the shape, not a mechanical pass/fail.")
    #
end
