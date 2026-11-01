#
#  [S22-RULE]  THE J = 1 -> 0 EQUAL-SHARING SELECTION RULE, and why the mixed-multipole MAGNITUDE is still not
#              quotable from this route.  Challenge S22's remaining half, measured 01-Nov-2026.
#
#  THE OBSERVABLE, and it is background-free, which is the whole reason for choosing this system.  He-like
#  1s2s 3S_1 -> 1s^2 1S_0 has J_i = 1 and J_f = 0, so Klist = {1} and NOTHING ELSE.  At omega1 = omega2 the two
#  photon time orderings are the SAME reduced amplitude, so  wk = U * (1 + (-1)^(L1+L2-K)):
#
#      pure E1E1 :  (-1)^(1+1-1) = -1   ->   wk = 0 EXACTLY, whatever the radial content
#      E1M2      :  (-1)^(1+2-1) = +1   ->   wk = 2U, not suppressed at all
#
#  That is the published rule -- J = 1 -> 0 two-photon decay is forbidden at equal photon energies ONLY when both
#  photons carry the same multipolarity (Phys. Rev. A 69, 062502) -- and it makes the MIDPOINT OF THE SPECTRUM a
#  measurement of the mixed-multipole channel with no pure-pair background to subtract.  It is also the reason
#  E1M2 reaches ~27 % of 2E1 in that system while contributing under 0.5 % to H-like 2s -> 1s at any Z: the
#  denominator is suppressed by the selection rule, not the numerator enhanced.
#
#  WHAT PASSES, AND IT IS THE STRONGEST END-TO-END CHECK THIS MACHINERY HAS HAD.  The rule is satisfied
#  IDENTICALLY -- 0.000000e+00, not merely small -- in BOTH gauges at Z = 2, 18 and 30, and every mirror pair of
#  sharings agrees to all printed digits.  It needs no literature number, no normalisation and none of the open
#  prefactor, and it tests the exchange phase, the per-ordering symmetry lists of `fee41f2f` and the whole K = 1
#  assembly at once.  A nonzero midpoint under pure E1E1 would mean one of those three is wrong.
#
#  WHAT DOES NOT PASS IS THE MAGNITUDE, and the obstruction is STRUCTURAL rather than "more configurations":
#
#      Z    pure E1E1 at x=1/2   with M2 at x=1/2 (Bab)   mixed/peak (Bab)   (Z alpha)^2   pure gauge ratio
#      2        0.000000e+00          1.3795e-25               0.951          2.1e-04            44
#      18       0.000000e+00          1.4158e-16               0.140          1.7e-02            1.9
#      30       0.000000e+00          9.3644e-15               0.095          4.8e-02            1.46
#
#  FOR J = 1 -> 0 THE SELECTION RULE SUPPRESSES THE 2E1 CHANNEL ITSELF, so what survives the cancellation is a
#  small RESIDUAL of a truncated intermediate sum -- a near-cancellation, which is exactly where a truncated sum
#  is worth nothing.  At Z = 2 that shows plainly: mixed/peak comes out 0.951 in Babushkin where (Z alpha)^2 puts
#  it near 2e-04, three orders too large, because the Babushkin PURE channel has been destroyed by truncation
#  (its two gauges sit a factor 44 apart).  At Z = 18 and 30 the mixed channel's own two gauges still disagree by
#  1.76 and 1.73.  This is priority item 66's finding in its sharpest form: the route is measured sound at Z = 1
#  and measured to fail by Z = 18, and the sharpening is that for a J = 1 -> 0 transition it fails at EVERY Z,
#  because the physics itself removes the large term.
#
#  AND A SECOND, INDEPENDENT BLOCKER FROM THE LITERATURE, which no amount of bound-state basis fixes: for mixed
#  two-photon channels the NEGATIVE-ENERGY (Dirac-sea) intermediate states are mandatory.  Restricting the
#  intermediate sum to the positive spectrum overestimates the E1M1 and E1E2 rates by factors of about 2 for
#  neutral hydrogen and about 2.5 for hydrogen-like uranium (arXiv:0908.2498, arXiv:0904.1503).  The intermediate
#  set used here is positive-energy BOUND states only.
#
#  SO THE CONDITION FOR S22'S LAST HALF IS NAMED RATHER THAN GUESSED: it needs a COMPLETE intermediate spectrum
#  including the negative-energy continuum, i.e. the chi route (`DalgarnoLewis.twoPhotonAmplitude`, which already
#  takes its two multipoles independently and is complete by construction) extended to TWO electrons.  It does
#  NOT need more bound configurations, and it does not need a better prefactor.  Until then JAC can demonstrate
#  that its mixed-multipole machinery obeys the right selection rule -- which this file does -- and cannot quote
#  a mixed-multipole rate.
#
#  See `tools/probeTwoPhotonMultipoleCoherence.jl` for the channel-GROUPING half of S22, which is settled.
#
#  RUN:  julia --project=. tools/probeTwoPhotonSelectionRule.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator
setDefaults("unit: energy", "eV")

for  (Z, nmax) in [(2.0, 10), (18.0, 10), (30.0, 10)]
    println("\n", "="^116)
    @printf("  He-like Z = %.0f:  1s2s 3S_1 -> 1s^2 1S_0,  Klist = {1} only\n", Z)
    println("="^116)
    ni     = Nuclear.Model(Z, Z < 2.5 ? PointNucleus() : FermiNucleus())
    g0     = Basics.recommendedGrid(Dict(Shell(1,0)=>1, Shell(nmax,1)=>1), Z; printout=false)
    grid   = Basics.recommendedGrid(Dict(Shell(1,0)=>1, Shell(nmax,1)=>1), Z;
                                    rbox=g0.tL[end], hp=g0.tL[end]/3000, printout=false)
    inter  = [Configuration("1s $(n)p") for n = 2:nmax]
    iRep   = Representation("1snp", ni, grid, inter, MeanFieldMultiplet(MeanFieldSettings()))
    iMp    = redirect_stdout(devnull) do;  generate(iRep, output=true)[ResultKeys.MeanFieldMultiplet];  end

    res = Dict{String,Any}()
    for  mps in ([E1], [E1, M2])
        set = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                  scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                               MultiPhotonTransition.AbstractMultiPhotonProperty[
                                   MultiPhotonTransition.EnergyDiffCs()], 9 ),   ## ODD: a node sits on x = 1/2
                  multipoles = mps, gauges = [UseCoulomb, UseBabushkin],
                  intermediateStates = iMp, calcOverview = false,
                  lineSelection = LineSelection(true,
                      symmetryPairs = [(LevelSymmetry(1, Basics.plus), LevelSymmetry(0, Basics.plus))]) )
        lines = redirect_stdout(devnull) do
            perform( Atomic.Computation(Atomic.Computation(), name="S22mag Z=$Z", grid=grid, nuclearModel=ni,
                                        initialConfigs=[Configuration("1s 2s")],
                                        finalConfigs=[Configuration("1s^2")],
                                        processSettings=set); output=true )[MultiPhotonTransition.Settings]
        end
        res[string(mps)] = lines
    end

    for  key in ("EmMultipole[E1]", "EmMultipole[E1, M2]")
        lines = res[key]
        for  line in lines
            Basics.twice(line.initialLevel.J) == 2  ||  continue
            ns = length(line.sharings);   mid = div(ns+1, 2)
            @printf("\n  %-20s  omega_tot = %10.4f eV   totalRate Cou %.6e  Bab %.6e\n", key,
                    (line.initialLevel.energy-line.finalLevel.energy)*27.211386245988,
                    line.totalRate.Coulomb, line.totalRate.Babushkin)
            peakC = maximum(s.differentialRate.Coulomb   for s in line.sharings)
            peakB = maximum(s.differentialRate.Babushkin for s in line.sharings)
            for  (j, s) in enumerate(line.sharings)
                @printf("      x = %.6f   Cou %.6e   Bab %.6e%s\n", s.omega1/(s.omega1+s.omega2),
                        s.differentialRate.Coulomb, s.differentialRate.Babushkin,
                        j == mid ? "   <== MIDPOINT: must VANISH for a pure pair" : "")
            end
            m = line.sharings[mid]
            @printf("      midpoint / peak :   Coulomb %.3e    Babushkin %.3e\n",
                    m.differentialRate.Coulomb/peakC, m.differentialRate.Babushkin/peakB)
        end
    end

    ## the mixed-multipole channel, read off the midpoint where the pure pair is identically zero
    lE1 = [l for l in res["EmMultipole[E1]"]      if Basics.twice(l.initialLevel.J) == 2]
    lMx = [l for l in res["EmMultipole[E1, M2]"]  if Basics.twice(l.initialLevel.J) == 2]
    if  !isempty(lE1) && !isempty(lMx)
        ns = length(lE1[1].sharings);   mid = div(ns+1, 2)
        a  = lE1[1].sharings[mid];      b = lMx[1].sharings[mid]
        pk = maximum(s.differentialRate.Babushkin for s in lMx[1].sharings)
        @printf("\n  >>> Z = %.0f   MIXED-MULTIPOLE MIDPOINT, background-free:\n", Z)
        @printf("        pure E1E1 at x=1/2      Cou %.6e   Bab %.6e   (the selection rule: must be ~0)\n",
                a.differentialRate.Coulomb, a.differentialRate.Babushkin)
        @printf("        with M2 at x=1/2        Cou %.6e   Bab %.6e\n",
                b.differentialRate.Coulomb, b.differentialRate.Babushkin)
        @printf("        mixed / spectrum peak   Bab %.4e      (Z*alpha)^2 = %.4e\n",
                pk > 0 ? b.differentialRate.Babushkin/pk : NaN, (Z/137.035999)^2)
    end
end
println()
