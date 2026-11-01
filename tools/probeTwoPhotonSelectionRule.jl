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
#  AND THE TWO CHEAP EXPLANATIONS WERE TESTED AND BOTH FAIL, 01-Nov-2026.  The question was whether the gauge
#  gap is ordinary incompleteness, curable by the PSEUDO-CONTINUUM route that `examples/example-Dh.jl` branch b
#  uses (nmax far up, gridStopper = false, box held fixed), which takes H 2s -> 1s from 9.7702 to 8.22902 /s
#  against an exact 8.22906 and closes ITS gauges from 6:1 to 4e-5.  Measured here at Z = 18, Bab/Cou of the
#  mixed midpoint, with the box and the mesh scaled together so the spline density per unit length is constant:
#
#      box [a.u.]   nmax 10   nmax 20   nmax 30     pure midpoint      negative-energy levels kept
#         21.5       1.7564    1.8128    1.9030     0.00e+00 all                  0
#         64.6       1.7564    1.7813    1.7908     0.00e+00 all                  0
#        129.2       1.7564    1.7813    1.7863     0.00e+00 all                  0
#
#  THE BOX WAS PARTLY GUILTY AND IS NOW EXONERATED.  At 21.5 a.u. the ratio marched AWAY from 1 (to 1.9030),
#  which is the signature CLAUDE.md Rule 12 documents for a too-small box -- and enlarging the box removed that
#  march, exactly as the rule says it would.  But 64.6 and 129.2 a.u. then agree to three digits, so the box is
#  converged;  and nmax 20 -> 30 moves the ratio by 0.3 %, so the BOUND pseudo-spectrum is converged as well.
#  The ratio converges to ~1.79 and NOT to 1.  So neither the box nor the basis size is the cause, and
#  "push nmax further" is not the answer.
#
#  WHAT THE ASYMMETRY SAYS, and it is the useful part of this negative.  Babushkin is rock stable -- 1.4142e-16
#  across all nine cells, varying in the fourth digit -- while COULOMB alone drifts (8.06e-17 -> 7.92e-17) and
#  alone fails to meet the other gauge.  The velocity gauge is the one sensitive to completeness of the
#  intermediate spectrum, its negative-energy part included, and the scan reports ZERO negative-energy levels
#  retained at every nmax and every box.  So the Dirac-sea contribution the literature calls mandatory is absent
#  BY CONSTRUCTION here rather than by truncation.  One suggestive number, NOT established: 1.79 is close to the
#  factor ~2 by which a positive-spectrum-only sum is said to overestimate these channels -- but that factor is
#  quoted for a RATE in H-like ions, not for a gauge ratio, so it is a hint and not a measurement.
#
#  SO THE CONDITION FOR S22'S LAST HALF IS NAMED, AND NARROWED BY MEASUREMENT RATHER THAN GUESSED: it needs the
#  NEGATIVE-ENERGY intermediates.  It does NOT need a bigger box, more bound configurations or a better
#  prefactor -- all three are now ruled out with numbers.  The quantitative target to close is the mixed
#  midpoint's Bab/Cou = 1.79 at Z = 18.  The chi route (`DalgarnoLewis.twoPhotonAmplitude`) remains the natural
#  candidate, since it takes its two multipoles independently and is complete by construction, but it is
#  ONE-ELECTRON and would have to be extended to two;  and whether its completeness reaches the negative-energy
#  branch has NOT been checked.  Until then JAC can demonstrate that its mixed-multipole machinery obeys the
#  right selection rule -- which this file does, exactly -- and cannot quote a mixed-multipole rate.
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
