#
#  [S22]  Does the two-photon emission prefactor and its CHANNEL GROUPING generalise to K != 0 and to
#         mixed multipoles?  Challenge S22.  THIS PROBE CHANGES NOTHING IN `src/`: it re-forms the very
#         same `wk` that `MultiPhotonTransition.computeEnergyDiffCs` forms, out of the reduced channels
#         that `Line_2pEmission` already carries, and then asks two questions the production code cannot
#         answer about itself.
#
#  WHAT IT FOUND, 29-Oct-2026 -- and the second half outgrew the question it was written for.
#
#  RUN A -- [S22-ODDK].  ONE QUESTION: does an ODD K vanish at EVERY energy sharing, or only at the
#  symmetric point omega1 = omega2?  The challenge's own record states the strong form ("odd K vanishes by
#  Bose exchange symmetry whenever BOTH photons carry the SAME multipole type, for ANY transition"), and
#  the suspicion was that this had been over-generalised from the midpoint.  MEASURED, H-like 3d_3/2 -> 1s
#  in pure E1E1 over ten sharings:  K = 1 sits at 1e-21 ... 1e-15 against a K = 2 of order 1e-12 at EVERY
#  sharing, not merely at the midpoint.  SO THE SUSPICION WAS WRONG AND THE RECORD STANDS;  the practical
#  conclusion it draws -- that a pure-E1E1 case cannot exercise odd K at all -- is confirmed.
#
#  RUN B -- [S22-CHAN].  ONE QUESTION: of the three M2-involving groups that came back EXACTLY zero in the
#  documented dead end, which are zero because the PHYSICS is zero and which because the code never finds
#  them?  ANSWER: ALL THREE WERE CODE, in three different ways, and fixing them is what this run now
#  demonstrates rather than diagnoses.  Before (H 2s -> 1s, multipoles [E1, M2], sharing 1):
#
#      group        stored     max|amplitude|   retrieved by the production lookup
#      (E1,E1)         4        1.20e-06        found
#      (E1,M2)         2        0.00e+00        found, and exactly zero
#      (M2,E1)         2        0.00e+00        found, and exactly zero
#      (M2,M2)         8        1.23e-17        NOT FOUND -- warning, returns 0
#
#  The three causes, none of which is a selection rule:
#    (1) `InteractionStrength.MabEmission` branched on the GAUGE ALONE and never looked at the multipole, so
#        the M2 factor of a mixed pair, filed under an electric gauge, was computed with the ELECTRIC E2
#        formula.  E2 and M2 carry OPPOSITE parity, so the substituted operator was forbidden for exactly
#        the orbital pairs the right one allows -- hence exactly zero rather than merely wrong.  Proof:
#        M2/Coulomb equalled E2/Coulomb to all 17 digits, and M1/Magnetic = +3.8150591874e-12 where
#        M1/Coulomb = 0, which calibrates the detector.
#    (2) `determineSharingsAndChannels` filed a pure-magnetic pair under the `Magnetic` gauge key, which
#        `getReducedAmplitudeEmission` never asks for -- its one-line `Magnetic` disjunct had been commented
#        out rather than finished.  Built, filled, then reported missing.
#    (3) The same pair was filed once per requested gauge, i.e. stored twice over (8 channels against 4).
#
#  AFTER (same case, same sharing) -- every group live, correctly keyed, and no duplicates:
#
#      gauge key    mp1  mp2   channels   max|amplitude|
#      Babushkin    E1   E1         4     3.037905e-06
#      Babushkin    E1   M2         2     1.648146e-11
#      Babushkin    M2   E1         2     2.093424e-12
#      Babushkin    M2   M2         4     1.759480e-17
#      Coulomb      E1   E1         4     9.516940e-07
#      Coulomb      E1   M2         2     3.798879e-12
#      Coulomb      M2   E1         2     2.298035e-12
#      Coulomb      M2   M2         4     1.759480e-17
#
#  READ THOSE MAGNITUDES AS A STRUCTURAL RESULT AND NOT AS PHYSICS.  This case's intermediate basis is np
#  only, so it does not span the j = 5/2 symmetries an M2 factor can reach -- the run says so itself, with a
#  "the intermediate basis does not span this symmetry" warning per missing path.  What is demonstrated is
#  that the channels EXIST, are keyed where the lookup reaches them, and are no longer duplicated;  a
#  mixed-multipole RATE needs nd/nf intermediates and a system where E1M2 is not a 1e-5 relativistic
#  afterthought.  The literature anchor for that is Derevianko & Johnson, PRA 56, 1288 (1997): He-like
#  2^3S_1 -> 1^1S_0, where E1M2 is 18 % of the 2E1 rate at Z = 30.
#
#  SO S22'S OWN QUESTION -- whether the prefactor's grouping is coherent or incoherent ACROSS (mp1, mp2) --
#  BECAME TESTABLE, because until that was fixed JAC computed no mixed-multipole two-photon amplitude at all.
#
#  RUN C -- [S22-COH], 01-Nov-2026.  THE QUESTION IS ANSWERED, AND A FOURTH DEFECT WAS FOUND ON THE WAY.
#
#  THE DEFECT FIRST, because it had to be repaired before any coherence statement could be measured.
#  `computeEnergyDiffCs` formed BOTH time orderings over ONE intermediate-symmetry list.
#  `AngularMomentum.allowedTotalSymmetries(symf, mpB, mpA, symi)` describes the single path
#  i -> (mpA) -> nu -> (mpB) -> f:  mpA couples to the INITIAL symmetry, mpB to the FINAL one.  So ordering 1
#  needs (symf, mp2, mp1, symi) and ordering 2 needs (symf, mp1, mp2, symi) -- THE SAME LIST ONLY WHEN
#  mp1 == mp2.  Measured, H-like 3d_3/2 -> 1s_1/2 with [E1, M2]:
#
#      mp1  mp2    ordering 1 (symf,mp2,mp1,symi)    ordering 2 (symf,mp1,mp2,symi)
#      E1   E1     1/2-  3/2-                        1/2-  3/2-        <- identical, hence nothing shipped moves
#      E1   M2     3/2-  5/2-                        1/2-  3/2-
#      M2   E1     1/2-  3/2-                        3/2-  5/2-
#      M2   M2     3/2-  5/2-                        3/2-  5/2-        <- identical
#
#  For a MIXED pair the single list therefore summed 5/2- in both orderings, never summed 1/2- in either, and
#  asked the stored channels for combinations that were never generated.  Through `perform` that was 80
#  "No stored U ... returning zero" warnings on one five-sharing run -- a zero that reads like a selection rule,
#  which is how the three defects of 29-Oct read too.  It is the same mistake on a third index.
#
#  AND THE CHECK IS A PROVABLE IDENTITY, not a reference number.  The two emitted photons are identical bosons,
#  so at every sharing x
#
#      wk(mp2,mp1)[x]  ==  (-1)^(L1+L2-K) * wk(mp1,mp2)[1-x]
#
#  and at x = 1/2 the ratio wk(M2,E1)/wk(E1,M2) must be exactly (-1)^(L1+L2-K) = -+1.  Gauss-Legendre nodes are
#  symmetric about the midpoint and an ODD node count puts a node ON it, which is why run C uses noSharings = 5.
#
#      rule              identity, K=1   identity, K=2   midpoint ratio K=1   K=2      channels missing
#      one list          3.2e-01         1.1e-01         Cou +3.106  Bab +2.319   -0.298 / -0.560     80
#      per ordering      0.0e+00         6.4e-17         Cou +1.000000000  -1.000000000               0
#
#  Both gauges independently, and the per-ordering midpoint ratios are +1.000000000 and -1.000000000 exactly.
#  THE FIX IS IN `src/module-MultiPhotonTransition-inc-2p-emission.jl` (01-Nov-2026): the Jsym loop is split in
#  two, one per ordering.  No channel is added -- ordering 2 looks a pair up with the multipole labels swapped,
#  so it lands on the COMPLEMENTARY pair's stored channels, which already carry the symmetries it needs.  The
#  sibling `-inc-2p-bichromatic.jl` had already recognised the same thing and solved it by the UNION of the two
#  lists, which it must, because it stores its channels under that union.
#
#  AND `-inc-2p-absorption.jl` IS NOT THE SAME DEFECT, checked rather than assumed.  It uses one list at five
#  sites and is consistent in doing so: it forms only ONE ordering -- `getReducedAmplitudeAbsorption` is called
#  with the single argument order (mp2, Jsym, omega, mp1) everywhere, and no exchange phase appears in the file
#  at all -- so the list it uses is the right list for the path it builds.  What is ABSENT there is the second
#  ordering itself.  For two monochromatic photons of a mixed pair that is a real and separate question, since
#  the ordered-pair loop supplies both label assignments but with no relative phase between them.  NOT touched.
#
#  NOW S22'S OWN QUESTION.  THE (mp1, mp2) PAIRS ADD INCOHERENTLY, and for the same reason the K do rather than
#  for the reason Jsym does not:  the multipolarity of each photon is a label of the FINAL state, the multipole
#  fields are an orthonormal basis of it, and a polarization-summed angle-integrated rate is a trace over that
#  basis -- so two different pairs are two different final states and cannot interfere.  Jsym is the opposite
#  kind of index: nothing observes the virtual intermediate.  So S22's worry that this was "the same shape of
#  bug as Jsym" is NOT borne out, and the distinction is the one the code's own K note already draws.
#  Two further confirmations, neither of them an appeal to the first:
#    * wk ALREADY CONTAINS BOTH ORDERINGS with the exchange phase, i.e. it IS the Bose-symmetric amplitude for
#      one final state;  adding wk(mp2,mp1) to it would symmetrise a second time.  The measured signature of
#      that double application is the midpoint identity above -- the two groups are mirror images of each other,
#      not independent degrees of freedom -- and the 1/2 at `computeProperties_2pEmission` is the matching
#      statistical factor for running the ordered pairs over the full sharing interval.
#    * the literature's own bookkeeping.  E1M2 is reported as a CONTRIBUTION whose RATIO to the 2E1 one is about
#      0.27 for He-like uranium;  a coherent cross term with 2E1 would instead be of order 2*sqrt(0.27) = 1.04,
#      i.e. larger than the thing it is quoted beside, which no paper reports.  THIS IS SUPPORTING AND NOT
#      DECISIVE:  whether those two numbers are added or interfered is not stated in what was read, and the
#      primary source was NOT checked.  The trace argument does not depend on it.
#
#  AND THE ANCHOR FOR S22'S REMAINING HALF IS A DIFFERENT PAPER FROM THE ONE RECORDED, found 01-Nov-2026.
#  Derevianko & Johnson, PRA 56, 1288 (1997) is the 2E1 BASELINE -- it tabulates two-photon rates for Z = 2..100
#  and has no E1M2 breakdown.  The E1M2 numbers are in PRA 69, 062502, which computes the 2^3S_1 rate with both
#  the 2E1 and the E1M2 amplitudes and puts their ratio near 0.27 for He-like uranium.  That paper also states
#  the structure measured in run C from the other side: THE SELECTION RULE FORBIDDING J = 1 -> 0 TWO-PHOTON
#  DECAY AT EQUAL PHOTON ENERGIES HOLDS ONLY WHEN BOTH PHOTONS CARRY THE SAME MULTIPOLARITY, and not for a mixed
#  amplitude such as E1M2 -- which is the exchange phase (-1)^(L1+L2-K) read as a physical statement, and it is
#  why E1M2 matters most at equal sharing.  Note for whoever takes that half up:  Savukov & Johnson
#  (arXiv:physics/0209086) record that their 2E1 expression differs from Derevianko & Johnson's by alpha^2/4, so
#  the two cannot be compared without settling that factor first.
#  MEASURED SIZE OF THE QUESTION on the case above: summing the four pairs coherently instead changes the
#  integrated rate by under 1e-6 in either gauge, so for THIS system the choice is academic.  It is not
#  academic at Z = 30, where E1M2 is 18 % of 2E1 -- and that case is a correlated two-electron relativistic
#  calculation, which is the part of S22 that remains and is EFFORT M.
#
#  RUN:  julia --project=. tools/probeTwoPhotonMultipoleCoherence.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")

## `wk` EXACTLY AS THE PRODUCTION CODE FORMS IT SINCE 01-Nov-2026: the Jsym sum is coherent and inside the
## modulus, the two photon orderings are added with the exchange phase (-1)^(L1+L2-K), and EACH ORDERING TAKES
## ITS OWN intermediate-symmetry list.  Reproduced here rather than called, because `computeEnergyDiffCs`
## returns only the summed cross section.
function wkOf(line, sharing, K, mp1, mp2, gauge)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)
    wk   = Complex(0.)
    for  Jsym  in  AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
        wk = wk + MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, sharing.omega2, mp2, Jsym,
                                        sharing.omega1, mp1, line.initialLevel, gauge, line.sharings)
    end
    for  Jsym  in  AngularMomentum.allowedTotalSymmetries(symf, mp1, mp2, symi)
        wk = wk + (-1.0)^( mp1.L + mp2.L - Basics.twice(K)/2 ) *
                  MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, sharing.omega1, mp1, Jsym,
                                        sharing.omega2, mp2, line.initialLevel, gauge, line.sharings)
    end
    return( wk )
end


## `wk` AS THE PRODUCTION CODE FORMED IT UNTIL 01-Nov-2026 -- ONE symmetry list for BOTH orderings.  Kept so
## that run C can show the defect and its repair side by side rather than only assert the repair;  it is NOT
## what `src/` does any more, and for mp1 == mp2 it is identical to `wkOf` above.
function wkSingleList(line, sharing, K, mp1, mp2, gauge)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)
    wk   = Complex(0.)
    for  Jsym  in  AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
        wk = wk + MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, sharing.omega2, mp2, Jsym,
                                        sharing.omega1, mp1, line.initialLevel, gauge, line.sharings) +
                  (-1.0)^( mp1.L + mp2.L - Basics.twice(K)/2 ) *
                  MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, sharing.omega1, mp1, Jsym,
                                        sharing.omega2, mp2, line.initialLevel, gauge, line.sharings)
    end
    return( wk )
end

function runCase(tag, nm, grid, initConfs, finalConfs, interConfs, mps, noSharings, scf)
    asf      = AsfSettings(AsfSettings(); scField=scf, gridStopper=false)
    interRep = Representation("intermediate levels", nm, grid, interConfs, MeanFieldMultiplet(MeanFieldSettings(scf, false)))
    interMp  = redirect_stdout(devnull) do
        generate(interRep, output=true)[ResultKeys.MeanFieldMultiplet]
    end
    mpSet = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                scheme = MultiPhotonTransition.TwoPhotonEmissionScheme(
                             MultiPhotonTransition.AbstractMultiPhotonProperty[
                                 MultiPhotonTransition.EnergyDiffCs()], noSharings ),
                multipoles = mps, gauges = [UseCoulomb, UseBabushkin],
                intermediateStates = interMp, calcOverview = false, lineSelection = LineSelection() )
    wa = redirect_stdout(devnull) do
        perform( Atomic.Computation(Atomic.Computation(), name=tag, grid=grid, nuclearModel=nm,
                                   initialConfigs=initConfs, finalConfigs=finalConfs,
                                   initialAsfSettings=asf, finalAsfSettings=asf,
                                   processSettings=mpSet); output=true )
    end
    return( wa[MultiPhotonTransition.Settings] )
end

############################################################  RUN A -- [S22-ODDK]
println("\n", "="^118, "\n  RUN A  [S22-ODDK]   H-like 3d_3/2 -> 1s_1/2, PURE E1E1, resolved in the energy sharing.")
println("  ONE QUESTION: is odd K zero everywhere, or only at omega1 = omega2?\n", "="^118)

nmA   = Nuclear.Model(1.0, PointNucleus())
gridA = Basics.recommendedGrid(Dict(Shell("1s") => 1, Shell("3d") => 1, Shell("8p") => 1), 1.0)
setDefaults("standard grid", gridA)
interA = [Configuration("$(n)p") for n = 2:8]
linesA = runCase("S22-oddK", nmA, gridA, [Configuration("3d")], [Configuration("1s")], interA, [E1], 10, Basics.NuclearField())

for  line in linesA
    Basics.twice(line.initialLevel.J) == 3  ||  continue
    @printf("\n  line %s -> %s,  total omega = %.4f eV,  %d sharings\n",
            string(LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)),
            string(LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)),
            (line.initialLevel.energy-line.finalLevel.energy)*27.211386245988, length(line.sharings))
    println("     omega1/omega_tot      |wk(K=1)|^2        |wk(K=2)|^2      odd-K share")
    for  sh  in  line.sharings
        tot = sh.omega1 + sh.omega2
        w1  = abs(wkOf(line, sh, AngularJ64(1), E1, E1, Basics.Coulomb))^2
        w2  = abs(wkOf(line, sh, AngularJ64(2), E1, E1, Basics.Coulomb))^2
        @printf("        %8.5f        %.6e      %.6e     %.3e\n", sh.omega1/tot, w1, w2, w1/(w1+w2))
    end
end

############################################################  RUN B -- [S22-CHAN]
println("\n\n", "="^118, "\n  RUN B  [S22-CHAN]   H 2s -> 1s with multipoles [E1, M2] -- the DOCUMENTED DEAD END, now instrumented.")
println("  ONE QUESTION: is every (mp1, mp2) group STORED, keyed where the production lookup reaches it, and stored")
println("  exactly once?  Compare against the two tables in this file's header.\n", "="^118)

nmB   = Nuclear.Model(1.0, PointNucleus())
gridB = Basics.recommendedGrid(Dict(Shell("1s") => 1, Shell("2s") => 1, Shell("8p") => 1), 1.0)
setDefaults("standard grid", gridB)
interB = vcat([Configuration("$(n)p") for n = 2:8], [Configuration("$(n)d") for n = 3:6],
              [Configuration("$(n)f") for n = 4:6])
linesB = runCase("S22-chan", nmB, gridB, [Configuration("2s")], [Configuration("1s")], interB, [E1, M2], 4,
                 Basics.NuclearField())

for  line in linesB
    sh = line.sharings[1]
    @printf("\n  line, total omega = %.4f eV;  sharing 1 at omega1/omega_tot = %.5f;  %d stored channels\n",
            (line.initialLevel.energy-line.finalLevel.energy)*27.211386245988,
            sh.omega1/(sh.omega1+sh.omega2), length(sh.channels))

    ## WHAT IS ACTUALLY STORED, grouped by the gauge KEY it was filed under
    println("\n    stored channels, by (gauge key, mp1, mp2):   count   max|amplitude|")
    seen = Dict{Tuple{String,String,String},Vector{Float64}}()
    for ch in sh.channels
        k = (string(ch.gauge), string(ch.multipole1), string(ch.multipole2))
        push!( get!(seen, k, Float64[]), abs(ch.amplitude) )
    end
    for k in sort(collect(keys(seen)))
        @printf("      %-12s %-4s %-4s                     %3d     %.6e\n", k[1], k[2], k[3],
                length(seen[k]), maximum(seen[k]))
    end

    ## WHAT computeEnergyDiffCs ACTUALLY RETRIEVES for the same groups
    println("\n    what the production lookup returns, |wk|^2:")
    println("      mp1   mp2      K        Coulomb          Babushkin         Magnetic")
    for mp1 in [E1, M2],  mp2 in [E1, M2]
        for K in [AngularJ64(0), AngularJ64(1)]
            v = map(g -> abs(wkOf(line, sh, K, mp1, mp2, g))^2,
                    [Basics.Coulomb, Basics.Babushkin, Basics.Magnetic])
            @printf("      %-4s  %-4s   %-4s    %.6e    %.6e    %.6e\n",
                    string(mp1), string(mp2), string(K), v[1], v[2], v[3])
        end
    end
end

############################################################  RUN C -- [S22-COH]
println("\n\n", "="^118, "\n  RUN C  [S22-COH]   H-like 3d_3/2 -> 1s_1/2 with [E1, M2] -- Klist = {1,2}, so a MIXED pair reaches an")
println("  even K and the exchange phase (-1)^(L1+L2-K) takes BOTH signs.  noSharings = 5 is odd, so one node sits")
println("  exactly on omega1 = omega2.")
println("  TWO QUESTIONS: (1) does the provable Bose identity hold -- it does not under one symmetry list;  (2) how")
println("  much would summing the (mp1,mp2) pairs coherently change the rate?\n", "="^118)

nmC   = Nuclear.Model(1.0, PointNucleus())
gridC = Basics.recommendedGrid(Dict(Shell("1s") => 1, Shell("3d") => 1, Shell("8p") => 1), 1.0)
setDefaults("standard grid", gridC)
## nf is needed as well as np: an M2 factor reaches j = 5/2, which no np level carries
interC = vcat([Configuration("$(n)p") for n = 2:8], [Configuration("$(n)f") for n = 4:7])
linesC = runCase("S22-coh", nmC, gridC, [Configuration("3d")], [Configuration("1s")], interC, [E1, M2], 5,
                 Basics.NuclearField())

pairsC = [(E1,E1), (E1,M2), (M2,E1), (M2,M2)]
for  line in linesC
    Basics.twice(line.initialLevel.J) == 3  ||  continue
    ns   = length(line.sharings)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)
    @printf("\n  %s -> %s,  omega_tot = %.4f eV,  %d sharings,  totalRate Cou %.9e  Bab %.9e\n",
            string(symi), string(symf), (line.initialLevel.energy-line.finalLevel.energy)*27.211386245988, ns,
            line.totalRate.Coulomb, line.totalRate.Babushkin)

    println("\n    THE TWO INTERMEDIATE-SYMMETRY LISTS, per (mp1,mp2) -- identical only for a PURE pair:")
    println("      mp1  mp2   ordering 1: (symf,mp2,mp1,symi)       ordering 2: (symf,mp1,mp2,symi)")
    for (a,b) in pairsC
        s1 = AngularMomentum.allowedTotalSymmetries(symf, b, a, symi)
        s2 = AngularMomentum.allowedTotalSymmetries(symf, a, b, symi)
        @printf("      %-4s %-4s  %-34s %s\n", string(a), string(b),
                join(string.(s1), "  "), join(string.(s2), "  "))
    end

    for  gauge in (Basics.Coulomb, Basics.Babushkin)
        println("\n    ### gauge $gauge")
        for (tag, f) in (("ONE LIST (pre-01-Nov-2026)", wkSingleList), ("PER ORDERING (production)  ", wkOf))
            println("      -- $tag --")
            for K in (AngularJ64(1), AngularJ64(2))
                worst = 0.;   ref = 0.
                for j = 1:ns
                    sj = line.sharings[j];    sm = line.sharings[ns+1-j]
                    ph  = (-1.0)^(1 + 2 - Basics.twice(K)/2)
                    lhs = f(line, sj, K, M2, E1, gauge)
                    rhs = ph * f(line, sm, K, E1, M2, gauge)
                    worst = max(worst, abs(lhs - rhs));    ref = max(ref, abs(rhs))
                end
                @printf("         Bose identity, K = %-4s  worst |lhs-rhs| / max|rhs| = %.3e   (must be 0)\n",
                        string(K), ref > 0 ? worst/ref : NaN)
            end
            mid = line.sharings[div(ns+1,2)]
            for K in (AngularJ64(1), AngularJ64(2))
                a = f(line, mid, K, E1, M2, gauge);    b = f(line, mid, K, M2, E1, gauge)
                @printf("         midpoint ratio, K = %-4s  wk(M2,E1)/wk(E1,M2) = %+.9f   must be %+.0f\n",
                        string(K), abs(a) > 0 ? real(b)/real(a) : NaN, (-1.0)^(3 - Basics.twice(K)/2))
            end
            ti = 0.;   tc = 0.
            for sh in line.sharings
                di = 0.;   dc = 0.
                for K in (AngularJ64(1), AngularJ64(2))
                    sm = Complex(0.)
                    for (a,b) in pairsC
                        w = f(line, sh, K, a, b, gauge);    di = di + abs(w)^2;    sm = sm + w
                    end
                    dc = dc + abs(sm)^2
                end
                ti = ti + sh.weight * di;    tc = tc + sh.weight * dc
            end
            @printf("         pair sum:  incoherent %.9e   coherent %.9e   coherent/incoherent = %.9f\n",
                    ti, tc, tc/ti)
        end
    end
end

println("\n   A nonzero Bose identity, or a midpoint ratio that is not exactly -+1, means the two time orderings are")
println("   not being summed over their own intermediate symmetries.  A coherent/incoherent ratio far from 1 would")
println("   mean S22's grouping question changes a number for this system;  measured, it does not.")
println()
