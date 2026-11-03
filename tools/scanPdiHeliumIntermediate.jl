#
#  IS THE INTERMEDIATE-ENERGY INTEGRAL WHAT MAKES HELIUM'S DOUBLE-IONIZATION CROSS SECTION MOVE?
#  03-Nov-2026.
#
#  WHAT THE BOX SCAN OF 02/03-Nov SETTLED, and it refuted its own hypothesis.  Five boxes from 12 to 40 a.u.
#  at 200 eV, with the Green space built two ways:
#
#    * the GREEN SPACE's ORBITALS ARE NOT THE CAUSE.  An SCF on the He+ configurations gives 2s and 2p at
#      -0.50003 Ha at EVERY ONE of the five boxes, and even the diffuse neutral-span orbitals settle by 24 a.u.
#      (-0.06506, -0.07163, -0.07219, -0.07223, -0.07223).  Both bases are box-converged; the cross section is
#      not.  The prediction written into tools/scanPdiHeliumBox.jl -- ion arm flat, neutral arm drifting -- is
#      therefore wrong, and that file says what to look at next: the number of intermediate states.
#    * THE B-SPLINE BASIS IS NOT STARVED either.  It GROWS with the box, 217 -> 620 functions per symmetry over
#      12 -> 40 a.u., so Rule 12's too-large failure does not apply here.
#    * THE ION ARM IS ALREADY NEARLY CONVERGED over 12 to 32 a.u.:  sigma2+ (velocity) = 3790, 3613, 3831,
#      3592 barn, a spread of 3.4 % about 3706.  Box 40 then falls to 2457 barn, 34 % below the plateau.
#    * AND THE ION ARM OBEYS THE ONE EXACT SYMMETRY of the problem.  Two identical electrons leaving a J=0
#      state make dsigma/d eps_1 symmetric under eps_1 <-> eps_2; the ion arm satisfies that to 1e-5 or better
#      at every box, while the NEUTRAL arm breaks it by 37 % at 32 a.u. and 4 % at 40 a.u.  The symmetry is a
#      free, exact correctness test, and it disqualifies the neutral-span Green space outright.
#
#  SO THE SUSPECT IS THE INTERMEDIATE-ENERGY INTEGRAL, and the box scan never varied it.  The whole scan ran at
#  the defaults: `NoIntermediateEnergies = 0`, which falls back to `NoEnergySharings` = 9, and
#  `intermediateEnergyFactor = 3`.  Nine quadrature points for a PRINCIPAL-VALUE integral with a pole inside the
#  interval is crude, and its error is carried by free-free matrix elements between continuum orbitals -- which
#  are the one class of quantity in this calculation that a wall at finite radius does change.  That is also why
#  single ionization is flat to 0.3 % across the same five boxes: it has no free-free integral in it at all.
#
#  THE HYPOTHESIS, stated before the run:  the box dependence of sigma2+ is a SYMPTOM of an unconverged
#  intermediate-energy integral, not of the box.  If so, raising NoIntermediateEnergies must (a) make sigma2+
#  settle at a fixed box and (b) bring boxes 24 and 40 into agreement, and the 34 % collapse at 40 a.u. will
#  turn out to be quadrature noise rather than physics.
#
#  IF INSTEAD sigma2+ is already converged in this integral at 9 points, the hypothesis is dead, boxes 24 and 40
#  will still disagree, and the next suspect is the partial-wave set maxKappa, which the box scan also held fixed.
#
#  ------------------------------------------------------------------------------------------------------------------
#  THE ANSWER, 03-Nov-2026:  THE HYPOTHESIS IS DEAD, AND THE WAY THE RUN DIED IS WHY.
#
#  Three rungs at box 24 a.u. and factor 3, in barn:
#
#      NoIntermediateEnergies      9          17          25
#      Coulomb                  4573.0      4520.7      4547.6        spread 1.2 %
#      Babushkin               14153.2     14150.3     14393.8        spread 1.7 %
#
#  So nine points already carry this integral to about one percent, and it CANNOT be what moves sigma2+ by 34 %
#  between two boxes.  The fourth rung then RAISED, and the message names the real defect:
#
#      energy = 3.5065e-02 Ha = 0.954 eV, kappa = 2, box 24 a.u.
#      "Pure-sine normalized continuum orbital with normalization constant N=3.3392e+01"
#      ERROR: Pure-sine normalized N=3.3392e+01==> enlarge box-size
#
#  A FINER MESH PUTS A QUADRATURE POINT AT A LOWER ENERGY, AND THE CONTINUUM MACHINERY REFUSES IT.  The advice in
#  that message is misdirected.  `normalizeOrbitalPureSine` fixes an orbital's scale from `max|P|` over the LAST
#  300 GRID POINTS -- 2.98 a.u. at hp = 1e-2 -- which is a fixed POINT COUNT, not a fixed fraction of an
#  oscillation.  At 0.954 eV the half-wavelength is 11.9 a.u., so the window spans 0.125 of a half-wave and need
#  contain no crest at all; the maximum it finds is then the wave's value at an arbitrary phase, N blows up, and
#  the guard fires.  ENLARGING THE BOX DOES NOT FIX THAT -- it only re-rolls where the window lands.  Measured on
#  the final-state field at 1.95 eV, the true amplitude divided by the one the window assumed is 1.000, 1.000,
#  1.000, 1.118, 1.006 at boxes 12, 18, 24, 32, 40:  right by luck at four boxes and 11.8 % high at the fifth.
#
#  AND THE SAME FILE ALREADY KNOWS THIS.  `Continuum.checkContinuumAmplitude`, written 22-Oct-2026, says a window
#  shorter than half an oscillation "samples the outgoing wave at one arbitrary phase" and sizes ITS window by the
#  wavelength for exactly that reason.  The lesson was applied to the check and never to the normalization it
#  checks.  That is the proposed repair, and the falsifiable prediction attached to it is that the box dependence
#  of sigma2+ largely collapses once the window is sized by the wavelength.
#
#  Two further negatives worth keeping, so they are not retried.  `CoulombSineNorm` is NOT the lever:  in JAC the
#  sine and Coulomb conventions differ only in the PHASE they assign, and give amplitude scales equal to four
#  decimals at every energy from 1.95 to 175 eV.  And initial-state correlation cannot be added by widening the
#  configuration list:  on one shared orbital set the averaged mean field RAISES the 1s^2 energy by 0.0425 Ha
#  while a 26-CSF treatment wins back only 0.0185 Ha, a net 0.0240 Ha WORSE than the single configuration, whose
#  own -2.85797 Ha against the exact -2.90372 Ha carries none of the correlation energy at all.
#
#  WHAT THIS RUN NOW DOES, with the rungs that died re-cast so a refusal cannot kill the rest:  the upper limit of
#  the integral at the largest mesh that survives, and the box-40 cross-check that the ladder never reached.
#
#  ------------------------------------------------------------------------------------------------------------------
#  AND THE RESULT, 03-Nov-2026.  Seven runs at box 24 a.u., velocity gauge, 3 sharings, in barn:
#
#      points        9      17      25      25      25      33
#      upper limit  3x      3x      3x      5x      8x      8x
#                 4573.0  4520.7  4547.6  4530.4  4489.1  4545.6      total spread 1.9 %
#
#  The sum over intermediate states is CONVERGED across a 3.7-fold range of points and a 2.7-fold range of upper
#  limit.  The same settings at box 40 a.u. give 2617.4 barn against 4573.0 -- 42.8 % lower -- with the Green
#  space identical at -0.50003 Ha in both.  So the box dependence survives everything, and by elimination it is
#  in the CONTINUUM ORBITALS.
#
#  THE CONTROLLED FALSIFICATION OF THE GUARD'S OWN ADVICE.  The same electron, energy 5.3244e-02 Ha = 1.449 eV,
#  kappa = -2, identical in every other respect:
#
#      box 24 a.u.    N = 1.7549e+01    accepted
#      box 40 a.u.    N = 3.1536e+01    REFUSED -- "enlarge box-size"
#
#  The sphere grew by 1.67 and the constant by 1.80, crossing the guard's fixed limit of 30.  ENLARGING THE BOX IS
#  WHAT CAUSED THE REFUSAL.  N absorbs the arbitrary scale of the Galerkin eigenvector, which grows with the box,
#  so a guard comparing N against an absolute number fires more readily the larger the box -- a second defect,
#  independent of the window, and it affects only the GUARD, since N cancels that scale in the orbital itself.
#
#  THE WINDOW DEFECT, MEASURED, AND IT SWITCHES ON EXACTLY WHERE PREDICTED.  300 points at hp = 1e-2 is 2.98 a.u.,
#  which is one half-wavelength at 14.9 eV.  The factor by which an orbital is scaled too high, against a
#  half-wavelength window that must contain a crest:
#
#      E [eV]     1.449   2.500   4.000   6.000   9.000  13.000 | 18.000  25.000  35.000  50.000  80.000 120.000
#      300pts/hw   0.31    0.41    0.51    0.63    0.77    0.93 |   1.09    1.29    1.52    1.82    2.30    2.82
#      box 24     1.2346  1.0000  1.0144  1.0000  1.0000  1.0000|  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
#      box 40     1.0638  1.1809  1.0011  1.0000  1.0109  1.0000|  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
#
#  Below the line the error reaches 23.5 % and differs between boxes with no pattern, which is the signature of a
#  window landing at an arbitrary phase.  ABOVE IT THE FACTOR IS EXACTLY 1.0000 AT EVERY ENERGY AND BOTH BOXES.
#
#  THAT LAST ROW IS WHAT DE-RISKS THE REPAIR.  With `nWin = min(max(nHalf,300), div(length(P),3))` the window can
#  only ever GROW, and it grows only where nHalf > 300, i.e. below 14.9 eV.  Every faster continuum electron in
#  JAC keeps a byte-identical window and therefore a byte-identical orbital, so the test suite can only move where
#  a SLOW continuum electron appears.  That is a measured statement about the blast radius, not a hope.
#  ------------------------------------------------------------------------------------------------------------------
#
#  THE DESIGN.  Three sharings, not nine, because the question is about the INNER integral:  the sharing mesh is
#  the outer one and a 3-point total compared against itself across settings is a valid relative measure, while
#  the middle sharing is EQUAL sharing exactly and is the cleanest single number to read.  One photon energy,
#  200 eV, where Kornberg & Miraglia, PRA 48, 3714 (1993), Table II tabulates sigma+ = 5.44e-02 Mb and
#  sigma2+ (C2) = 1.69e-03 Mb rather than something to interpolate.
#
#    ladder 1   box 24 a.u.,  factor 3,  NoIntermediateEnergies = 9, 17, 25, 33, 49
#    ladder 2   box 24 a.u.,  NoIntermediateEnergies = 33,  factor = 3, 5, 8
#    check      box 40 a.u.,  the best setting of ladder 1, against the same setting at box 24
#
#  Run as:  julia --project=. tools/scanPdiHeliumIntermediate.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")

const OMEGA  = 200.0
const NSHARE = 3
const MAXKAP = 2
const nm     = Nuclear.Model(2.0)
const ionConfigs = [Configuration("1s"), Configuration("2s"), Configuration("2p")]

##  (box, NoIntermediateEnergies, intermediateEnergyFactor).  25 is the largest mesh that survives at factor 3:
##  33 places a point at 0.954 eV, where the pure-sine normalization refuses.  A HIGHER factor raises the whole
##  range and so lifts the lowest point away from the refusal, which is why 33 is kept for factors 5 and 8.
const RUNS = [ (24.0, 25, 5.0), (24.0, 25, 8.0), (24.0, 33, 8.0),
               (40.0,  9, 3.0), (40.0, 25, 3.0) ]

println("="^118)
println("PhotoDoubleIonization, He at $(OMEGA) eV:  is the box dependence an unconverged INTERMEDIATE integral?")
println("$(NSHARE) sharings, maxKappa = $(MAXKAP), ION Green space;  K&M C2 target sigma2+ = 1.69e+03 barn")
println("="^118)

grids = Dict{Float64,Radial.Grid}()
greens = Dict{Float64,Multiplet}()

for (rbox, nInt, fac) in RUNS
    if !haskey(grids, rbox)
        grids[rbox] = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = rbox)
        greens[rbox] = redirect_stdout(devnull) do
            SelfConsistent.performSCF(ionConfigs, nm, grids[rbox], AsfSettings(); printout=false)
        end
        o = greens[rbox].levels[1].basis.orbitals
        @printf("\n  Green space at box %.1f a.u.:  2s %9.5f  2p %9.5f Ha,  %d levels\n",
                rbox, o[Subshell("2s_1/2")].energy, o[Subshell("2p_1/2")].energy, length(greens[rbox].levels))
    end
    grid = grids[rbox];   gMp = greens[rbox]
    t = time()
    dSet = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
               multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
               NoEnergySharings=NSHARE, NoIntermediateEnergies=nInt, intermediateEnergyFactor=fac,
               maxKappa=MAXKAP, calcDifferentialCs=false, printBefore=false,
               lineSelection=LineSelection(true, indexPairs=[(1,0)]),
               eeInteraction=CoulombInteraction(), gMultiplet=gMp )
    println("  >>>> MARK  box=$(rbox)  nInt=$(nInt)  factor=$(fac)")
    ##  A REFUSAL MUST NOT KILL THE REST OF THE LADDER.  The run of 03-Nov lost four rungs and the box-40
    ##  cross-check to one "enlarge box-size" raise on the fourth, which is a harness fault and not a result.
    try
        perform(Atomic.Computation(Atomic.Computation(), name="s2+ n$(nInt) f$(fac)", grid=grid, nuclearModel=nm,
                    initialConfigs=[Configuration("1s^2")], finalConfigs=[Configuration("1s^0")],
                    orbitalSource=Basics.SeparateOrbitals(), processSettings=dSet))
        @printf("  box %.1f  nInt %3d  factor %.1f  computed in %7.1f s\n", rbox, nInt, fac, time()-t)
    catch err
        @printf("  box %.1f  nInt %3d  factor %.1f  REFUSED after %6.1f s:  %s\n",
                rbox, nInt, fac, time()-t, sprint(showerror, err))
    end
    flush(stdout)
end
println("\n" * "="^118)
println("DONE.  Read the per-run tables above.")
println("="^118)
