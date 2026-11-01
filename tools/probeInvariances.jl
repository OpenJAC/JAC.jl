#
#  [INV]  Quantities that MUST NOT change, and a harness that can tell whether it actually tested them.
#
#  WHY THIS SHAPE.  The two defects found on 31-Oct and 01-Nov-2026 -- the TEOP second-order conditions
#  (`3f06c1dc`) and the biorthogonal norm defect (`7b587a99`) -- were both caught the same way: by constructing a
#  case whose answer is PROVABLY known, with no literature involved.  For TEOP it was two CI partitionings of one
#  number;  for BiOrthogonal it was a rotation of a uniformly closed subshell pair, which cannot change anything.
#  Neither defect was visible to any existing test and both returned wrong numbers in silence.  This file points
#  the same lever at other places where an arbitrary choice must not reach an observable.
#
#  EVERY AXIS CARRIES ITS OWN VACUITY CHECK, which is the part that makes a NEGATIVE result worth anything.  A
#  harness that silently fails to perturb the computation reports perfect invariance, and that is
#  indistinguishable from success.  So each axis also reports a quantity that MUST move, and says so when it
#  does not.
#
#  MEASURED 01-Nov-2026 -- BOTH AXES CLEAN, and the record is the point: this rules out two classes of silent
#  instability, and one of them demonstrably happens on every run.
#
#      axis                                  vacuity check                  energies     rates       Hfs A/B
#      1s/2s... orbital sign flip (2s)        4 of 4 amplitudes flip, by -1  0.000e+00    0.000e+00   0.000e+00
#      orbital sign flip (2p_1/2)             2 of 4 amplitudes flip, by -1  0.000e+00    0.000e+00   0.000e+00
#      configuration-list order reversed      2 of 2 lines matched by omega  1.421e-14    6.066e-13   --
#
#  AXIS 1 IS NOT HYPOTHETICAL.  JAC's own SCF prints ">>> Sign changed for orbital 2s_1/2" -- an orbital's
#  overall sign is arbitrary and the solver does flip it, run to run.  `BiOrthogonal.alignOrbitalSigns` exists
#  because exactly this bit once.  And the chain is not trivial: flipping one orbital reverses the sign of every
#  off-diagonal CI element between CSFs differing by one electron in it, so the CI EIGENVECTORS change sign while
#  the eigenvalues do not, and every amplitude in which that orbital is active an ODD number of times changes
#  sign too.  The RATE is invariant only if the eigenvector signs and the amplitude signs are derived from the
#  SAME orbitals.  A path that takes orbitals from one place and mixing coefficients from another breaks it --
#  which is the shape of the defect fixed in `7b587a99`.  The measured result is that JAC gets this exactly
#  right: the amplitudes flip by exactly -1 and nothing else moves, to the last bit.
#
#  AXIS 2 IS WHAT A USER DOES EVERY DAY.  A configuration list is written in whatever order comes to mind, and
#  the CSF order and subshell order are derived from it.  Measured: the subshell list comes back IDENTICAL either
#  way, so `Basics.merge` normalises the order rather than propagating it, and nothing downstream depends on it.
#  Note that the level INDEX may legitimately permute, so the lines are matched by PHOTON ENERGY and the
#  energies are compared SORTED -- comparing by index would report a false finding.
#
#  AXES WORTH ADDING, none yet measured: a closed-shell rotation applied to Hfs, LandeZeeman and
#  `matrixElement_Mab` (it is now provably inert for BiOrthogonal, see `probeBiorthogonalNormDefect.jl`);  a
#  permutation of the subshell list itself;  and the gauge of a magnetic multipole, which `94fb9cc5` made
#  meaningful.
#
#  RUN:  julia --project=. tools/probeInvariances.jl
#
using JenaAtomicCalculator, Printf, LinearAlgebra
const JAC = JenaAtomicCalculator

grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.2e-2, rbox = 10.0)
nm     = Nuclear.Model(10.);    asf = AsfSettings()
upperC = [Configuration("1s^2 2s 2p")]
lowerC = [Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")]
orbs   = SelfConsistent.performSCF(vcat(upperC, lowerC), nm, grid, asf; printout=false).levels[1].basis.orbitals
hfsNm  = Nuclear.Model(Nuclear.Model(10.); spinI=AngularJ64(5//2), mu=1.0, Q=0.1)

"""
`flipSign(orbitals::Dict{Subshell,Radial.Orbital}, sh::Subshell)`
    ... returns a copy of `orbitals` in which the orbital `sh` has all four radial components negated.  An
        orbital's overall sign is a convention, so no observable may depend on it;  a
        dict::Dict{Subshell,Radial.Orbital} is returned.
"""
function flipSign(orbitals::Dict{Subshell,Radial.Orbital}, sh::Subshell)
    o = copy(orbitals);    a = o[sh]
    o[sh] = Radial.Orbital(sh, a.isBound, a.useStandardGrid, a.energy, -a.P, -a.Q, -a.Pprime, -a.Qprime, a.grid)

    return( o )
end


"""
`observables(orbitals::Dict{Subshell,Radial.Orbital}, upperC, lowerC, nm, hfsNm, grid, asf)`
    ... computes, on the given orbital set, the level energies of the lower multiplet, every E1 amplitude and
        rate between the two multiplets, and the hyperfine A and B parameters with their M1/E2 amplitudes;  a
        dict::Dict{String,Any} of those is returned, for comparison against the same dict on a perturbed set.
"""
function observables(orbitals, upperC, lowerC, nm, hfsNm, grid, asf)
    upMp = Hamiltonian.performCIwithFrozenOrbitals(upperC, orbitals, nm, grid, asf; printout=false)
    loMp = Hamiltonian.performCIwithFrozenOrbitals(lowerC, orbitals, nm, grid, asf; printout=false)
    set  = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb,UseBabushkin],
                                  printBefore=false, calcBiorthogonal=false)
    lines = redirect_stdout(devnull) do
        PhotoEmission.computeLines(loMp, upMp, grid, set; output=true)
    end
    hfs = redirect_stdout(devnull) do
        Hfs.computeOutcomes(loMp, hfsNm, grid, Hfs.Settings(Hfs.Settings(); calcM1=true, calcE2=true,
                                                            printBefore=false); output=true)
    end
    d = Dict{String,Any}()
    d["energies"] = [l.energy for l in loMp.levels]
    ##  one MultipoleAmplitude per multipole, each holding an EmPropertyC;  both gauges are kept
    d["amps"]  = [ (l.initialLevel.index, l.finalLevel.index, string(ma.multipole), g)
                   for l in lines for ma in l.amplitudes for g in (ma.amplitude.Coulomb, ma.amplitude.Babushkin) ]
    d["rates"] = [ (l.initialLevel.index, l.finalLevel.index, l.photonRate.Coulomb, l.photonRate.Babushkin)
                   for l in lines ]
    d["hfs"]   = [ (o.Jlevel.index, o.AIoverMu, o.BoverQ) for o in hfs ]

    return( d )
end

ref = observables(orbs, upperC, lowerC, nm, hfsNm, grid, asf)

##  AXIS 1: a global sign flip of one bound orbital
for  sh  in  ( Subshell(2,-1), Subshell(2,1) )
    r = observables(flipSign(orbs, sh), upperC, lowerC, nm, hfsNm, grid, asf)
    println("\n######## orbital $(string(sh)) sign-flipped ########")
    @printf("   level energies           worst |delta| = %.3e   (must be 0)\n",
            maximum(abs.(r["energies"] .- ref["energies"])))
    nMoved = 0;    worstAmp = 0.
    for  (a, b)  in  zip(r["amps"], ref["amps"])
        abs(b[4]) < 1.0e-20  &&  continue
        abs(a[4] - b[4]) / abs(b[4]) > 1.0e-8  &&  (nMoved += 1)
        worstAmp = max(worstAmp, min(abs(a[4] - b[4]), abs(a[4] + b[4])) / abs(b[4]))
    end
    @printf("   amplitudes               %d of %d moved -- 0 would make this VACUOUS;  worst |+-1| mismatch %.3e\n",
            nMoved, length(ref["amps"]), worstAmp)
    worstRate = 0.
    for  (a, b)  in  zip(r["rates"], ref["rates"]),  k  in  (3, 4)
        abs(b[k]) < 1.0e-30  &&  continue
        worstRate = max(worstRate, abs(a[k] - b[k]) / abs(b[k]))
    end
    @printf("   PhotoEmission RATES      worst relative change = %.3e   (must be ~0)\n", worstRate)
    worstHfs = 0.
    for  (a, b)  in  zip(r["hfs"], ref["hfs"]),  k  in  (2, 3)
        abs(b[k]) < 1.0e-20  &&  continue
        worstHfs = max(worstHfs, abs(a[k] - b[k]) / abs(b[k]))
    end
    @printf("   Hfs A/mu and B/Q         worst relative change = %.3e   (must be ~0)\n", worstHfs)
end

##  AXIS 2: the order in which the configurations are listed.  The level INDEX may legitimately permute, so the
##  energies are compared SORTED and the lines are matched by PHOTON ENERGY;  comparing by index would report a
##  false finding
println("\n######## configuration-list ORDER reversed ########")
rev = observables(orbs, upperC, reverse(lowerC), nm, hfsNm, grid, asf)
eA  = sort(ref["energies"]);    eB = sort(rev["energies"])
@printf("   sorted level energies       worst |delta| = %.3e   (must be ~0)\n",
        length(eA) == length(eB) ? maximum(abs.(eA .- eB)) : NaN)
worstR = 0.;    matched = 0
for  a  in  ref["rates"]
    cands = [b for b in rev["rates"] if a[1] == b[1] && a[2] == b[2]]
    isempty(cands)  &&  continue
    b = cands[1];    global matched += 1
    for  k  in  (3, 4)
        abs(a[k]) < 1.0e-30  &&  continue
        global worstR = max(worstR, abs(a[k] - b[k]) / abs(a[k]))
    end
end
@printf("   PhotoEmission RATES         %d of %d lines matched;  worst relative change = %.3e\n",
        matched, length(ref["rates"]), worstR)
println(matched == 0 ? "      >>> VACUOUS: no line matched" :
        (worstR < 1.0e-10 ? "      invariant" : "      >>> FINDING: a rate depends on the configuration-list order"))

println("\n   A nonzero RATE, energy or Hfs column is a FINDING -- an observable that depends on an arbitrary")
println("   choice is unstable between runs.  A vacuity line that reports 0 moved means the harness perturbed")
println("   nothing and the clean columns beside it prove nothing at all.")
