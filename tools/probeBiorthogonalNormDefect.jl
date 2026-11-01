#
#  [BIO-NORM]  What the biorthogonal transformation does to the ORBITAL NORMS, and what that cost in an
#              N-vs-(N-1) Auger rate.  The demonstration that found the defect repaired on 01-Nov-2026.
#
#  WHY THIS EXISTS.  `BiOrthogonal.computeTransformation` returns correctly biorthogonal orbitals --
#  `<newLeft_i|newRight_j> = delta_ij` -- and that condition does NOT fix the pair of transformation matrices.
#  The unpivoted LU used until 01-Nov-2026 chose the asymmetric pair:  `Cleft = inv(L)'` is unit triangular and
#  leaves every norm alone, while `Cright = inv(U)` carries the WHOLE scaling.  A rotation that mixes two
#  subshells which are BOTH closed and BOTH occupation 2 in EVERY CSF of BOTH multiplets is a relabelling with
#  provably zero physical consequence -- and it came back with right-hand norms of 1/cos(theta).
#
#  THAT WAS NOT COSMETIC, because the N-vs-(N-1) consumers REBUILD A MEAN FIELD from the transformed ionic core
#  in order to generate their continuum orbital.  The same-N consumers (PhotoEmission, PhotoExcitation,
#  GeneralizedOscillatorStrength) never do, which is exactly why they validated and this did not.
#
#  MEASURED on Mg, Z = 12, initial [Ne] 3p^2 + [Ne] 3s^2 (N = 12) against final [Ne] 3s (N = 11), 1s/2s mixed by
#  theta = 0.37 on the final side, every multiplet on ONE frozen orbital set so that the overlap is EXACTLY the
#  identity without the rotation and EXACTLY a rotation with it:
#
#      control line 2+ -> 1/2+        reference    rot, no bio   no rot, bio    rot + bio
#      unpivoted LU (before)        7.6421837e-2  7.6421837e-2  7.6421837e-2  1.4118321e-1   (+84.7 %)
#      symmetric split (after)      7.6421837e-2  7.6421837e-2  7.6421837e-2  7.6421837e-2   (+5.6e-10)
#
#  THE CONTROL IS CHOSEN, NOT CONVENIENT: that line's initial CSF is pure 3p^2 (J = 2, so the J = 0 padding CSF
#  cannot mix into it), so by the ordinary two-orbital Slater-Condon rule it CANNOT depend on the 1s/2s
#  spectators -- and it is bit-identical under the rotation whenever the transformation is off.
#
#  AND THE ATTRIBUTION WAS FORCED RATHER THAN ARGUED, by three facts printed below.  (1) `Mleft` and `Mright`,
#  the counter-rotating CI matrices, are BIT-IDENTICAL with and without the rotation -- the +-log(cos theta) that
#  the LU leaves on the diagonal of `log(Cright)` cancels exactly, because the two closed subshells carry equal
#  one-particle coefficients.  So `exp(-Ghat)` was never the problem.  (2) With a frozen set and the rotation
#  confined to (1s,2s), `Cright[1,3] = Cright[2,3] = 0`, so `newRight_3s` is EXACTLY unchanged, and the 3p
#  orbitals sit in other kappa blocks -- none of the control line's ACTIVE orbitals move.  (3) The transformed
#  core's density carried 11.300876 electrons instead of 11 and its DFS screening shifted by ~0.3 charge units at
#  every radius.  The continuum orbital was therefore the only channel left.
#
#  THE DECIDING MEASUREMENT WAS THE SCALING.  The orbital MIXING is O(theta);  the norm defect 1/cos(theta) - 1 is
#  O(theta^2).  The rate error followed theta^2 -- err/theta^2 = 6.19, 5.63, 5.42, 5.37 at theta = 0.37, 0.185,
#  0.0925, 0.04625, converging rather than doubling -- and err/(1/cos-1) converged to 10.7.  The error is linear
#  in the NORM defect, not in the mixing.
#
#  TWO DEAD ENDS, recorded so they are not retried.  The hypothesis of record (27-Oct-2026) was that a
#  single-CSF final multiplet gives `generateCounterRotatingCiMatrices` no CSF to absorb one-particle mixing into,
#  so `exp(-Ghat)` produces a spurious scalar on a 1x1 block:  REFUTED by (1) above.  And the balanced symmetric
#  split `Cleft = U*Sigma^-1/2`, `Cright = V*Sigma^-1/2` is NOT usable here:  the singular vectors are fixed only
#  up to `(U, V) -> (U W, V W)` for orthogonal W commuting with Sigma, which for a rotation is complete freedom,
#  so the factors come back as an arbitrary orthogonal matrix and RELABEL the orbitals -- measured, that moved the
#  same rate by a factor 31000 and produced a reflection whose logarithm is complex.  Appending `U'` to both makes
#  each factor invariant under W;  that is what shipped.
#
#  THE SUITE GUARDS THE INVARIANT, NOT THIS RATE:  `TestFrames.testModule_BiOrthogonalRotation` asserts the norm
#  preservation in milliseconds, where this file costs two continuum computations.  Keep both -- the cheap test
#  cannot drift away from the cause, and this file is what makes the cause legible.
#
#  RUN:  julia --project=. tools/probeBiorthogonalNormDefect.jl
#
using JenaAtomicCalculator, Printf, LinearAlgebra
const JAC = JenaAtomicCalculator

grid  = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
nm    = Nuclear.Model(12.);    asf = AsfSettings();    theta = 0.37

##  ONE frozen orbital set for every multiplet, so the overlap is exactly I without the rotation and exactly a
##  rotation with it;  the SCF must run on one electron number, so the 12-electron set supplies the orbitals and
##  the 11-electron CI reuses them
orbitals = SelfConsistent.performSCF([Configuration("[Ne] 3p^2"), Configuration("[Ne] 3s^2")], nm, grid, asf;
                                     printout=false).levels[1].basis.orbitals
iMp = Hamiltonian.performCIwithFrozenOrbitals([Configuration("[Ne] 3p^2"), Configuration("[Ne] 3s^2")],
                                              orbitals, nm, grid, asf; printout=false)
fMp = Hamiltonian.performCIwithFrozenOrbitals([Configuration("[Ne] 3s")], orbitals, nm, grid, asf; printout=false)

"""
`rotateClosedPair(mp::Multiplet, th::Float64, grid::Radial.Grid)`
    ... mixes the 1s and 2s orbitals of the given multiplet's basis by the angle `th`.  Both are kappa = -1 and
        both carry occupation 2 in every CSF, so this is a relabelling with no physical consequence whatever;
        a multiplet::Multiplet on the rotated basis is returned.
"""
function rotateClosedPair(mp::Multiplet, th::Float64, grid::Radial.Grid)
    b  = mp.levels[1].basis;    o = copy(b.orbitals)
    s1 = Subshell(1,-1);        s2 = Subshell(2,-1);    a = o[s1];    c = o[s2]
    n  = min(length(a.P), length(c.P))
    mix(x, y, ca, cb) = ca .* x[1:n] .+ cb .* y[1:n]
    o[s1] = Radial.Orbital(s1, a.isBound, a.useStandardGrid, a.energy, mix(a.P,c.P,cos(th),sin(th)),
                mix(a.Q,c.Q,cos(th),sin(th)), mix(a.Pprime,c.Pprime,cos(th),sin(th)), mix(a.Qprime,c.Qprime,cos(th),sin(th)), grid)
    o[s2] = Radial.Orbital(s2, c.isBound, c.useStandardGrid, c.energy, mix(a.P,c.P,-sin(th),cos(th)),
                mix(a.Q,c.Q,-sin(th),cos(th)), mix(a.Pprime,c.Pprime,-sin(th),cos(th)), mix(a.Qprime,c.Qprime,-sin(th),cos(th)), grid)
    nb = Basis(b.isDefined, b.NoElectrons, b.subshells, b.csfs, b.coreSubshells, o)

    return( Multiplet(mp.name, [Level(l.J, l.M, l.parity, l.index, l.energy, l.relativeOcc, l.hasStateRep, nb, l.mc)
                                for l in mp.levels]) )
end

##  (1) the counter-rotating CI matrices, and the orbital norms, with and without the rotation
println("\n######## the transformation itself ########")
for  (tag, fM)  in  ( ("no rotation", fMp), ("rotated 1s/2s by theta = $theta", rotateClosedPair(fMp, theta, grid)) )
    lB = iMp.levels[1].basis
    (rB, _, _) = BiOrthogonal.alignOrbitalSigns(lB, fM.levels[1].basis, grid)
    tr = BiOrthogonal.computeTransformationMatrices(lB, rB, grid)
    Ml = BiOrthogonal.generateCounterRotatingCiMatrices(lB, tr, :left)
    Mr = BiOrthogonal.generateCounterRotatingCiMatrices(rB, tr, :right)
    (_, nR)  = BiOrthogonal.generateBiorthogonalShellMatrices(lB, rB, grid)
    worst    = 0.;    label = ""
    for  (sh, orb)  in  nR
        dev = abs( sqrt(JAC.RadialIntegrals.overlap(orb, orb, grid)) - 1. )
        dev > worst  &&  (worst = dev;   label = string(sh))
    end
    nEff = 0.
    for  (sh, occ)  in  Basics.extractMeanOccupation(fM.levels[1].basis)
        nEff = nEff + occ * JAC.RadialIntegrals.overlap(nR[sh], nR[sh], grid)
    end
    @printf("   %-34s  Mright = %+.9f   worst right-hand norm departure %.3e (%s)   N_eff = %.8f\n",
            tag, Mr[1,1], worst, label, nEff)
end

##  (2) the four cells of the Auger comparison
println("\n######## the control line, four cells ########")
@printf("   %-7s %-7s %-12s %-24s\n", "rot", "bio", "line", "total rate [a.u.]")
for  (rot, bio)  in  ( (false,false), (true,false), (false,true), (true,true) )
    fM  = rot ? rotateClosedPair(fMp, theta, grid) : fMp
    set = AutoIonization.Settings(AutoIonization.Settings(); calcBiorthogonal=bio, printBefore=false)
    lines = redirect_stdout(devnull) do
        AutoIonization.computeLines(fM, iMp, nm, grid, set; output=true, printout=false)
    end
    for  l  in  lines
        symi = string(LevelSymmetry(l.initialLevel.J, l.initialLevel.parity))
        symf = string(LevelSymmetry(l.finalLevel.J,   l.finalLevel.parity))
        @printf("   %-7s %-7s %-12s %.12e\n", rot, bio, "$symi->$symf", l.totalRate)
    end
end
println("\n   The 2+ -> 1/2+ line must read the same in all four cells.  A rot/bio cell that departs from the")
println("   reference says the transformation is no longer norm-preserving, i.e. the repair of 01-Nov-2026 is gone.")
