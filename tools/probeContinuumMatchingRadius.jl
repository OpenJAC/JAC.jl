#
# tools/probeContinuumMatchingRadius.jl
#
# WHAT THIS IS FOR.  `Continuum.normalizeOrbitalCoulombFunctions` -- the SHIPPED DEFAULT normalization -- pins the
# numerical orbital to the exact Coulomb (or Bessel) functions at ONE radius, `pot.grid.r[mtp]`, by matching the
# value and the derivative there.  Nothing checks whether that one radius was a FAIR place to do it:  if the
# potential has not yet settled onto its long-range -Zbar/r form at that point, or if the orbital is corrupted
# there, the match still succeeds and returns a confident number.  The register records this as the one
# continuum convention with no form check, guarded only by the interior/asymptotic ratio tier -- which is
# measured BLIND to this class of failure.
#
# This harness serves two purposes, in this order:
#
#   1. BASELINE.  Printed first, as a table of the three quantities the routine itself reports -- the
#      normalization constant N, the inner phase iPhase and the Coulomb phase cPhase -- over a spread of cases.
#      A restructuring of that routine must reproduce this table CHARACTER FOR CHARACTER.  That is the control:
#      without it, "the refactor changed nothing" is an opinion rather than a measurement.
#
#   2. CALIBRATION.  The same cases, with the radius-stability measure reported per case, so that a threshold
#      can be set from what SOUND orbitals actually do rather than chosen.  The companion mutation run
#      (--mutate) ruins orbitals deliberately and reports the same measure, which is how the existing continuum
#      guards were calibrated (interior/asymptotic over 264 cases; the asymptotic-fit residual over 16 sound
#      and 7 ruined).
#
# Usage:   julia --project=. tools/probeContinuumMatchingRadius.jl [--mutate]
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf
using JenaAtomicCalculator: Basics, Continuum, Defaults, Nuclear, Radial, Subshell, Orbital

MUTATE = "--mutate" in ARGS

Defaults.setDefaults("method: continuum, Galerkin")
## NOTE THE KEY.  "pure Coulomb" selects CoulombSineNorm, a DIFFERENT convention; the one wanted here
## is "Coulomb functions".  The two names are one word apart and the wrong one runs clean.
Defaults.setDefaults("method: normalization, Coulomb functions")

# THE CASE LIST spans what the convention is actually used on:  a near-neutral system where Zbar ~ 1 and the
# Bessel branch is taken, up to a bare heavy nucleus where eta is large and the asymptotic sine form is
# meaningless -- which is precisely the regime this convention exists to cover and the sine conventions do not.
# Energies run from just above threshold to fast, since the number of mesh points per oscillation falls with
# energy and that is the dominant failure mode of the whole continuum machinery.
cases = NamedTuple[]
for  Z        in  [2.0, 18.0, 36.0, 92.0]
    for  en   in  [0.05, 0.5, 5.0, 50.0]            ## Hartree
        for  kappa  in  [-1, 1, -2]
            push!( cases, (Z=Z, energy=en, kappa=kappa) )
        end
    end
end

println("#")
println("# probeContinuumMatchingRadius  --  ", MUTATE ? "MUTATION run (orbitals deliberately ruined)" : "BASELINE run (sound orbitals)")
println("#")
@printf("%6s %9s %6s %18s %15s %14s %10s %9s %9s %6s\n", "Z", "E [Ha]", "kap", "N", "phase",
        "verdict", "spread", "jitter", "sDrift", "ppo")
println(repeat("-", 125))


# ---------------------------------------------------------------------------------------------------------------
# THE MUTATION PASS.  A threshold is only defensible if it SEPARATES, so sound orbitals alone are not enough:  the
# same measure is applied to orbitals ruined on purpose, in the four ways a real solve goes wrong.  This is the
# discipline the sibling continuum guards were calibrated with (interior/asymptotic over 264 cases; the
# asymptotic-fit residual over 16 sound and 7 ruined).
#
# `spike` is included KNOWING this guard may miss it, and that is the point:  an interior spike is what
# `Continuum.checkContinuumAmplitude` exists to catch, and the register records the two tiers as complementary
# rather than alternatives.  A guard that claimed to catch everything would be the suspicious one.
function mutateOrbital(o::Orbital, how::Symbol, grid::Radial.Grid)
    P = copy(o.P)  ;  Pp = copy(o.Pprime)  ;  n = length(P)
    if      how == :ramp                     ## a slow multiplicative distortion of the outer half
        for i = div(n,2):n
            f = 1.0 + 0.30 * (i - div(n,2)) / (n - div(n,2))
            P[i] *= f  ;  Pp[i] *= f
        end
    elseif  how == :noise                    ## point-to-point noise, as an under-resolved mesh gives
        for i = 1:n
            f = 1.0 + 0.02 * sin(37.0 * i)   ## deterministic, so the calibration is reproducible
            P[i] *= f  ;  Pp[i] *= f
        end
    elseif  how == :spike                    ## one enormous interior value
        P[div(n,20)] = 1.0e3 * maximum(abs.(P))
    elseif  how == :derivative               ## P and Pprime made mutually inconsistent, growing outward
        for i = 1:n
            Pp[i] *= 1.0 + 0.20 * i / n
        end
    end
    return( Orbital(o.subshell, o.isBound, o.useStandardGrid, o.energy, P, o.Q, Pp, o.Qprime, o.grid) )
end

if  MUTATE
    @printf("%6s %9s %6s %12s %14s %10s %9s %9s\n", "Z", "E [Ha]", "kap", "mutation", "verdict", "spread",
            "jitter", "sDrift")
    println(repeat("-", 92))
    for  c  in  [ (Z=2.0, energy=5.0, kappa=-1), (Z=18.0, energy=5.0, kappa=-1), (Z=18.0, energy=50.0, kappa=1),
                  (Z=36.0, energy=5.0, kappa=-2), (Z=92.0, energy=50.0, kappa=-1) ]
        q      = sqrt(2c.energy)  ;  lambda = 2pi / q
        rbox   = max(25.0, 10.0 * lambda)
        grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = lambda/40, rbox = rbox)
        nm     = Nuclear.Model(c.Z)  ;  pot = Nuclear.nuclearPotential(nm, grid)
        sh     = Subshell(101, c.kappa)  ;  st = Continuum.Settings(false, grid.NoPoints - 10)
        l      = Basics.subshell_l(sh)   ;  Zbar = -Radial.determineZbar(pot)
        raw    = redirect_stdout(devnull) do;  Continuum.generateOrbitalGalerkin(c.energy, sh, pot, st)  end
        for  how  in  [:none, :ramp, :noise, :spike, :derivative]
            o  = how == :none ? raw : mutateOrbital(raw, how, grid)
            wa = Continuum.checkMatchingRadius(o, pot, l, Zbar, size(o.P,1))
            @printf("%6.1f %9.3f %6d %12s %14s %10.2e %9.1f %9.2e\n", c.Z, c.energy, c.kappa, string(how),
                    string(wa.verdict), wa.spread, wa.jitter, wa.scaleDrift)
        end
    end
    exit(0)
end

for  c  in  cases
    nm   = Nuclear.Model(c.Z)
    # BOTH THE BOX AND THE MESH ARE SIZED TO THE WAVELENGTH, which is the only sizing that means anything for a
    # continuum orbital:  the box must hold enough oscillations for an asymptotic statement to be possible, and
    # the mesh must resolve each one.  40 points per oscillation is the measured sound value (22-Oct-2026:
    # 15.65 fails catastrophically, ~40 is sound), so hp follows from the wavelength rather than being fixed.
    # A FIXED hp IS WHAT MADE THE FIRST VERSION OF THIS HARNESS DIE:  hp = 5e-3 with a 497 a.u. box at 0.05 Ha
    # is ~100 000 mesh points, and the B-spline solve ran the machine out of memory.
    q      = sqrt(2c.energy)
    lambda = 2pi / q
    rbox   = max(25.0, 10.0 * lambda)
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = lambda/40, rbox = rbox)
    pot  = Nuclear.nuclearPotential(nm, grid)
    sh   = Subshell(101, c.kappa)
    st   = Continuum.Settings(false, grid.NoPoints - 10)
    try
        ## The normalization is called DIRECTLY rather than through generateOrbitalLocalPotential, for two
        ## reasons:  that wrapper hardcodes normFactor = 1.0 for this convention, so its returned factor
        ## carries no information;  and the direct call hands back N and the phase as Float64 at full
        ## precision, which is a stricter baseline than re-parsing the %.4e line the routine prints.
        raw  = redirect_stdout(devnull) do
            Continuum.generateOrbitalGalerkin(c.energy, sh, pot, st)
        end
        ## the guard's own verdict is reported beside the numbers, since THAT is the calibration datum:  a
        ## threshold must be read off what SOUND orbitals do, never chosen in advance.
        l    = Basics.subshell_l(sh)
        Zbar = -Radial.determineZbar(pot)
        mtp  = size(raw.P, 1)
        wa   = Continuum.checkMatchingRadius(raw, pot, l, Zbar, mtp)
        cOrb, phase, N = redirect_stdout(devnull) do
            Continuum.normalizeOrbitalCoulombFunctions(raw, pot, st)
        end
        @printf("%6.1f %9.3f %6d %18.10e %15.7e %14s %10.2e %9.2f %9.2e %6.1f\n",
                c.Z, c.energy, c.kappa, N, phase, string(wa.verdict), wa.spread, wa.jitter, wa.scaleDrift,
                wa.ppo)
    catch e
        msg = e isa ErrorException ? first(split(e.msg, "\n")) : string(typeof(e))
        @printf("%6.1f %10.3f %7d   RAISED: %s\n", c.Z, c.energy, c.kappa, first(msg, 48))
    end
end
