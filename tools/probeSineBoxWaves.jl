#
# tools/probeSineBoxWaves.jl
#
# WHAT THIS IS FOR.  It is the calibration behind `CONT_BOX_WAVES_FATAL` and `CONT_BOX_WAVES_WARN` in
# `src/module-Continuum.jl`:  how wrong a SINE normalization gets when too few oscillations of the outgoing
# wave fit inside the radial box.  Kept so the two numbers can be re-derived rather than trusted.
#
# THE DESIGN, and both halves were learned by getting them wrong first:
#
#   * A BARE POINT PROTON, so the exact solution IS the Coulomb function.  `CoulombFunctionNorm` is then correct
#     at any radius and serves as the CONTROL -- if its answer also moved with the box, the box, the mesh or the
#     solve would be at fault and the experiment would say nothing about the sine fit.  In series B that control
#     reproduces itself to 7e-8 across a factor seven in box size.
#   * A SCALE-FREE STATISTIC.  The first version compared the normalization CONSTANT N between two boxes, which
#     is meaningless:  the Galerkin eigenvector carries an arbitrary scale set by the basis and N is precisely
#     the number that absorbs it -- its "control" drifted by factors of two.  What is compared instead is the
#     norm of the NORMALIZED large component over a radial window both boxes contain.  After normalization the
#     orbital is on an absolute scale, so the same physical orbital must come out of any adequate box.
#   * ONE `hp` IN BOTH BOXES OF A PAIR.  With the step equal the grids NEST exactly, so the window carries
#     identical mesh points AND identical weights and nothing is interpolated;  the run checks this bit for bit.
#     A first attempt used a different step in the two boxes and the control showed that alone moved the answer
#     by 27 %.
#
# THIS FILE IS SERIES B -- the energy FIXED and only the box varied, so the Sommerfeld parameter eta is constant
# and the oscillations-in-box is the only variable.  Series A varied the energy in a fixed 12 a.u. box, which
# moves eta with it; its result is recorded beside the constants in `src/module-Continuum.jl`.  Measured
# 06-Nov-2026, eta = 1.818:
#
#     waves  0.54   0.71   0.97   1.23   1.76   2.46   3.51   5.27
#     error  74.6%  57.7%  26.2%  10.2%  19.3%   2.6%  13.2%   6.0%
#
# BELOW ONE WAVE THE ERROR IS 26-75 %, which squares into a cross section wrong by a factor of 1.3 to 3, and
# that is where the fatal level sits.  Above one wave the error is an erratic 2.6-19.3 % and does not settle by
# 5.3 waves, and the 240 a.u. reference itself carries a fit residual of 0.337 -- so the warn level of 4 marks
# where the box can no longer be shown to be the cause, not where the answer is known good.
#
# AND THE EXISTING FIT RESIDUAL CANNOT REPLACE THIS, which is why a box tier exists at all:  it is ANTI-correlated
# at the dangerous end -- 0.038 at the 74.6 % error, below its own 0.20 warning, while 2.6 % reads 0.387 and
# warns.  With less than one wave in the box the fit has almost nothing to satisfy, so it matches the stub nearly
# perfectly;  a small residual there means unconstrained, not right.
#
# Usage:   julia --project=. tools/probeSineBoxWaves.jl
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf
using JenaAtomicCalculator: Basics, Continuum, Defaults, Nuclear, Radial, Subshell, Orbital

Defaults.setDefaults("method: continuum, Galerkin")

# SERIES B -- the confound removed.  Series A varied the ENERGY to vary boxWaves, so the Sommerfeld parameter
# eta = Zbar/q varied with it (5.0 down to 0.33) and a reader could fairly say the trend was an eta trend.
# Here the ENERGY IS FIXED and only the box is varied, so eta is constant at 1.82 and boxWaves is the only
# variable.  Same scale-free statistic, same hp throughout, same bare point proton, same control.

const Q      = 0.55                   ## eta = 1/q = 1.818, fixed
const HP     = 0.03
const RBOXES = [6.0, 8.0, 11.0, 14.0, 20.0, 28.0, 40.0, 60.0]
const RREF   = 240.0                  ## 21 wavelengths
const W1     = 3.0
const W2     = 6.0                    ## the outer half of the SMALLEST box, so every grid contains it

function solveOn(grid, energy, kappa)
    nm  = Nuclear.Model(1.0, Nuclear.PointNucleus())
    pot = Nuclear.nuclearPotential(nm, grid)
    sh  = Subshell(101, kappa);    st = Continuum.Settings(false, grid.NoPoints - 10)
    orb = redirect_stdout(devnull) do;  Continuum.generateOrbitalGalerkin(energy, sh, pot, st)  end
    return( orb, pot, st )
end

function windowNorm(P, grid, r1, r2)
    s = 0.0
    for k = 1:min(length(P), length(grid.r))
        if  r1 <= grid.r[k] <= r2    s += P[k]^2 * grid.wr[k]    end
    end
    return( sqrt(s) )
end

function oneConvention(conv, grid, energy, kappa, r1, r2)
    orb, pot, st = solveOn(grid, energy, kappa)
    Defaults.setDefaults("method: normalization, " * conv)
    if  conv == "pure sine"
        Zbar = -Radial.determineZbar(pot)
        o, phi, nf = redirect_stdout(devnull) do
            Continuum.normalizeOrbitalPureSine(orb, grid, st; Zbar = Zbar)   end
        wa = Continuum.asymptoticAmplitude(orb.P, grid, size(orb.P,1), energy, Continuum.asymptoticCharge(Zbar))
        return( (S = windowNorm(o.P, grid, r1, r2), resid = wa.residual) )
    else
        o, ph, N = redirect_stdout(devnull) do
            Continuum.normalizeOrbitalCoulombFunctions(orb, pot, st)   end
        return( (S = windowNorm(o.P, grid, r1, r2), resid = NaN) )
    end
end

energy = Q*Q/2;    lam = 2pi/Q
@printf("#\n# SERIES B:  energy FIXED at %.5f Ha (lambda = %.3f, eta = %.3f);  only the box varies.\n", energy, lam, 1/Q)
@printf("# hp = %.3f throughout;  reference box %.0f a.u. (%.1f waves);  window %.1f-%.1f a.u.\n#\n",
        HP, RREF, RREF/lam, W1, W2)
gR = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = HP, rbox = RREF)
bS = oneConvention("pure sine",         gR, energy, -1, W1, W2)
bC = oneConvention("Coulomb functions", gR, energy, -1, W1, W2)
@printf("reference:  S_sine = %.5e   S_CF = %.5e   (resid %.3f, %d pts)\n\n", bS.S, bC.S, bS.resid, gR.NoPoints)
@printf("%8s %7s %12s %10s %8s %10s %6s\n", "rbox", "boxWav", "S_sine", "err_sine", "resid", "err_CF", "nest")
println(repeat("-", 70))
for rb in RBOXES
    g = Radial.Grid(Radial.Grid(false); rnt = 2.0e-6, h = 5.0e-2, hp = HP, rbox = rb)
    nest = all(g.r[k] == gR.r[k] for k = 1:g.NoPoints) && all(g.wr[k] == gR.wr[k] for k = 1:g.NoPoints)
    try
        aS = oneConvention("pure sine",         g, energy, -1, W1, W2)
        aC = oneConvention("Coulomb functions", g, energy, -1, W1, W2)
        @printf("%8.1f %7.3f %12.5e %10.3e %8.3f %10.3e %6s\n", g.r[g.NoPoints], g.r[g.NoPoints]/lam, aS.S,
                abs(aS.S-bS.S)/abs(bS.S), aS.resid, abs(aC.S-bC.S)/abs(bC.S), nest ? "yes" : "NO")
    catch e
        msg = e isa ErrorException ? first(split(e.msg, "\n")) : string(typeof(e))
        @printf("%8.1f %7.3f   RAISED: %s %6s\n", g.r[g.NoPoints], g.r[g.NoPoints]/lam, first(msg, 44),
                nest ? "yes" : "NO")
    end
end
