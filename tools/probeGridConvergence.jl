#
# tools/probeGridConvergence.jl
#
# WHAT THIS IS FOR.  It is the calibration behind `Basics.checkGridConvergence`'s default tolerance, kept so that the
# number can be re-derived rather than trusted.  That check re-performs a computation with the outer step `hp` refined
# and the box held, and calls the grid converged when the compared quantity moves by less than `tolerance`.
#
# THE THRESHOLD IS READ OFF A TABLE ALREADY IN THE REPOSITORY, not chosen.  The second measured table in the docstring
# of `Basics.recommendedGrid` gives, for RR into 1s 4s of H-like Bi with the box fixed at 142 a.u.,
#
#     hp        0.472        0.200        0.100        0.050        0.025
#     sigma_C   3.8338e-10   5.0133e-09   9.7472e-03   1.1324e-02   1.1324e-02   barn
#
# Read as this check reads it -- each grid against its OWN halving -- those are relative changes of 1.00, 1.00, 0.139
# and 0.000.  So every threshold in (0, 0.139) separates the three unconverged grids from the converged one, and none
# outside it does;  1.0e-2 sits inside that interval.  Note also what the same table rules out:  the two gauges agree
# to better than 2 % at EVERY row, so a gauge comparison is blind here while the answer is wrong by seven orders.
#
# WHAT THIS SCRIPT ADDS is the demonstration that arithmetic on a recorded table cannot give:  that the check actually
# FLAGS a grid which passes every a-priori guard JAC has, and PASSES one that is converged.  Ne 2p photoionization at a
# 10 Hartree electron energy, lambda = 1.4050 a.u., box 30 a.u.:
#
#     hp 0.080 -> 0.040    17.6 points per oscillation    4.27999909e-03 -> 4.21612976e-03    1.49e-02   NOT CONVERGED
#     hp 0.030 -> 0.015    47   points per oscillation    4.21651553e-03 -> 4.21659879e-03    1.98e-05   CONVERGED
#
# The coarse grid carries 17.6 points per oscillation against `Continuum.gridConsistency`'s refusal at 15, and 693 mesh
# points against its refusal at 600, so it passes everything -- and its cross section is 1.5 % high.  That gap, three
# orders wide, is what makes the verdict usable.
#
# Usage:   julia --project=. tools/probeGridConvergence.jl
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf
using JenaAtomicCalculator: Basics, Radial, Nuclear, Atomic, Configuration, PhotoIonization,
                            ResultKeys, UseCoulomb, UseBabushkin
using JenaAtomicCalculator: setDefaults

setDefaults("unit: energy", "Hartree")

# One pair:  the given grid against its own halving, with the box held.  The Babushkin cross sections of all
# lines are compared, so the verdict is the worst of them.
function onePair(hp::Float64, label::String)
    ps   = PhotoIonization.Settings(PhotoIonization.Settings(); gauges=[UseCoulomb, UseBabushkin],
                                    electronEnergies=[10.0], calcAnisotropy=false, calcPartialCs=false,
                                    printBefore=false)
    grid = Radial.Grid(Radial.Grid(false); rnt = 4.0e-6, h = 5.0e-2, hp = hp, rbox = 30.0)
    comp = Atomic.Computation(Atomic.Computation(); name = label, grid = grid,
                              nuclearModel = Nuclear.Model(10.),
                              initialConfigs = [Configuration("1s^2 2s^2 2p^6")],
                              finalConfigs   = [Configuration("1s^2 2s^2 2p^5")],
                              processSettings = ps)
    println("\n################  ", label, "   hp = ", hp, "  ################")
    res = Basics.checkGridConvergence(comp, r -> [ln.crossSection.Babushkin for ln in r[PhotoIonization.Settings]])
    @printf(">>>> %-22s converged = %-5s  relative = %.4e  (%d values)\n",
            label, string(res.converged), res.relative, length(res.coarse))
    return( res )
end

wa = onePair(0.080, "COARSE 17.6 pts/osc")
wb = onePair(0.030, "FINER  47 pts/osc")

println("\n================  SUMMARY  ================")
@printf("coarse pair:  relative = %.4e   converged = %s\n", wa.relative, string(wa.converged))
@printf("fine   pair:  relative = %.4e   converged = %s\n", wb.relative, string(wb.converged))
println("The two must differ by orders, or the verdict is not usable.")
