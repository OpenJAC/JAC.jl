#
#  IS THE 3.868 % RATIO CONVERGED?  The two truncations that have never been scanned together.
#
#  `tools/comparePdiHeliumKornberg.jl` put sigma(2+)/sigma(+) at 3.868 % for He at 200 eV, INSIDE the
#  measured 3-4 % band.  That was obtained at maxKappa = 2 and NoIntermediateEnergies = 9, and NEITHER
#  is converged, so the agreement is an order-and-sign statement and not yet an accuracy figure.  This
#  scan is what would turn it into one -- or show that it was luck.
#
#  TWO ARMS, ONE VARIABLE EACH, around the common point (maxKappa = 2, NINT = 9) already measured:
#
#    arm KAPPA   maxKappa = 2, 3, 4          at NINT = 9      ... the partial-wave pairs of the amplitude
#    arm INT     NINT     = 9, 17, 25        at maxKappa = 2  ... the second-order intermediate sum
#
#  THE BOX AND THE SHARING COUNT ARE HELD FIXED at 40 a.u. and 9, because both were scanned already:
#  sigma(2+) moves 0.05 % over 24 -> 40 a.u. once the asymptotic amplitude is fitted, and +0.58 % over
#  5 -> 21 sharings once the intermediate mesh is held fixed.  Varying four things would make a drift
#  unattributable, which is exactly how the withdrawn "+29 %" claim was produced.
#
#  sigma(+) IS NOT RECOMPUTED:  it does not depend on either variable.  Measured 04-Nov-2026 on this box
#  and these orbitals, 5.9499e-02 Mb (Coulomb).  The ratio is formed against that.
#
#  EACH RUN PRINTS ITS OWN RESULT LINE AS IT FINISHES rather than only a table at the end, so a job that
#  dies halfway still leaves usable evidence -- the runs are ~40 min each and maxKappa = 4 may be far
#  longer.  The order below deliberately alternates the two arms so that a partial run covers both.
#
#  Run as:  julia --project=. tools/scanPdiHeliumConvergence.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Coulomb functions")

const RBOX    = 40.0
const grid    = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = RBOX)
const nm      = Nuclear.Model(2.0)
const OMEGA   = 200.0
const NSHARE  = 9
const SIG1P   = 5.9499e-02        ##  Mb, Coulomb, same orbitals/grid/prefactor (04-Nov-2026)

##  (label, maxKappa, NoIntermediateEnergies);  the (2,9) point is the one already measured
const RUNS = [ ("kappa 3", 3,  9), ("int 17",  2, 17), ("kappa 4", 4,  9), ("int 25",  2, 25) ]

println("="^128)
@printf("PhotoDoubleIonization convergence: He at %.0f eV, box %.1f a.u. (%d points), %d sharings FIXED\n",
        OMEGA, RBOX, grid.NoPoints, NSHARE)
println("Reference point already measured: maxKappa = 2, NINT = 9  ->  2301.4539 barn Coulomb, ratio 3.868 %")
println("="^128)

spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")]
wa          = Representation("He span", nm, grid, spanConfigs, MeanFieldMultiplet(MeanFieldSettings()) )
wb          = redirect_stdout(devnull) do
    generate(wa, output=true)
end
orbitals    = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
gMultiplet  = redirect_stdout(devnull) do
    Hamiltonian.performCIwithFrozenOrbitals([Configuration("1s"), Configuration("2s"), Configuration("2p")],
                                             orbitals, nm, grid, AsfSettings(); printout=false)
end
@printf("\n>> shared orbital set: %s;  Green space %d levels\n\n",
        join(sort(string.(collect(keys(orbitals)))), ", "), length(gMultiplet.levels))

tob(x) = Defaults.convertUnits("cross section: from atomic to barn", x)
results = Tuple{String,Int64,Int64,Float64,Float64,Float64}[]

for (label, mk, ni) in RUNS
    println("-"^128)
    @printf("RUN %s:  maxKappa = %d, NoIntermediateEnergies = %d\n", label, mk, ni)
    println("-"^128)
    t = time()
    dSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                     multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                     NoEnergySharings=NSHARE, NoIntermediateEnergies=ni, maxKappa=mk,
                     calcDifferentialCs=true, printBefore=false,
                     lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                     eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )
    comp = Atomic.Computation(Atomic.Computation(), name="He $label", grid=grid, nuclearModel=nm,
                     initialConfigs  = [Configuration("1s^2")],
                     finalConfigs    = [Configuration("1s^0")],
                     orbitalSource   = Basics.GivenOrbitals(orbitals),
                     processSettings = dSettings )
    res  = perform(comp; output=true)
    line = res[PhotoDoubleIonization.Settings][1]
    csC  = tob(line.crossSection.Coulomb);   csB = tob(line.crossSection.Babushkin)
    ratio = 100 * csC * 1.0e-6 / SIG1P
    push!(results, (label, mk, ni, csC, csB, ratio))
    ##  PRINTED NOW, not only in the final table:  a job that dies here still leaves this line behind.
    @printf("\n>>> RESULT  %-8s  maxKappa=%d  NINT=%2d   sigma2+ = %10.4f barn (Cou), %10.4f (Bab)   ratio = %6.3f %%   [%.0f s]\n\n",
            label, mk, ni, csC, csB, ratio, time()-t)
end

println("\n" * "="^128)
println("CONVERGENCE -- read each arm against the common (maxKappa 2, NINT 9) point: 2301.4539 barn, 3.868 %")
println("="^128)
@printf("%-10s %9s %6s %16s %16s %10s %12s\n", "run", "maxKappa", "NINT", "sigma2+ Cou", "sigma2+ Bab", "ratio", "vs (2,9)")
println("-"^128)
@printf("%-10s %9d %6d %16.4f %16.4f %9.3f %% %11s\n", "reference", 2, 9, 2301.4539, 4234.7389, 3.868, "--")
for (label, mk, ni, csC, csB, ratio) in results
    @printf("%-10s %9d %6d %16.4f %16.4f %9.3f %% %10.2f %%\n",
            label, mk, ni, csC, csB, ratio, 100*(csC-2301.4539)/2301.4539)
end
println("="^128)
println(">> MEASURED near 200 eV: 3-4 %.  A ratio that stays in the band across both arms makes 3.868 % an")
println(">> accuracy statement;  one that walks out of it makes the agreement at (2,9) a coincidence.")
