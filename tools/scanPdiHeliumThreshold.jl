#
#  IS THE MODEL RIGHT WHERE ITS OWN MECHANISM IS THE DOMINANT ONE?
#
#  `module-PhotoDoubleIonization.jl` computes BOTH TIME ORDERINGS OF THE KNOCK-OUT (TS1) MECHANISM and nothing
#  else:  no shake-off, and a single-configuration 1s^2 initial state carries no ground-state correlation.  The
#  literature splits helium's double photoionization into exactly those two mechanisms, and they dominate in
#  OPPOSITE regimes -- knock-out just above the 79 eV threshold, shake-off at high energy, where it alone sets
#  the asymptotic ratio (1.645-1.67 % theory; 1.72 +- 0.12 % measured at 7 keV).
#
#  SO EVERY NUMBER THIS MODULE HAS EVER BEEN JUDGED ON -- 200, 400 AND 800 eV -- LIES AT OR BEYOND THE EDGE OF
#  ITS OWN DOMAIN, and the one region where it is physically complete has never been run.  This script runs it.
#
#  THE TEST IS THE SHAPE, NOT THE TOTAL, which is what makes it a test of the MECHANISM rather than of the
#  prefactor.  The two mechanisms write different energy-sharing distributions:
#
#    * KNOCK-OUT shares the excess energy almost evenly, because the second electron is ejected in a
#      half-collision.  Near threshold the distribution is close to FLAT -- measured flat at ~20 eV excess.
#    * SHAKE-OFF is sudden and strongly asymmetric:  one fast electron and one slow one, so it piles weight at
#      the two ENDS and deepens the U.
#
#  Measured at omega = 200 eV (excess 122 eV), this module gives dsigma/d eps1 = 117.2 barn/eV at the outermost
#  node against 18.84 at the midpoint, a U-ratio of 6.22.  (Those two figures were first read off the module's
#  own table as 3189.9 and 512.7, which that table labels [barn]/[eV] while carrying barn per HARTREE; the
#  U-ratio is a ratio and so was never affected.)  THE PREDICTION IS THEREFORE FALSIFIABLE AND
#  DIRECTIONAL:  that ratio must FALL towards 1 as omega approaches threshold.  If it does, the module reproduces
#  the knock-out shape where knock-out is what the physics does, and its validity domain is established from
#  below rather than asserted.  If the U stays deep at 100 eV, the shape is not knock-out at all and the
#  agreement at 200 eV would be a coincidence.
#
#  120, 150 and 200 eV are the three energies at which Kornberg & Miraglia's Fig. 2 PLOTS d sigma / d eps1 in
#  Mb/eV.  That is a figure to be read off, not a table, and it has NOT been read here -- so this scan tests the
#  SHAPE against the knock-out expectation, not digit-for-digit against them.  Their Table II at 200 eV is the
#  tabulated entry, and tools/comparePdiHeliumKornberg.jl is where that comparison is made.
#
#  AND K&M ARE NOT A NEAR-THRESHOLD BENCHMARK, by their own account:  their abstract says the velocity form agrees
#  well with measurement AT HIGH photon energies and that at LOW photon energies the approximation "fails to
#  account for the data in magnitude and shape".  Their 3C final state is a high-energy construction, whereas a
#  knock-out-only amplitude is at its best near threshold.  The two models are good in OPPOSITE regimes, so
#  disagreement at 100-150 eV would say little about either.
#
#  TWO CAVEATS THAT THIS SCRIPT DOES NOT HIDE.  Near threshold the shared excess is small, so the OUTERMOST
#  Gauss-Legendre node is a very SLOW electron -- 0.35 eV at omega = 100 eV with nine sharings, a de Broglie
#  wavelength near 39 a.u.  A 40 a.u. box then holds about one wavelength, and a Galerkin orbital that has not
#  reached its asymptotic form inside the box would be fitted against a form it does not yet have.  That is
#  exactly what `Continuum.asymptoticAmplitude`'s residual tiers exist to catch, so THE GUARD WARNINGS ARE PART
#  OF THE RESULT and are reported rather than filtered.  And the total now carries the identical-particle 1/2
#  supplied on 03-Nov-2026, so these totals are NOT comparable with any recorded before that date.
#
#  Run as:  julia --project=. tools/scanPdiHeliumThreshold.jl  [rbox]
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Coulomb functions")
setDefaults("unit: rate", "a.u.")

const RBOX    = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 40.0
const grid    = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = RBOX)
const nm      = Nuclear.Model(2.0)
const OMEGAS  = [100.0, 120.0, 150.0, 200.0]
const MAXKAP  = 2
const NINT    = 9           ##  held fixed, as in tools/scanPdiHeliumSharings.jl
const NSHARE  = 9           ##  converged to 0.6 % as of 03-Nov-2026; nine nodes also resolve the shape

println("="^128)
println("PhotoDoubleIonization: does the sharing distribution FLATTEN towards threshold, as knock-out requires?")
@printf("box %.1f a.u., %d grid points;  maxKappa = %d;  NoEnergySharings = %d;  NoIntermediateEnergies = %d;  omega = %s eV\n",
        RBOX, grid.NoPoints, MAXKAP, NSHARE, NINT, string(OMEGAS))
println("="^128)

##  ----------------------------------------------------------------------------------------------------------
##  The shared orbital set and the (N-1)-electron Green space, exactly as tools/benchmarkPdiHelium.jl builds
##  them:  one orbital set for all three multiplets, which is the precondition of any explicit second-order
##  amplitude and is why the span list is WIDER than any single multiplet's CI.
##  ----------------------------------------------------------------------------------------------------------
spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")]
wa          = Representation("He span", nm, grid, spanConfigs, MeanFieldMultiplet(MeanFieldSettings()) )
wb          = redirect_stdout(devnull) do
    generate(wa, output=true)
end
orbitals    = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
println("\n>> shared orbital set: ", join(sort(string.(collect(keys(orbitals)))), ", "))

gMultiplet = redirect_stdout(devnull) do
    Hamiltonian.performCIwithFrozenOrbitals([Configuration("1s"), Configuration("2s"), Configuration("2p")],
                                             orbitals, nm, grid, AsfSettings(); printout=false)
end
@printf(">> Green space: %d levels, %d electron(s)\n\n", length(gMultiplet.levels),
        gMultiplet.levels[1].basis.NoElectrons)

##  ----------------------------------------------------------------------------------------------------------
##  One run per photon energy
##  ----------------------------------------------------------------------------------------------------------
results = Dict{Float64,Any}()
for (i, omega) in enumerate(OMEGAS)
    println("-"^128)
    @printf("RUN %d of %d:  omega = %.1f eV\n", i, length(OMEGAS), omega)
    println("-"^128)
    t = time()
    dSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                     multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[omega],
                     NoEnergySharings=NSHARE, NoIntermediateEnergies=NINT, maxKappa=MAXKAP,
                     calcDifferentialCs=true, printBefore=false,
                     lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                     eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )
    comp = Atomic.Computation(Atomic.Computation(), name="He sigma2+ omega=$omega", grid=grid, nuclearModel=nm,
                     initialConfigs  = [Configuration("1s^2")],
                     finalConfigs    = [Configuration("1s^0")],
                     orbitalSource   = Basics.GivenOrbitals(orbitals),
                     processSettings = dSettings )
    res = perform(comp; output=true)
    ##  A TYPED result key, not the retired string one:  it is checked where it is WRITTEN (priority item 39).
    results[omega] = res[PhotoDoubleIonization.Settings]
    @printf(">> omega = %.1f eV done in %.1f s\n\n", omega, time()-t)
end

##  ----------------------------------------------------------------------------------------------------------
##  THE SHAPE TABLE.  The U-ratio is the outermost node against the midpoint one, both in the Coulomb gauge.
##  Collected into arrays FIRST:  a top-level `for` body is its own soft scope, so assigning a running value
##  inside it makes that name local and reading it back raises an UndefVarError rather than carrying.
##  ----------------------------------------------------------------------------------------------------------
tob(x) = Defaults.convertUnits("cross section: from atomic to barn", x)

shs    = [results[om][1].sharings                                     for om in OMEGAS]
mid    = [div(length(s), 2) + 1                                       for s in shs]
dEdge  = [tob(s[1].differentialCs.Coulomb)   / 27.211386              for s in shs]
dMid   = [tob(shs[k][mid[k]].differentialCs.Coulomb) / 27.211386      for k in eachindex(OMEGAS)]
eEdge  = [s[1].epsilon1  * 27.211386                                  for s in shs]
total  = [tob(results[om][1].crossSection.Coulomb)                    for om in OMEGAS]
totalB = [tob(results[om][1].crossSection.Babushkin)                  for om in OMEGAS]

println("="^128)
println("THE SHAPE:  a U-ratio falling towards 1 is the knock-out signature; a deep U near threshold refutes it")
println("="^128)
@printf("%10s  %10s  %12s  %14s  %14s  %9s   %13s  %13s\n",
        "omega[eV]", "excess[eV]", "eps1 edge[eV]", "dsig edge", "dsig mid", "U-ratio", "sigma2+ Cou", "sigma2+ Bab")
@printf("%10s  %10s  %12s  %14s  %14s  %9s   %13s  %13s\n",
        "", "", "", "[barn/eV]", "[barn/eV]", "", "[barn]", "[barn]")
println("-"^128)
for k in eachindex(OMEGAS)
    @printf("%10.1f  %10.2f  %12.4f  %14.2f  %14.2f  %9.3f   %13.4f  %13.4f\n",
            OMEGAS[k], sum(x -> x.epsilon1 + x.epsilon2, shs[k][1:1]) * 27.211386, eEdge[k],
            dEdge[k], dMid[k], dEdge[k]/dMid[k], total[k], totalB[k])
end
println("="^128)
println(">> U-ratio at 200 eV before this scan, 24 a.u. box:  6.22  (117.2 against 18.84 barn/eV)")
println(">> NOTE:  the module's OWN energy-sharing table labels this column [barn]/[eV] and prints barn per")
println(">>        HARTREE -- a factor 27.2114.  Proven unit-blind by its quadrature identity:  the weights of")
println(">>        Basics.determineEnergySharings sum to the excess in a.u. (4.53028 against 4.53046 at 200 eV),")
println(">>        and 0.5*sum(w*d) with d AS DISPLAYED reproduces the printed total to all 8 digits.  So any")
println(">>        comparison with Kornberg & Miraglia Fig. 2 must divide that column by 27.2114 first.")
println(">> Totals carry the identical-particle 1/2 of 03-Nov-2026 and do NOT compare with earlier records.")
