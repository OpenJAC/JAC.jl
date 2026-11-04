#
#  DOES THE SHARING INTEGRAL CONVERGE, NOW THAT NEITHER CONFOUND IS PRESENT?
#
#  `module-PhotoDoubleIonization.jl` names the sharing quadrature as its LEADING error:  "Measured at 200 eV:
#  the total moves +2.3 %, +12.4 % and +14.1 % as the Gauss-Legendre points go 5 -> 9 -> 15 -> 21, i.e. +29 %
#  over the whole scan and still rising."  That measurement is of 15-Sep-2026 and it varied TWO things at once,
#  which is why it is worth repeating rather than believing:
#
#    1.  THE INTERMEDIATE MESH WAS ENTANGLED WITH THE SHARING MESH until 30-Sep-2026 -- `NoIntermediateEnergies`
#        did not exist, so asking for more sharings silently made the second-order sum more complete as well.
#        The module's own docstring records a factor 2.7 at FIXED sharing between 3 and 41 sharings from this
#        alone.  This scan fixes `NoIntermediateEnergies` at a constant value, so the sharing count is the only
#        variable.
#
#    2.  THE ASYMPTOTIC AMPLITUDE WAS READ AS A MAXIMUM OVER THE LAST 300 GRID POINTS until 03-Nov-2026
#        (`6d1d5700`).  That is a fixed POINT COUNT:  at hp = 1e-2 it spans one half-wavelength only above
#        ~14 eV, so below that the window need contain no crest and the maximum samples the outgoing wave at an
#        arbitrary phase -- measured against JAC's own analytic wave, exact at 0.272 eV, 1.252 at 0.544 eV and
#        2.24 at 0.136 eV.  A cross section goes as the SQUARE of the normalization.
#
#  AND (2) BEARS ON A SHARING SCAN IN PARTICULAR, which is the reason for this script.  `Basics.determineEnergy-
#  Sharings` is plain Gauss-Legendre mapped onto [0, excess], so adding points pushes the OUTERMOST node towards
#  the endpoints -- towards a very slow electron.  For He at 200 eV (excess 121 eV) the lowest node marches
#
#        n =  5   ->   5.676 eV        n = 15   ->   0.726 eV
#        n =  9   ->   1.926 eV        n = 21   ->   0.378 eV
#
#  straight into the band where the old normalization was worst.  At 200 eV the sharing distribution is U-shaped,
#  so those outermost nodes carry the LARGEST differential cross sections.  A quadrature whose added points land
#  on progressively more over-normalized orbitals rises with the point count and never settles:  exactly the
#  reported signature.  THE PREDICTION IS THEREFORE FALSIFIABLE AND DIRECTIONAL -- the old normalization was too
#  LARGE at low energy, so the drift should shrink AND the total should come DOWN.
#
#  THE PER-NODE DIFFERENTIAL CROSS SECTIONS ARE PRINTED, NOT ONLY THE TOTALS, because they are what tests the
#  MECHANISM rather than the outcome:  if the above is right, the nodes shared between two meshes agree while
#  the newly added low-energy ones are where any remaining movement sits.
#
#  Run as:  julia --project=. tools/scanPdiHeliumSharings.jl  [rbox]
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Coulomb functions")
setDefaults("unit: rate", "a.u.")

##  24 a.u. is box-converged for this observable as of 03-Nov-2026:  sigma(2+) moved 0.05 % over 24 -> 40 a.u.
##  (2347.98 vs 2346.84 barn) once the normalization was fitted rather than windowed, against 43 % before.
const RBOX   = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 24.0
const grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = RBOX)
const nm     = Nuclear.Model(2.0)
const OMEGA  = 200.0
const MAXKAP = 2
const NINT   = 9            ##  FIXED, so that the sharing count is the only variable
const NSHARES = [5, 9, 15, 21]

println("="^128)
println("PhotoDoubleIonization: does the sharing integral converge with the intermediate mesh HELD FIXED?")
@printf("omega = %.0f eV;  box %.1f a.u., %d grid points;  maxKappa = %d;  NoIntermediateEnergies = %d FIXED;  sharings = %s\n",
        OMEGA, RBOX, grid.NoPoints, MAXKAP, NINT, string(NSHARES))
println("Baseline to beat, 15-Sep-2026 (BOTH confounds present):  +2.3 %, +12.4 %, +14.1 % step to step, +29 % overall.")
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
##  One run per sharing count
##  ----------------------------------------------------------------------------------------------------------
results = Dict{Int64,Any}()
for nShare in NSHARES
    println("-"^128)
    @printf("RUN %d of %d:  NoEnergySharings = %d  (NoIntermediateEnergies = %d, fixed)\n",
            findfirst(isequal(nShare), NSHARES), length(NSHARES), nShare, NINT)
    println("-"^128)
    t = time()
    dSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                     multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                     NoEnergySharings=nShare, NoIntermediateEnergies=NINT, maxKappa=MAXKAP,
                     calcDifferentialCs=true, printBefore=false,
                     lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                     eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )
    comp = Atomic.Computation(Atomic.Computation(), name="He sigma2+ nShare=$nShare", grid=grid, nuclearModel=nm,
                     initialConfigs  = [Configuration("1s^2")],
                     finalConfigs    = [Configuration("1s^0")],
                     orbitalSource   = Basics.GivenOrbitals(orbitals),
                     processSettings = dSettings )
    res = perform(comp; output=true)
    ##  A TYPED result key, not the retired string one:  it is checked where it is WRITTEN (priority item 39).
    results[nShare] = res[PhotoDoubleIonization.Settings]
    @printf(">> nShare = %d done in %.1f s\n\n", nShare, time()-t)
end

##  ----------------------------------------------------------------------------------------------------------
##  THE TWO TABLES.  Totals first, then the per-node differential cross sections that test the mechanism.
##  ----------------------------------------------------------------------------------------------------------
tob(x) = Defaults.convertUnits("cross section: from atomic to barn", x)

println("="^128)
println("TABLE 1 -- THE TOTAL, which is what the +29 % claim is about")
println("="^128)
@printf("%8s   %18s %10s   %18s %10s\n", "nShare", "Coulomb [barn]", "step", "Babushkin [barn]", "step")
##  Collected into arrays FIRST:  a top-level `for` body is its own soft scope, so assigning a running
##  total inside it makes that name local and reading it back raises an UndefVarError rather than carrying.
csC = [tob(results[n][1].crossSection.Coulomb)    for n in NSHARES]
csB = [tob(results[n][1].crossSection.Babushkin)  for n in NSHARES]
for (i, nShare) in enumerate(NSHARES)
    sC = i == 1 ? "--" : @sprintf("%+8.2f %%", 100*(csC[i]-csC[i-1])/csC[i-1])
    sB = i == 1 ? "--" : @sprintf("%+8.2f %%", 100*(csB[i]-csB[i-1])/csB[i-1])
    @printf("%8d   %18.4f %10s   %18.4f %10s\n", nShare, csC[i], sC, csB[i], sB)
end
@printf("\nOVERALL, %d -> %d:  Coulomb %+.2f %%,  Babushkin %+.2f %%      (15-Sep-2026, both confounds present:  +29 %%)\n",
        NSHARES[1], NSHARES[end], 100*(csC[end]-csC[1])/csC[1], 100*(csB[end]-csB[1])/csB[1])

println("\n" * "="^128)
println("TABLE 2 -- PER-NODE DIFFERENTIAL CROSS SECTIONS, the mechanism test")
println("The lowest node marches towards zero energy as points are added;  that is where the old 300-point")
println("normalization window was worst, and so where any surviving movement must sit.")
println("="^128)
for nShare in NSHARES
    line = results[nShare][1]
    @printf("\n  nShare = %d\n", nShare)
    @printf("  %4s  %12s %12s  %14s  %18s\n", "node", "eps1 [eV]", "eps2 [eV]", "weight", "d(sigma)/d(eps) C")
    for (i,sh) in enumerate(line.sharings)
        @printf("  %4d  %12.4f %12.4f  %14.6e  %18.6e\n", i,
                Defaults.convertUnits("energy: from atomic to eV", sh.epsilon1),
                Defaults.convertUnits("energy: from atomic to eV", sh.epsilon2),
                sh.weight, sh.differentialCs.Coulomb)
    end
end
println("\n" * "="^128)
