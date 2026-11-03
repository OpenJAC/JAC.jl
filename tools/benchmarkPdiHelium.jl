#
#  THE ONLY EXTERNAL ANCHOR PhotoDoubleIonization HAS EVER HAD, re-measured.
#
#  Helium's double-to-single photoionization ratio sigma(2+)/sigma(+) is measured, and it is the one number
#  against which this module's absolute scale has ever been checked:  on 15-Sep-2026 it came out 2.70 / 3.18 /
#  2.28 % at 200 / 400 / 800 eV against 3-4 % measured near 200 eV falling towards an asymptotic 1.68 %.  That
#  measurement is quoted in the module's own STATUS block and is the sole evidence that the scale is no longer
#  "three orders of magnitude too small".
#
#  FIVE REPAIRS HAVE LANDED SINCE, TWO OF THEM SEVERE, and the anchor has not been re-measured once:
#    01-Oct  a pole sitting exactly ON the upper quadrature endpoint, where no principal value exists;
#    05-Oct  the same pole, properly excluded;
#    01-Nov  the Green space was the WRONG ELECTRON NUMBER -- two-electron cores against a two-electron
#            initial state, so every junction of the second-order sum paired unequal electron counts;
#    02-Nov  the Green orbitals were read on the WRONG RADIAL MESH, with norms wrong by factors of 3 to 9;
#    02-Nov  all three multiplets can now stand on ONE orbital set, so the dropped -E <r|s> overlap term of
#            the H-E identity is no longer carried at one-orbital junctions.
#
#  So this script answers ONE question:  does the 15-Sep helium ratio survive those five repairs, and what does
#  the shared orbital set do to it?  Three runs, one variable between them:
#
#    sigma(+)     single photoionization, needed for the denominator
#    sigma(2+)    SEPARATE orbital sets  -- the 15-Sep configuration
#    sigma(2+)    ONE shared orbital set -- what 02-Nov made possible
#
#  Run as:  julia --project=. tools/benchmarkPdiHelium.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")

##  THE BOX IS MATCHED TO THE GREEN SPACE, NOT TO THE INITIAL STATE (Rule 12).  Neutral He's 1s^2 sees
##  Z_eff ~ 1.69 and needs only r ~ 3, but the (N-1)-electron cores carry He+ 2s, whose outer turning point is
##  r_+ = (n^2/Z_eff)(1 + sqrt(1 - l(l+1)/n^2)) = (4/2)(1+1) = 4 a.u., so 2.5 r_+ = 10.  A box of 12 covers it
##  with margin;  a box chosen for the initial state alone would misrepresent the Green space's 2s and 2p.
##  SETTABLE, because the first run of 02-Nov-2026 got it wrong:  12.0 a.u. was sized for He+ cores (Z_eff = 2),
##  but the SCF over [1s^2, 1s 2s, 1s 2p] returns NEUTRAL-He Rydberg orbitals, which see Z_eff ~ 1 and need
##  2.5 r_+ ~ 20 a.u.  The symptom was quiet:  all four norms read 1.0000000000 and both precondition checks
##  behaved, while 2s came out at -0.065 Ha where even neutral He's 1s2s lies near -0.15 Ha.
const RBOX   = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 24.0
const grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = RBOX)
const nm     = Nuclear.Model(2.0)
const OMEGAS = [200.0, 400.0, 800.0]
const NSHARE = 9
const MAXKAP = 2

println("="^118)
println("PhotoDoubleIonization: the helium anchor re-measured on today's code")
@printf("box %.1f a.u., %d grid points;  %d sharings, maxKappa = %d;  omega = %s eV\n",
        RBOX, grid.NoPoints, NSHARE, MAXKAP, string(OMEGAS))
println("="^118)

##  ----------------------------------------------------------------------------------------------------------
##  The shared orbital set, over a list WIDER than any multiplet's CI so that it spans the Green space
##  ----------------------------------------------------------------------------------------------------------
spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")]
wa          = Representation("He span", nm, grid, spanConfigs, MeanFieldMultiplet(MeanFieldSettings()) )
wb          = redirect_stdout(devnull) do
    generate(wa, output=true)
end
orbitals    = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
println("\n>> shared orbital set: ", join(sort(string.(collect(keys(orbitals)))), ", "))
for sh in sort(collect(keys(orbitals)), by=string)
    o = orbitals[sh]
    nrm = RadialIntegrals.overlap(o, o, grid)
    @printf("   %-8s  e = %12.6f Ha   <a|a> = %14.10f   length(P) = %5d of %5d\n",
            string(sh), o.energy, nrm, length(o.P), grid.NoPoints)
end

##  the (N-1)-electron ionic cores, as a CI on that SAME set
gMultiplet = redirect_stdout(devnull) do
    Hamiltonian.performCIwithFrozenOrbitals([Configuration("1s"), Configuration("2s"), Configuration("2p")],
                                             orbitals, nm, grid, AsfSettings(); printout=false)
end
@printf(">> Green space: %d levels, %d electron(s)\n", length(gMultiplet.levels),
        gMultiplet.levels[1].basis.NoElectrons)

##  ----------------------------------------------------------------------------------------------------------
##  sigma(+)  -- single photoionization, the denominator of the measured ratio
##  ----------------------------------------------------------------------------------------------------------
println("\n" * "-"^118)
println("RUN 1 of 3:  sigma(+), single photoionization of He")
println("-"^118)
t1 = time()
piSettings = PhotoIonization.Settings(PhotoIonization.Settings();
                  multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=OMEGAS,
                  calcPartialCs=false, printBefore=false, lineSelection=LineSelection(true, indexPairs=[(1,1)]) )
compSingle = Atomic.Computation(Atomic.Computation(), name="He sigma+", grid=grid, nuclearModel=nm,
                  initialConfigs  = [Configuration("1s^2")],
                  finalConfigs    = [Configuration("1s")],
                  processSettings = piSettings )
resSingle  = perform(compSingle; output=true)
@printf(">> sigma(+) done in %.1f s\n", time()-t1)

##  ----------------------------------------------------------------------------------------------------------
##  sigma(2+)  -- twice, with the orbital source as the ONE variable
##  ----------------------------------------------------------------------------------------------------------
function doubleRun(label, source)
    println("\n" * "-"^118)
    println(label)
    println("-"^118)
    t = time()
    dSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                     multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=OMEGAS,
                     NoEnergySharings=NSHARE, maxKappa=MAXKAP, calcDifferentialCs=false, printBefore=false,
                     lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                     eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )
    comp = Atomic.Computation(Atomic.Computation(), name=label, grid=grid, nuclearModel=nm,
                     initialConfigs  = [Configuration("1s^2")],
                     finalConfigs    = [Configuration("1s^0")],
                     orbitalSource   = source,
                     processSettings = dSettings )
    res  = perform(comp; output=true)
    @printf(">> %s done in %.1f s\n", label, time()-t)
    return( res )
end

resSepar  = doubleRun("RUN 2 of 3:  sigma(2+), SEPARATE orbital sets (the 15-Sep-2026 configuration)",
                      Basics.SeparateOrbitals())
resShared = doubleRun("RUN 3 of 3:  sigma(2+), ONE shared orbital set (possible since 02-Nov-2026)",
                      Basics.GivenOrbitals(orbitals))

println("\n" * "="^118)
println("READ THE THREE TABLES ABOVE FOR THE CROSS SECTIONS AND FORM THE RATIO sigma(2+)/sigma(+) BY HAND.")
println("Against 15-Sep-2026, COULOMB:  2.70 % at 200 eV, 3.18 % at 400 eV, 2.28 % at 800 eV,")
println("and against experiment:  3-4 % near 200 eV, falling towards an asymptotic 1.68 %.")
println("="^118)
