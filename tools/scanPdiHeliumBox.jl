#
#  WHY THE HELIUM DOUBLE-IONIZATION CROSS SECTION STILL MOVES WITH THE RADIAL BOX, 02-Nov-2026.
#
#  Measured that day on He at 200 / 400 / 800 eV:  doubling the box from 12 to 24 a.u. left SINGLE
#  photoionization unchanged to 0.3 %, while DOUBLE photoionization moved by -19 %, -25 % and +49 % with
#  separate orbital sets and by +3 %, -18 % and +23 % with one shared set.  A cross section may not depend
#  on the radius of an imaginary sphere, so the calculation is not converged, and no statement about the
#  normalization can be made until it is.
#
#  THE HYPOTHESIS, stated before the run so that it can be refuted:  it is the GREEN SPACE's orbitals that
#  feel the wall, not the machinery.  The intermediate state of this process is a helium ION, but both runs
#  of 02-Nov built the Green space from a NEUTRAL-atom average field, where an outer electron sees the
#  nucleus almost fully screened and 2s/2p come out at -0.072 / -0.036 Ha instead of -0.5.  Such orbitals are
#  bound by the box rather than by the atom.  Ion orbitals are not:  measured the same day, an SCF on the
#  He+ configurations gives 1s -2.00011, 2s -0.50003, 2p -0.50003 Ha, IDENTICAL to five decimals at 12 and at
#  24 a.u.
#
#  So:  ONE VARIABLE, the Green space's orbitals, across five boxes.  The initial state is its own SCF in
#  BOTH arms, so nothing else differs.
#
#    arm ION      the Green space is its own SCF on the He+ configurations -- compact, box-independent
#    arm NEUTRAL  the Green space is a frozen-orbital CI on the neutral span set -- diffuse
#
#  PREDICTION:  arm ION is flat in the box, arm NEUTRAL drifts.  If BOTH drift, the orbitals are not the
#  cause and the number of intermediate states is the next suspect.
#
#  One photon energy only, 200 eV, because that is where Kornberg & Miraglia's Table II has a TABULATED
#  entry rather than something to interpolate:  sigma+ = 5.44e-02 Mb and sigma2+ (C2) = 1.69e-03 Mb, i.e. a
#  target ratio of 3.11 %.  Dropping from three energies to one pays for the extra boxes.
#
#  Run as:  julia --project=. tools/scanPdiHeliumBox.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")

const BOXES  = [12.0, 18.0, 24.0, 32.0, 40.0]
const OMEGA  = 200.0
const NSHARE = 9
const MAXKAP = 2
const nm     = Nuclear.Model(2.0)
const ionConfigs  = [Configuration("1s"), Configuration("2s"), Configuration("2p")]
const spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")]

println("="^118)
println("PhotoDoubleIonization, He at $(OMEGA) eV:  is the box dependence the GREEN SPACE's orbitals?")
println("boxes $(BOXES) a.u.;  $(NSHARE) sharings, maxKappa = $(MAXKAP);  target sigma2+/sigma+ = 3.11 % (K&M C2)")
println("="^118)

results = []

for rbox in BOXES
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = rbox)
    @printf("\n%s\nBOX %.1f a.u., %d grid points\n%s\n", "-"^118, rbox, grid.NoPoints, "-"^118)

    ##  sigma(+), the denominator;  known box-stable but recomputed so the ratio is self-consistent.
    ##  The tables are left to PRINT rather than captured:  `perform` exposes no ResultKey for
    ##  PhotoIonization or PhotoDoubleIonization lines, so the log is the record and is parsed afterwards.
    t0 = time()
    piSet = PhotoIonization.Settings(PhotoIonization.Settings();
                multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                calcPartialCs=false, printBefore=false, lineSelection=LineSelection(true, indexPairs=[(1,1)]) )
    println("  >>>> MARK  box=$(rbox)  arm=SINGLE")
    perform(Atomic.Computation(Atomic.Computation(), name="s+", grid=grid, nuclearModel=nm,
                initialConfigs=[Configuration("1s^2")], finalConfigs=[Configuration("1s")],
                processSettings=piSet))
    @printf("  sigma(+) computed in %.1f s\n", time()-t0)

    ##  the two Green spaces
    ionMp = redirect_stdout(devnull) do
        SelfConsistent.performSCF(ionConfigs, nm, grid, AsfSettings(); printout=false)
    end
    spanOrbs = redirect_stdout(devnull) do
        SelfConsistent.performSCF(spanConfigs, nm, grid, AsfSettings(); printout=false)
    end.levels[1].basis.orbitals
    neuMp = redirect_stdout(devnull) do
        Hamiltonian.performCIwithFrozenOrbitals(ionConfigs, spanOrbs, nm, grid, AsfSettings(); printout=false)
    end
    io = ionMp.levels[1].basis.orbitals
    @printf("  arm ION      2s %9.5f  2p %9.5f Ha\n", io[Subshell("2s_1/2")].energy, io[Subshell("2p_1/2")].energy)
    @printf("  arm NEUTRAL  2s %9.5f  2p %9.5f Ha\n",
            spanOrbs[Subshell("2s_1/2")].energy, spanOrbs[Subshell("2p_1/2")].energy)

    for (arm, gMp) in (("ION", ionMp), ("NEUTRAL", neuMp))
        t = time()
        dSet = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                   multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                   NoEnergySharings=NSHARE, maxKappa=MAXKAP, calcDifferentialCs=false, printBefore=false,
                   lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                   eeInteraction=CoulombInteraction(), gMultiplet=gMp )
        println("  >>>> MARK  box=$(rbox)  arm=$(arm)")
        perform(Atomic.Computation(Atomic.Computation(), name="s2+ $arm", grid=grid, nuclearModel=nm,
                    initialConfigs=[Configuration("1s^2")], finalConfigs=[Configuration("1s^0")],
                    orbitalSource=Basics.SeparateOrbitals(), processSettings=dSet))
        @printf("  arm %-8s sigma(2+) computed in %6.1f s\n", arm, time()-t)
        flush(stdout)
    end
end
println("\n" * "="^118)
println("DONE.  Read the per-run tables above.")
println("="^118)
