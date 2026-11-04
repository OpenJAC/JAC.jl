#
#  THE ONE COMPARISON THIS MODULE HAS A TABULATED BENCHMARK FOR:  He at omega = 200 eV.
#
#  Kornberg & Miraglia, Phys. Rev. A 48, 3714 (1993), Table II, at E_gamma = 200 eV:
#      sigma(+) = 5.44e-02 Mb,   sigma(2+)_C2 = 1.69e-03 Mb,   sigma(2+)_C3 = 1.3e-02 Mb
#
#  TABLE II HAS TWO DOUBLE-IONIZATION COLUMNS AND THE SMALLER ONE IS NOT THEIR MODEL -- read 04-Nov-2026.
#  C3 is the correlated three-Coulomb-wave final state the paper is ABOUT;  C2 is the independent-electron
#  product of two Coulomb waves.  They differ by 7.7x at 200 eV.  Earlier notes in this repository quoted
#  1.69e-03 as "the K&M value" without saying it was C2, so both are printed here.
#
#  AND NEITHER IS A BENCHMARK AT 200 eV -- the paper says so itself.  "At intermediate energies, the
#  calculated ratio is up to five times larger than the data", and C3 gives 23.9 % here against a measured
#  3-4 %.  At the other end, "at E_gamma = 84 eV the theory gives a ratio of 0.004 % while the data is
#  0.6 %".  Their velocity form is trustworthy only ABOVE about 1 keV, where it tends to 1.75 % against the
#  1.68 % asymptote (C2 gives 1.60 %).  SO THE COMPARISON THAT COUNTS IS WITH THE MEASUREMENT;  K&M supply
#  sigma(+) and the collected data, and are a target only above 1 keV.
#
#  WHAT THIS SCRIPT ADDS is sigma(+) from THE SAME orbitals, grid, normalization and prefactor, because
#  the ratio is the quantity that cancels a common normalization error while the absolute number does
#  not.  Quoting a ratio built from JAC's sigma(2+) and somebody else's sigma(+) would hide exactly the
#  error the ratio exists to expose.
#
#  sigma(2+) is NOT recomputed here:  it is read in below from the 03-Nov-2026 threshold scan at the same
#  40 a.u. box and the same 9 sharings, which cost 2435 s for this one energy.  Recomputing it to print a
#  table would be four hours for numbers already in hand.
#
#  Run as:  julia --project=. tools/comparePdiHeliumKornberg.jl  [rbox]
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Coulomb functions")

const RBOX   = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 40.0
const grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = RBOX)
const nm     = Nuclear.Model(2.0)
const OMEGA  = 200.0

##  Measured 03-Nov-2026, box 40 a.u., 9 sharings, maxKappa = 2, NoIntermediateEnergies = 9, AFTER the
##  identical-particle 1/2 was supplied.  Coulomb (velocity) and Babushkin (length), in barn.
const SIG2P_COU = 2301.4539
const SIG2P_BAB = 4234.7389

##  Kornberg & Miraglia, PRA 48, 3714 (1993), Table II, in Mb
const KM_SIG1P  = 5.44e-02
const KM_SIG2P  = 1.69e-03        ##  C2:  independent-electron, two Coulomb waves, SH initial state
const KM_SIG2P3 = 1.30e-02        ##  C3:  the correlated final state the paper is about

println("="^128)
@printf("sigma(+) for He at omega = %.0f eV, box %.1f a.u., %d grid points -- the denominator of the K&M ratio\n",
        OMEGA, RBOX, grid.NoPoints)
println("="^128)

##  ----------------------------------------------------------------------------------------------------------
##  The SAME shared orbital set as the double-ionization runs:  one orbital set spanning all three multiplets.
##  ----------------------------------------------------------------------------------------------------------
spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")]
wa          = Representation("He span", nm, grid, spanConfigs, MeanFieldMultiplet(MeanFieldSettings()) )
wb          = redirect_stdout(devnull) do
    generate(wa, output=true)
end
orbitals    = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
println("\n>> shared orbital set: ", join(sort(string.(collect(keys(orbitals)))), ", "), "\n")

pSettings = PhotoIonization.Settings(PhotoIonization.Settings();
                 multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                 printBefore=false, lineSelection=LineSelection(true, indexPairs=[(1,1)]) )
comp = Atomic.Computation(Atomic.Computation(), name="He sigma+", grid=grid, nuclearModel=nm,
                 initialConfigs  = [Configuration("1s^2")],
                 finalConfigs    = [Configuration("1s")],
                 orbitalSource   = Basics.GivenOrbitals(orbitals),
                 processSettings = pSettings )
res   = perform(comp; output=true)
lines = res[PhotoIonization.Settings]

toMb(x) = Defaults.convertUnits("cross section: from atomic to Mbarn", x)
s1C = sum(toMb(l.crossSection.Coulomb)   for l in lines)
s1B = sum(toMb(l.crossSection.Babushkin) for l in lines)
tMb(barn) = barn * 1.0e-6

println("\n" * "="^128)
println("THE COMPARISON   ... Coulomb = velocity, Babushkin = length;  K&M name velocity the trustworthy form")
println("="^128)
@printf("%-26s %14s %14s %14s %14s\n", "", "JAC Coulomb", "JAC Babushkin", "K&M C2", "K&M C3")
println("-"^128)
@printf("%-26s %14.4e %14.4e %14.4e %14.4e\n", "sigma(+)      [Mb]", s1C, s1B, KM_SIG1P, KM_SIG1P)
@printf("%-26s %14.4e %14.4e %14.4e %14.4e\n", "sigma(2+)     [Mb]", tMb(SIG2P_COU), tMb(SIG2P_BAB), KM_SIG2P, KM_SIG2P3)
@printf("%-26s %13.3f %% %13.3f %% %13.3f %% %13.3f %%\n", "sigma(2+)/sigma(+)",
        100*tMb(SIG2P_COU)/s1C, 100*tMb(SIG2P_BAB)/s1B, 100*KM_SIG2P/KM_SIG1P, 100*KM_SIG2P3/KM_SIG1P)
println("-"^128)
println("MEASURED near 200 eV:  sigma(2+)/sigma(+) = 3-4 %, falling to the 1.68 % asymptote above ~1 keV.")
println("Landing between C2 and the data is the right place;  agreeing with C3 (23.9 %) would be a FAILURE.")
println("-"^128)
@printf("%-26s %13.3f  %13.3f\n", "sigma(+)  / K&M", s1C/KM_SIG1P, s1B/KM_SIG1P)
@printf("%-26s %13.3f  %13.3f\n", "sigma(2+) / K&M", tMb(SIG2P_COU)/KM_SIG2P, tMb(SIG2P_BAB)/KM_SIG2P)
@printf("%-26s %13.3f  %13.3f\n", "ratio     / K&M",
        (tMb(SIG2P_COU)/s1C)/(KM_SIG2P/KM_SIG1P), (tMb(SIG2P_BAB)/s1B)/(KM_SIG2P/KM_SIG1P))
println("="^128)
println(">> sigma(2+) is the 03-Nov-2026 threshold-scan value at this box and 9 sharings, not recomputed here.")
println(">> maxKappa = 2 and NoIntermediateEnergies = 9 in that run:  NEITHER is converged, so read the")
println(">> comparison as an order-and-sign statement, not as an accuracy figure.")
