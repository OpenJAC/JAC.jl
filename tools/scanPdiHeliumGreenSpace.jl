#
#  IS THE GREEN SPACE WHAT MAKES sigma(2+) OVERSHOOT, AND WHAT OPENS THE GAUGE GAP?
#
#  By 04-Nov-2026 every numerical truncation in this calculation has been measured and every one is small:
#  the radial box moves sigma(2+) by 0.05 % over 24 -> 40 a.u., the intermediate mesh by 0.4 % over 9 -> 25
#  points, the partial-wave sum by 0.8 % over maxKappa 3 -> 4.  The disagreement with the measured ratio is
#  20-50 %.  So the remaining error is NOT numerical, and the one truncation never varied is the (N-1)-electron
#  GREEN SPACE, which has held four levels -- He+ 1s_1/2, 2s_1/2, 2p_1/2, 2p_3/2 -- since the module was built.
#
#  THE GAUGE RATIO IS THE THERMOMETER, and it already says something sharp.  Across the intermediate-energy
#  mesh it is DEAD FLAT -- 1.840, 1.836, 1.842 at 9, 17, 25 points -- while across the partial-wave sum it
#  walks 1.840 -> 2.346 -> 2.395.  Gauge agreement in a second-order amplitude is a statement about the
#  COMPLETENESS of the intermediate sum, so a gap that is blind to how finely that sum is SAMPLED but
#  sensitive to which final channels are open points at the BOUND states available to it, not at the mesh.
#
#  WHY THIS IS TWO RUNS AND NOT ONE.  There are no 3s, 3p or 3d orbitals to put into the Green space:  the
#  span that generates the shared orbital set stops at n = 2.  Widening it also shifts the 1s, 2s and 2p
#  orbitals themselves, so a single run would move the Green space AND the orbital basis together -- the
#  same entanglement that made the 15-Sep sharing scan unreadable for six weeks and was only separated on
#  04-Nov.  So the widened orbital set is built ONCE and used by BOTH runs:
#
#      RUN B   widened orbitals, Green space 1s,2s,2p          (4 levels)     the ANCHOR
#      RUN C   widened orbitals, Green space 1s,2s,2p,3s,3p,3d (9 levels)     the QUESTION
#
#  B -> C is then the Green space alone.  And B against the maxKappa = 3 run on the OLD orbital set
#  (2797.6991 barn Coulomb, 6563.4778 Babushkin, 04-Nov-2026) is a FREE control on what the widening of the
#  orbital basis did by itself -- a number we need either way, and it costs nothing because we already have it.
#
#  THE PREDICTION IS FALSIFIABLE AND DIRECTIONAL, which is the only reason the run is worth five hours:
#  enlarging the Green space should drive sigma(2+) DOWN towards 2.3e-03 Mb, drive Babushkin/Coulomb DOWN from
#  2.395 towards 1, and concentrate the drop in the MIDDLE of the U, where the maxKappa = 3 channels live.
#  If the total RISES, or the gauge ratio stays at 2.4, the Green space is exonerated and the next suspect is
#  the second-order amplitude itself.
#
#  THE BOX IS 40 a.u., AND A DETOUR TO 60 IS RECORDED HERE BECAUSE THE DETECTOR THAT CAUSED IT WAS WRONG.
#
#  The orbital diagnostic printed below reports, for each generated orbital, the radius at which |P| has
#  fallen to 1e-6 of its maximum. On the first launch at 40 a.u. it flagged SEVEN of nine orbitals as
#  reaching the wall, the five n = 3 members at exactly 1.00 of the box, and that was read as truncation.
#  It is not. Re-running at 60 a.u. moved the energies by almost nothing:
#
#        orbital   box 40       box 60       shift        relative
#        3s        -0.025033    -0.025103    -0.000070    -0.28 %
#        3p        -0.017342    -0.017415    -0.000073    -0.42 %
#        3d        -0.010734    -0.010821    -0.000087    -0.81 %
#        2p        -0.084737    -0.084869    -0.000132    -0.16 %
#        1s        -1.210324    -1.210469    -0.000145    -0.01 %
#
#  The 1e-6 criterion simply cannot tell a truncated orbital from a DIFFUSE one: at a binding of 0.025 a.u.
#  the tail decays as exp(-sqrt(2E) r) with sqrt(2E) = 0.22 /a.u., so reaching 1e-6 of the maximum takes
#  ln(1e6)/0.22 = 62 a.u. beyond the turning point, and a perfectly converged diffuse orbital trips the flag
#  at ANY box. So the flag is uninformative -- but the table above does NOT settle the question either, and
#  the first version of this header claimed it did. The claim was: a hard wall moved from 40 to 60 a.u. sheds
#  pi^2/2 (1/40^2 - 1/60^2) = 0.00171 a.u., which is 7 % of the 3s energy, whereas the 3s shed 0.00007. That
#  arithmetic is a FREE-PARTICLE yardstick, and it overstates by a large factor the energy a COULOMB-BOUND
#  orbital gives up when the wall is pushed out beyond its turning point -- so "the energies barely moved"
#  is consistent with truncation as well, and proves nothing. BOTH detectors were miscalibrated, in opposite
#  directions, and the second was built to refute the first.
#
#  WHAT DOES DECIDE IT was already settled in the project's own register on 17-Oct-2026 and amended
#  21-Oct-2026, and was found only after this run was launched twice: `extent/box` OVER-WARNS because it
#  measures where a tail EXISTS and not its WEIGHT, and the calibrated discriminator is Z_eff TOGETHER WITH
#  n >= 3. Both hold here, so the n = 3 orbitals ARE genuinely under-boxed at 40 a.u. and at 60 a.u. alike.
#  That is a caveat on the SIZE of any B -> C movement and not on its SIGN, for the reason in the next
#  paragraph; a wide-box re-run of RUN C is the follow-up if B -> C turns out large.
#
#  WHAT THE DIAGNOSTIC DOES SHOW is worth more than what it was built for. Read as effective charges, the
#  energies give Z_eff = 0.672 (3s), 0.560 (3p) and 0.441 (3d), where an outer electron of neutral He should
#  see ~1 and a He+ Green-space state should see 2. The configuration-averaged mean field OVER-SCREENS these
#  orbitals, so the n = 3 functions entering the Green space are far more diffuse than the He+ states they
#  stand for. That is the price of the one-shared-orbital-set precondition, which the 2s and 2p already pay,
#  and it is a caveat on INTERPRETING B -> C rather than a defect to repair here: part of any movement may be
#  these over-diffuse functions rather than genuine n = 3 intermediate structure.
#
#  AND 40 a.u. IS THE RIGHT BOX FOR THIS COMPARISON FOR A REASON INDEPENDENT OF ALL THAT: 2s and 2p sit in
#  BOTH Green spaces, so whatever the box does to them is common to B and C and cancels in the difference,
#  while 40 a.u. keeps RUN B a SINGLE-VARIABLE control against the 2797.6991 barn point -- the orbital span
#  alone. At 60 a.u. that control carries the box as well, and it costs twice the machine time: measured on
#  the aborted launch, 1360 of ~11754 continuum orbitals in 25 min, i.e. ~3.6 h for RUN B against ~1.5 h.
#
#  Run as:  julia --project=. tools/scanPdiHeliumGreenSpace.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Coulomb functions")

const RBOX   = 40.0              ##  see the Rule 12 paragraph in the header;  60 was tried and reverted
const grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = RBOX)
const nm     = Nuclear.Model(2.0)
const OMEGA  = 200.0
const NSHARE = 9
const NINT   = 9
const MAXKAP = 3                  ##  converged to 0.8 %;  maxKappa 4 costs 29 % more for that 0.8 %

##  The maxKappa = 3 run on the NARROW orbital set, 04-Nov-2026 -- the free control for RUN B
const OLD_C  = 2797.6991
const OLD_B  = 6563.4778

println("="^128)
println("PhotoDoubleIonization: does the GREEN SPACE explain the overshoot and the gauge gap?")
@printf("He at %.0f eV;  box %.1f a.u., %d grid points;  %d sharings;  NoIntermediateEnergies = %d;  maxKappa = %d\n",
        OMEGA, RBOX, grid.NoPoints, NSHARE, NINT, MAXKAP)
println("RUN B = widened orbitals + 4-level Green space (anchor);  RUN C = widened orbitals + 9-level Green space")
@printf("Control, narrow orbitals + 4-level Green space at this maxKappa:  %.4f barn (Cou), %.4f (Bab), ratio %.3f\n",
        OLD_C, OLD_B, OLD_B/OLD_C)
println("="^128)

##  ----------------------------------------------------------------------------------------------------------
##  THE WIDENED SHARED ORBITAL SET, built once and used by both runs and by sigma(+).
##  ----------------------------------------------------------------------------------------------------------
spanConfigs = [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p"),
               Configuration("1s 3s"), Configuration("1s 3p"), Configuration("1s 3d")]
wa          = Representation("He span", nm, grid, spanConfigs, MeanFieldMultiplet(MeanFieldSettings()) )
wb          = redirect_stdout(devnull) do
    generate(wa, output=true)
end
orbitals    = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
println("\n>> widened orbital set: ", join(sort(string.(collect(keys(orbitals)))), ", "))

##  THE RULE 12 DIAGNOSTIC:  where does each orbital actually end, against the wall at RBOX?  An orbital whose
##  last significant amplitude sits at the box edge is squeezed, and the whole point is to SEE that rather than
##  to discover it later in a number.  Wrapped, because a diagnostic must never kill a five-hour job.
try
    println("\n  ", rpad("subshell",12), rpad("energy [a.u.]",18), rpad("r at |P| < 1e-6 max",22), "fraction of box")
    for sh in sort(collect(keys(orbitals)), by=string)
        orb = orbitals[sh];   np = min(length(orb.P), length(grid.r))
        pmx = maximum(abs, view(orb.P, 1:np))
        idx = findlast(i -> abs(orb.P[i]) > 1.0e-6*pmx, 1:np)
        rr  = idx === nothing ? NaN : grid.r[idx]
        @printf("  %-12s%-18.6f%-22.3f%.2f %s\n", string(sh), orb.energy, rr, rr/RBOX,
                rr/RBOX > 0.9 ? "  <-- reaches the wall;  diffuse OR truncated, see header" : "")
    end
catch e
    println("  >> orbital-extent diagnostic unavailable: ", e)
end

##  ----------------------------------------------------------------------------------------------------------
##  sigma(+) on the WIDENED orbital set -- the ratio's denominator must come from the same orbitals as its
##  numerator, or the ratio is not a ratio.  Cheap, so it goes first and the job has a number within minutes.
##  ----------------------------------------------------------------------------------------------------------
pSettings = PhotoIonization.Settings(PhotoIonization.Settings();
                 multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                 printBefore=false, lineSelection=LineSelection(true, indexPairs=[(1,1)]) )
pComp = Atomic.Computation(Atomic.Computation(), name="He sigma+ widened", grid=grid, nuclearModel=nm,
                 initialConfigs  = [Configuration("1s^2")],
                 finalConfigs    = [Configuration("1s")],
                 orbitalSource   = Basics.GivenOrbitals(orbitals),
                 processSettings = pSettings )
pLines = perform(pComp; output=true)[PhotoIonization.Settings]
toMb(x) = Defaults.convertUnits("cross section: from atomic to Mbarn", x)
const SIG1P_C = sum(toMb(l.crossSection.Coulomb)   for l in pLines)
const SIG1P_B = sum(toMb(l.crossSection.Babushkin) for l in pLines)
@printf("\n>>> RESULT  sigma(+) on the widened orbitals:  %.4e Mb (Cou),  %.4e (Bab)    [narrow set gave 5.9499e-02]\n\n",
        SIG1P_C, SIG1P_B)

##  ----------------------------------------------------------------------------------------------------------
##  The two double-ionization runs.  Each prints its own result line the moment it lands, so a job that dies
##  halfway still leaves evidence -- the lesson of the /tmp reboot of 03-Sep-2026.
##  ----------------------------------------------------------------------------------------------------------
tob(x) = Defaults.convertUnits("cross section: from atomic to barn", x)
RUNS = [ ("B anchor", [Configuration("1s"), Configuration("2s"), Configuration("2p")]),
         ("C green",  [Configuration("1s"), Configuration("2s"), Configuration("2p"),
                       Configuration("3s"), Configuration("3p"), Configuration("3d")]) ]
results = Tuple{String,Int64,Float64,Float64,Any}[]

for (label, greenConfigs) in RUNS
    println("-"^128)
    @printf("RUN %s:  Green space from %s\n", label, join(string.(greenConfigs), ", "))
    println("-"^128)
    t = time()
    gMultiplet = redirect_stdout(devnull) do
        Hamiltonian.performCIwithFrozenOrbitals(greenConfigs, orbitals, nm, grid, AsfSettings(); printout=false)
    end
    @printf(">> Green space: %d levels, %d electron(s)\n", length(gMultiplet.levels),
            gMultiplet.levels[1].basis.NoElectrons)
    dSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                     multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                     NoEnergySharings=NSHARE, NoIntermediateEnergies=NINT, maxKappa=MAXKAP,
                     calcDifferentialCs=true, printBefore=false,
                     lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                     eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )
    comp = Atomic.Computation(Atomic.Computation(), name="He sigma2+ $label", grid=grid, nuclearModel=nm,
                     initialConfigs  = [Configuration("1s^2")],
                     finalConfigs    = [Configuration("1s^0")],
                     orbitalSource   = Basics.GivenOrbitals(orbitals),
                     processSettings = dSettings )
    line = perform(comp; output=true)[PhotoDoubleIonization.Settings][1]
    cC   = tob(line.crossSection.Coulomb);    cB = tob(line.crossSection.Babushkin)
    push!(results, (label, length(gMultiplet.levels), cC, cB, line.sharings))
    ##  ONE literal format string:  @printf rejects a concatenation, and splitting this line to stay inside
    ##  140 columns is what killed the 16:47 launch after it had already computed sigma(+).
    @printf("\n>>> RESULT %-8s %d-level Green  sigma2+ %10.4f (Cou) %10.4f (Bab) barn  Bab/Cou %5.3f  ratio %6.3f %%  [%.0f s]\n\n",
            label, length(gMultiplet.levels), cC, cB, cB/cC, 100*cC*1.0e-6/SIG1P_C, time()-t)
end

##  ----------------------------------------------------------------------------------------------------------
##  The two tables.  Totals carry the answer;  the per-node column says WHERE the Green space acts, which is
##  what distinguishes a missing intermediate channel from a scale factor.
##  ----------------------------------------------------------------------------------------------------------
println("="^128)
println("TABLE 1 -- THE ANSWER.  Prediction: C below B, Bab/Cou towards 1.  Either falsifies or confirms it.")
println("="^128)
@printf("%-12s %8s %16s %16s %10s %12s %12s\n", "run", "levels", "sigma2+ Cou", "sigma2+ Bab", "Bab/Cou", "ratio", "vs B")
@printf("%-12s %8d %16.4f %16.4f %10.3f %11s %12s\n", "narrow ctl", 4, OLD_C, OLD_B, OLD_B/OLD_C, "--", "--")
for (label, nl, cC, cB, _) in results
    vsB = label == "B anchor" ? "--" : @sprintf("%+8.2f %%", 100*(cC-results[1][3])/results[1][3])
    @printf("%-12s %8d %16.4f %16.4f %10.3f %10.3f %% %12s\n", label, nl, cC, cB, cB/cC,
            100*cC*1.0e-6/SIG1P_C, vsB)
end
@printf("\nOrbital-basis control, narrow -> widened at a 4-level Green space:  Coulomb %+.2f %%,  Babushkin %+.2f %%\n",
        100*(results[1][3]-OLD_C)/OLD_C, 100*(results[1][4]-OLD_B)/OLD_B)
@printf("Measured near 200 eV: 3-4 %%.  sigma(+) used as the denominator: %.4e Mb (Coulomb, widened orbitals).\n", SIG1P_C)

println("\n" * "="^128)
println("TABLE 2 -- WHERE IT ACTS.  The partial-wave sum fills the MIDDLE of the U;  if the Green space is the")
println("missing support for those channels, it must take the middle back out again.")
println("="^128)
##  `differentialCs` is held in ATOMIC UNITS (a_0^2 per Hartree).  The run of 04-Nov-2026 printed the raw
##  field under a "[barn/eV]" header and at %.4f, so the middle of the U came out as 0.0000 and the table
##  was unreadable -- the same unit-label defect that was fixed in this module's own displayResults that
##  morning, reproduced here within hours.  Convert with the module's OWN conversion rather than a
##  hand-rolled factor, so that the two cannot drift apart.  The CHANGE column and the U-ratio were and
##  remain unaffected, being ratios of like quantities.
@printf("  %12s %16s %16s %12s\n", "eps1 [eV]", "B [barn/eV]", "C [barn/eV]", "change")
for i in 1:length(results[1][5])
    e1 = Defaults.convertUnits("energy: from atomic to eV", results[1][5][i].epsilon1)
    dB = Defaults.convertUnits("energy-diff. cross section: from atomic", results[1][5][i].differentialCs.Coulomb)
    dC = Defaults.convertUnits("energy-diff. cross section: from atomic", results[2][5][i].differentialCs.Coulomb)
    @printf("  %12.3f %16.6e %16.6e %11.2f %%\n", e1, dB, dC, 100*(dC-dB)/dB)
end
@printf("\n  U-ratio, outermost/equal sharing:  B = %.3f,  C = %.3f   (narrow orbitals at this maxKappa: 4.481)\n",
        results[1][5][1].differentialCs.Coulomb / results[1][5][cld(NSHARE,2)].differentialCs.Coulomb,
        results[2][5][1].differentialCs.Coulomb / results[2][5][cld(NSHARE,2)].differentialCs.Coulomb)
println("\n" * "="^128)
