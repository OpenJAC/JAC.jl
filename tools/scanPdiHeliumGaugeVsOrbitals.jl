#
#  IS THE SECOND-ORDER GAUGE GAP A FUNCTION OF THE FIRST-ORDER ONE?
#
#  By 05-Nov-2026 the PhotoDoubleIonization overshoot has survived every numerical truncation (all under
#  1 %) and the enlargement of the (N-1) Green space (-0.01 % from 4 levels to 9).  What has NOT been
#  varied is the shared ORBITAL SET, and a first-order scan that morning showed it moves the FIRST-ORDER
#  gauge agreement of sigma(+) over a wide range -- 4.19 % to 35.67 % across thirteen sets -- while the
#  SECOND-ORDER gap on the production set is 139 %.
#
#  SO THE QUESTION IS NOT "WHICH SET IS BEST" BUT WHETHER THE TWO GAPS ARE RELATED AT ALL:
#
#      if the second-order gap TRACKS the first-order one, the second-order construction merely AMPLIFIES
#          an orbital inconsistency, and item 6 is about the orbital design;
#      if it sits near 140 % whatever the orbitals do, the construction is defective INDEPENDENTLY of
#          them, and item 6 is about the amplitude.
#
#  Five points, and two of them are FREE because they were already measured at this maxKappa:
#
#      orbital set            eps(1s)     1st-order gauge dev    sigma(+) vs K&M    2nd-order
#      DFS, span to n=2       -1.1008           13.33 %               +9.4 %        2.346  (04-Nov)
#      DFS, span to n=3       -1.2103           13.13 %              +10.4 %        2.377  (04-Nov)
#      HS,  span to n=2       -1.0033           17.47 %               -0.0 %        this run
#      AL,  span to n=3       -1.7144            6.99 %              +21.9 %        this run
#      V^(N-1), He+ span      -2.0001            4.19 %              +29.6 %        this run
#
#  THE GREEN SPACE IS HELD AT FOUR LEVELS FOR ALL OF THEM, AND THAT IS WHAT MAKES THIS AFFORDABLE.  The
#  run of 04-Nov-2026 opened it from 4 levels to 9 and moved sigma(2+) by -0.01 %, so paying for 9 levels
#  here would triple the machine time to measure a quantity already shown to be insensitive.  One 3.7 h
#  run becomes three of ~1.4 h, which is the whole reason a null result is worth having.
#
#  WHAT IS READ, AND WHAT MUST NOT BE:  read Bab/Cou.  DO NOT QUOTE THE RATIO sigma(2+)/sigma(+) FROM THE
#  V^(N-1) OR AL ROWS -- their sigma(+) is 22-30 % away from the literature BY CONSTRUCTION, since a set
#  chosen for gauge consistency is not the set that reproduces the one-electron observable.  The ratio is
#  printed for the HS row, whose sigma(+) is anchored to 0.0 %, and is printed as a dash for the others.
#
#  Run as:  julia --project=. tools/scanPdiHeliumGaugeVsOrbitals.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Coulomb functions")

const RBOX   = 40.0
const grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = RBOX)
const nm     = Nuclear.Model(2.0)
const OMEGA  = 200.0
const NSHARE = 9
const NINT   = 9
const MAXKAP = 3
const KM1P   = 5.44e-02           ##  Mb, Kornberg & Miraglia's own sigma(+) at 200 eV

##  The He ground state is ALWAYS 1s^2 and the Green space is ALWAYS these three configurations;  only the
##  ORBITALS they are expanded on change from row to row.
const GREEN  = [Configuration("1s"), Configuration("2s"), Configuration("2p")]

##  (label, span configurations, mean field, first-order gauge deviation measured 05-Nov-2026, quote ratio?)
const SETS = [
    ("HS n=2",   [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p")],
                 Basics.HSField(), 17.47, true),
    ("AL n=3",   [Configuration("1s^2"), Configuration("1s 2s"), Configuration("1s 2p"),
                  Configuration("1s 3s"), Configuration("1s 3p"), Configuration("1s 3d")],
                 Basics.ALField(),  6.99, false),
    ("V^(N-1)",  [Configuration("1s"), Configuration("2s"), Configuration("2p"),
                  Configuration("3s"), Configuration("3p"), Configuration("3d")],
                 Basics.ALField(),  4.19, false) ]

##  A label given on the command line restricts the run to that set, so that two sets can run as two
##  processes (two at a time is the standing limit -- each Julia process holds a couple of gigabytes).
const WANTED = isempty(ARGS) ? String[] : ARGS
const TODO   = isempty(WANTED) ? SETS : [s for s in SETS if s[1] in WANTED]
isempty(TODO) && error("no set matches $(WANTED);  labels are $(join([s[1] for s in SETS], ", "))")

toMb(x) = Defaults.convertUnits("cross section: from atomic to Mbarn", x)
tob(x)  = Defaults.convertUnits("cross section: from atomic to barn",  x)

println("="^128)
println("PhotoDoubleIonization: is the SECOND-ORDER gauge gap a function of the FIRST-ORDER one?")
@printf("He at %.0f eV;  box %.1f a.u., %d grid points;  %d sharings;  NINT = %d;  maxKappa = %d\n",
        OMEGA, RBOX, grid.NoPoints, NSHARE, NINT, MAXKAP)
println("Green space held at FOUR levels throughout (1s, 2s, 2p) -- licensed by the 4 -> 9 null of 04-Nov-2026")
println("Free points at this maxKappa:  DFS n=2  Bab/Cou 2.346 (1st order 13.33 %);  DFS n=3  2.377 (13.13 %)")
println("="^128)

results = Tuple{String,Float64,Float64,Float64,Float64,Float64,Bool}[]

for (label, spanConfigs, field, dev1, quoteRatio) in TODO
    println("-"^128)
    @printf("SET %s:  %d span configurations, %s\n", label, length(spanConfigs), string(typeof(field)))
    println("-"^128)
    t = time()

    wa  = Representation("He span", nm, grid, spanConfigs, MeanFieldMultiplet(MeanFieldSettings(field)))
    wb  = redirect_stdout(devnull) do
        generate(wa, output=true)
    end
    orbitals = wb[ResultKeys.MeanFieldMultiplet].levels[1].basis.orbitals
    e1s = orbitals[[s for s in keys(orbitals) if string(s) == "1s_1/2"][1]].energy
    @printf(">> orbital set: %s\n>> eps(1s) = %.6f Ha\n", join(sort(string.(collect(keys(orbitals)))), ", "), e1s)

    ##  sigma(+) on THESE orbitals -- the denominator must share the numerator's orbitals or it is not a ratio
    pSettings = PhotoIonization.Settings(PhotoIonization.Settings();
                     multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                     printBefore=false, lineSelection=LineSelection(true, indexPairs=[(1,1)]) )
    pComp = Atomic.Computation(Atomic.Computation(), name="He sigma+", grid=grid, nuclearModel=nm,
                     initialConfigs  = [Configuration("1s^2")],
                     finalConfigs    = [Configuration("1s")],
                     orbitalSource   = Basics.GivenOrbitals(orbitals),
                     processSettings = pSettings )
    pLines = redirect_stdout(devnull) do
        perform(pComp; output=true)[PhotoIonization.Settings]
    end
    s1C = sum(toMb(l.crossSection.Coulomb)   for l in pLines)
    s1B = sum(toMb(l.crossSection.Babushkin) for l in pLines)
    @printf(">> sigma(+) = %.4e Mb (Cou), %.4e (Bab);  1st-order gauge dev %.2f %%;  vs K&M %+.1f %%\n",
            s1C, s1B, 100*abs(s1B-s1C)/s1C, 100*(s1C-KM1P)/KM1P)

    gMultiplet = redirect_stdout(devnull) do
        Hamiltonian.performCIwithFrozenOrbitals(GREEN, orbitals, nm, grid, AsfSettings(); printout=false)
    end
    dSettings = PhotoDoubleIonization.Settings(PhotoDoubleIonization.Settings();
                     multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], photonEnergies=[OMEGA],
                     NoEnergySharings=NSHARE, NoIntermediateEnergies=NINT, maxKappa=MAXKAP,
                     calcDifferentialCs=true, printBefore=false,
                     lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                     eeInteraction=CoulombInteraction(), gMultiplet=gMultiplet )
    dComp = Atomic.Computation(Atomic.Computation(), name="He sigma2+", grid=grid, nuclearModel=nm,
                     initialConfigs  = [Configuration("1s^2")],
                     finalConfigs    = [Configuration("1s^0")],
                     orbitalSource   = Basics.GivenOrbitals(orbitals),
                     processSettings = dSettings )
    dLines = perform(dComp; output=true)[PhotoDoubleIonization.Settings]
    cC = sum(tob(l.crossSection.Coulomb)   for l in dLines)
    cB = sum(tob(l.crossSection.Babushkin) for l in dLines)
    push!(results, (label, dev1, e1s, cC, cB, s1C, quoteRatio))

    ##  ONE literal format string:  a concatenated one raises at macro expansion, which cost a launch on 04-Nov.
    @printf("\n>>> RESULT %-9s 1st-order %5.2f %%   sigma2+ %10.4f (Cou) %10.4f (Bab) barn   Bab/Cou %6.3f   2nd-order %7.2f %%   [%.0f s]\n\n",
            label, dev1, cC, cB, cB/cC, 100*(cB-cC)/cC, time()-t)
end

println("\n" * "="^128)
println("THE ANSWER.  Does the second-order gauge gap follow the first-order one?")
println("="^128)
@printf("  %-9s %12s %12s %14s %10s %12s %12s\n", "set", "eps(1s)", "1st order", "2nd order", "Bab/Cou",
        "sigma2+ Cou", "ratio")
println("  " * "-"^88)
@printf("  %-9s %12.4f %11.2f %% %13.2f %% %10.3f %12.4f %12s\n", "DFS n=2", -1.1008, 13.33, 132.7, 2.346,
        2797.6991, "4.702 %")
@printf("  %-9s %12.4f %11.2f %% %13.2f %% %10.3f %12.4f %12s\n", "DFS n=3", -1.2103, 13.13, 137.7, 2.377,
        2902.2322, "4.831 %")
for (label, dev1, e1s, cC, cB, s1C, quoteRatio) in results
    rs = quoteRatio ? @sprintf("%.3f %%", 100*cC*1.0e-6/s1C) : "--"
    @printf("  %-9s %12.4f %11.2f %% %13.2f %% %10.3f %12.4f %12s\n", label, e1s, dev1,
            100*(cB-cC)/cC, cB/cC, cC, rs)
end
println()
println("  A ratio is printed only where sigma(+) is anchored to the literature;  elsewhere it is a dash,")
println("  because a set chosen for gauge consistency does not reproduce the one-electron observable.")
println("="^128)
