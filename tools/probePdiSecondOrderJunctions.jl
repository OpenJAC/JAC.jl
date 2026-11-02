#
#  [PDI-J]  Which operator belongs at the junctions of PhotoDoubleIonization's second-order sum, and why the
#           answer is "leave it alone".  The measurement that settled it on 02-Nov-2026.
#
#  WHY THIS EXISTS.  A one-orbital junction of a second-order sum needs H - E and not V^(e-e) alone:  the
#  one-body part contributes <b|h|a>, which cancels the spectator mean field that V_ee counts a second time.
#  That argument is right, it is written out in ManyElectron.matrixElement_HminusE, and it is what sent
#  TwoElectronOnePhoton to H - E.  PhotoDoubleIonization's junction 1 reaches exactly such pairs, so it looked
#  like the same repair.  IT IS NOT, and the reason is a term that neither operator shows at the call site.
#
#  MEASURED 02-Nov-2026, He-like Be as examples/example-Dj.jl computes it -- Green space [1s, 2s, 2p] of Be^3+
#  and initial state 1s^2, each from its OWN SCF, both on the one grid that branch uses since `947880e4`:
#
#      junction 1, Green core 1s   nDiff 1    V_ee +2.470622e-01   H - E +1.560546e-01   ratio 1.583
#      junction 1, Green core 2s   nDiff 2    V_ee -6.731512e-02   H - E -6.731512e-02   ratio 1.000
#      one-body part at the nDiff = 1 pair                                   -9.100758e-02   (37 % of V_ee)
#
#  SO THE ONE-BODY TERM BEHAVES:  right sign, sensible size, and the nDiff = 2 pairs are untouched to 1.000
#  exactly -- which is also the control proving a swap could reach only the one-orbital junctions.  The ratio
#  1.583 SUPERSEDES a 1.690 measured before the grid repair, which was therefore never a result.
#
#  WHAT KILLS THE SWAP IS THE DROPPED OVERLAP TERM.  matrixElement_HminusE omits -E <r|s> on a stated
#  orthonormality assumption.  Here the two orbital sets are different, so it does not vanish:
#
#      <1s_green|1s_initial>  = +0.994787        <2s_green|1s_initial> = +0.081996   (worst departure 8.2e-02)
#      <101s_1/2(eps = 0.5076)|1s_initial> = -0.023524,  E_i = -13.684606 Ha
#      |E_i <r|s>| ~ 3.2e-01  =  2.05 x the whole H - E element  =  3.5 x the correction a swap would add
#
#  The continuum orbital is generated in the field of the GREEN level, hence orthogonal to the GREEN bound set
#  and NOT to the initial state's -- and that single overlap of -0.0235, times the initial energy, is twice the
#  matrix element it would correct.  A swap would ADD 0.091 and SILENTLY OMIT 0.32.  Unlike the
#  TwoElectronOnePhoton case this is not a sign flip;  it fails the other way, and a known factor 1.583 would be
#  traded for an unquantified larger one.  PhotoDoubleIonization therefore KEEPS V_ee, checks the precondition in
#  `checkSecondOrderConditions`, and waits for a Green space built on the initial state's own orbitals.  Building
#  one needs Hamiltonian.performCIwithFrozenOrbitals, which Atomic.Computation + perform cannot express today.
#
#  SECTION 4 IS THE VACUITY CHECK, and the negatives above are worth nothing without it:  a guard that cannot
#  fail is decoration, so the same case is re-run with the Green multiplet rebuilt on the initial state's
#  orbitals, where the check MUST fall silent.
#
#  RUN:  julia --project=. tools/probePdiSecondOrderJunctions.jl
#
using JenaAtomicCalculator, Printf, LinearAlgebra
const JAC = JenaAtomicCalculator

grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
nm     = Nuclear.Model(4.01);           asf    = AsfSettings()
nucPot = Nuclear.nuclearPotential(nm, grid)
prims  = Bsplines.generatePrimitives(grid)

##  the Green space exactly as example-Dj.jl branch 1 builds it, and the initial state as the Atomic.Computation does
wa  = Representation("He-like Be", nm, grid, [Configuration("1s"), Configuration("2s"), Configuration("2p")],
                     MeanFieldMultiplet(MeanFieldSettings()) )
gMp = (redirect_stdout(devnull) do;   generate(wa, output=true);   end)[ResultKeys.MeanFieldMultiplet]
iMp = redirect_stdout(devnull) do
    SelfConsistent.performSCF([Configuration("1s^2")], nm, grid, asf; printout=false)
end
iLevel0 = iMp.levels[1];        symi = LevelSymmetry(iLevel0.J, iLevel0.parity)

##  SECTION 1: do the Green space and the initial state share one orbital set?
println("\n######## condition (b): one orbital set for the Green space and the initial state ########")
@printf("   Green basis   %s\n", join(sort(string.(collect(keys(gMp.levels[1].basis.orbitals)))), ", "))
@printf("   initial basis %s        E_i = %+.9f Ha\n\n", join(sort(string.(collect(keys(iLevel0.basis.orbitals)))), ", "),
        iLevel0.energy)
worst = 0.;     worstName = ""
for  sha  in  sort(collect(keys(gMp.levels[1].basis.orbitals)), by=string),
     shb  in  sort(collect(keys(iLevel0.basis.orbitals)), by=string)
    sha.kappa == shb.kappa  ||  continue
    ovl = JAC.RadialIntegrals.overlap(gMp.levels[1].basis.orbitals[sha], iLevel0.basis.orbitals[shb], grid)
    dev = sha == shb ? abs(abs(ovl) - 1.) : abs(ovl)
    @printf("   <%-8s _green | %-8s _initial> = %+.9f     departure from one orbital set %.3e\n",
            string(sha), string(shb), ovl, dev)
    dev > worst  &&  (global worst = dev;    global worstName = "<$(string(sha))|$(string(shb))>")
end
@printf("\n   WORST departure %.4e at %s   -- the condition needs this below 1e-4\n", worst, worstName)

##  SECTION 2: the two operators at junction 1, one sharing at omega = 400 eV
println("\n######## the two operators at junction 1 ########")
omega  = Defaults.convertUnits("energy: to atomic", 400.)
e1     = e2 = (omega + iLevel0.energy) / 2.
@printf("   omega = %.6f Ha,   eps1 = eps2 = %.6f Ha\n\n", omega, e1)
contSettings = Continuum.Settings(false, grid.NoPoints-100)
sh1 = Subshell(101, -1);        sh2 = Subshell(102, -1)
subshells = copy(iLevel0.basis.subshells);      push!(subshells, sh1);      push!(subshells, sh2)
newiLevel = Basics.generateLevelWithSymmetryReducedBasis(iLevel0, subshells)
nLevels, nWeights, nOnShell, nBlock, nEnergy = redirect_stdout(devnull) do
    PhotoDoubleIonization.generateIntermediateLevels(symi, gMp, nm, grid,
            Radial.GridGL(Radial.GridGaussLegendreFinite(), 0.01, 1., 2; printout=false), contSettings, 2;
            nuclearPot=nucPot, primitives=prims, spectators=[(sh1, e1), (sh2, e2)])
end
##  the same subshell merge that PhotoDoubleIonization.amplitude performs before evaluating the sum
subs = copy(newiLevel.basis.subshells)
for  nL  in  nLevels    global subs = Basics.merge(subs, nL.basis.subshells)    end
iLevel = Level(newiLevel, subs)
@printf("   %-26s %-5s %14s %14s %14s %11s\n", "intermediate CSF", "nDiff", "V_ee", "H - E", "one-body", "V_ee/(H-E)")
for  nL0  in  nLevels
    nL = Level(nL0, subs)
    for  tp = 1:length(nL.basis.csfs),  s = 1:length(iLevel.basis.csfs)
        occN = nL.basis.csfs[tp].occupation
        nD   = sum(abs.(occN .- iLevel.basis.csfs[s].occupation)) / 2
        vee  = ManyElectron.matrixElement_Vee(CoulombInteraction(), nL.basis, tp, iLevel.basis, s, grid)
        hme  = ManyElectron.matrixElement_HminusE(Basics.HminusE(CoulombInteraction()), nL.basis, tp,
                                                  iLevel.basis, s, grid, nucPot)
        lab  = join([ "$(string(sh))^$(o)" for (sh, o) in zip(subs, occN) if o != 0 ], " ")
        @printf("   %-26s %-5.0f %+14.6e %+14.6e %+14.6e %11.3f\n", lab, nD, vee, hme, hme - vee,
                abs(hme) < 1.0e-30 ? NaN : vee/hme)
    end
end
println("\n   `one-body` is H - E minus V_ee.  The nDiff = 2 rows MUST read 1.000 exactly:  a one-body operator")
println("   cannot connect CSFs differing in two orbitals, so they bound what any operator swap could reach.")

##  SECTION 3: the term matrixElement_HminusE drops, and its size against the element it would correct
println("\n######## the dropped overlap term -E <r|s> ########")
gL = Basics.generateLevelWithSymmetryReducedBasis(gMp.levels[1], gMp.levels[1].basis.subshells)
cOrbital, cPhase = redirect_stdout(devnull) do
    Continuum.generateOrbitalForLevel(e1, sh1, gL, nm, grid, contSettings; nuclearPot=nucPot, primitives=prims)
end
ovBound = JAC.RadialIntegrals.overlap(gMp.levels[1].basis.orbitals[Subshell(1,-1)], iLevel0.basis.orbitals[Subshell(1,-1)], grid)
ovCont  = JAC.RadialIntegrals.overlap(cOrbital, iLevel0.basis.orbitals[Subshell(1,-1)], grid)
estimate = abs( iLevel0.energy * ovBound * ovCont )
@printf("   <1s_green|1s_initial>                   = %+.9f\n", ovBound)
@printf("   <%s(eps = %.4f)|1s_initial>   = %+.9f   -- zero only within ONE orbital set\n",
        string(sh1), e1, ovCont)
@printf("\n   |E_i <r|s>| ~ |%.6f| * %.6f * %.6f = %.6e\n", iLevel0.energy, ovBound, ovCont, estimate)
@printf("   against the junction's own H - E = 1.560546e-01 and one-body part -9.100758e-02:\n")
@printf("      dropped / (H - E)  = %.3f        dropped / one-body = %.3f\n", estimate/1.560546e-01,
        estimate/9.100758e-02)
println("\n   This is an ESTIMATE: the true CSF overlap carries an angular factor and a determinant of orbital")
println("   overlaps, and only the leading product is formed here.  It SIZES the dropped term, it does not")
println("   evaluate it -- and a term twice the element it would correct is already decisive.")

##  SECTION 4: the vacuity check, without which none of the above means anything
println("\n######## calibration by mutation of checkSecondOrderConditions ########")
gFixed = Multiplet(gMp.name,
                   [Level(l.J, l.M, l.parity, l.index, l.energy, l.relativeOcc, l.hasStateRep,
                          Basis(l.basis.isDefined, l.basis.NoElectrons, l.basis.subshells, l.basis.csfs,
                                l.basis.coreSubshells, iLevel0.basis.orbitals), l.mc)
                    for l in gMp.levels])
for  (tag, gm)  in  ( ("separate SCF sets, as example-Dj.jl", gMp), ("Green space on the initial orbitals", gFixed) )
    problems = PhotoDoubleIonization.checkSecondOrderConditions(iMp, gm, grid)
    @printf("\n   %-38s %d problem(s) reported\n", tag, length(problems))
    for  problem in problems    println("      ", problem)    end
end
println("\n   The first case MUST report and the second MUST NOT.  If both fall silent the check has stopped")
println("   testing anything, and the measurement above would no longer be reproducible from this file.")
