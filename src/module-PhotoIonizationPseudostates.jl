
"""
`module  JAC.PhotoIonizationPseudostates`
... a submodel of JAC that treats resonant photoionization by DIAGONALIZING ONCE a CI matrix that mixes the
    bound autoionizing ("doubly-excited") configurations together with the ejected electron's BOX-DISCRETIZED
    continuum -- a ladder of pseudostate orbitals of one fixed partial wave kappa, all mutually orthogonal
    members of ONE generalized one-electron eigenproblem (`Bsplines.diagonalizeLocalMatrix`) -- rather than
    computing the direct and the resonant amplitude SEPARATELY and adding them through Fano's closed-form
    q-formula, as `PhotoIonizationInterference` does.

    THIS IS THE DIRECT ALTERNATIVE TO `PhotoIonizationInterference`, NOT A REPLACEMENT OF IT, built for the
    SAME physical observable after that module's own q-extraction got stuck on a phase-convention mismatch
    among three separately-computed amplitudes (challenges-list S29). Diagonalizing once removes the need for
    the q-formula, and with it the phase convention it assumes, entirely: the interference between the direct
    and the resonant route lives IN THE EIGENVECTORS of the combined CI, not in a sum of two amplitudes that
    must be put on a common phase footing by hand.

    WHERE THE BOUND, AUTOIONIZING CHARACTER COMES FROM -- deliberately NOT a `Settings` field here. A caller
    hands `computeSpectra` an `initialLevel::Level` whose own `.basis` may ALREADY be a multi-configuration CI
    spanning both the ordinary bound configuration and any doubly-excited ("autoionizing") ones -- exactly how
    `PhotoIonizationInterference`'s own `intermediateMultiplet` is built, just folded into ONE multiplet
    instead of two. Everything already sitting in `initialLevel.basis.csfs` is carried, unchanged, into the
    combined basis; only the PSEUDOSTATE part is new. This is simpler than an earlier sketch that gave this
    module its own `autoionizingConfigs` field: that field would have duplicated what `initialLevel.basis`
    already knows, and merging it at the `Configuration`/`Shell` level would have reintroduced exactly the
    kappa-splitting problem the pseudostate ladder itself had to avoid (see below).

    WHY THE PSEUDOSTATE LADDER CANNOT USE `AtomicState.RasStep`/`Basics.Configuration`. A pseudostate ladder is,
    by construction, ONE fixed kappa -- the ejected electron's own partial wave. `Basics.shellSplitOccupation`
    always splits a non-relativistic `l > 0` shell into BOTH kappa branches (`j = l-1/2` and `j = l+1/2`), so
    routing the ladder through the ordinary `Configuration`/`Shell` machinery would silently double every
    pseudostate into two unwanted kappa families. The ladder is therefore built directly at the `ConfigurationR`
    level -- one relativistic configuration per retained pseudostate orbital, each the final ion's own core
    occupation plus one electron in `Subshell(1000+k, kappa)` -- merged with `initialLevel.basis`'s own CSFs
    through the SAME `Basics.generateSubshellList`/`Basics.generateCsfRs`/`Basis(...)` sequence
    `Basics.generateBasis` itself uses, so no new CSF machinery is written, only composed differently.

    ONE CONSISTENT ONE-ELECTRON BASIS, SO NO `BiOrthogonal` IS NEEDED. The pseudostate ladder is generated in
    the FROZEN local potential of `finalIonLevel` (nuclear + `Basics.DFSField(1.0)` mean field, the same
    combination `Continuum.generateOrbitalForLevel` already uses), via the SAME generalized eigenproblem
    (`Bsplines.diagonalizeLocalMatrix`) that an ordinary hydrogenic SCF start already solves -- the pseudostate
    orbitals are mutually orthogonal to each other by construction. They are explicitly Gram-Schmidt
    orthogonalized against any OCCUPIED bound orbital of the SAME kappa already in `finalIonLevel.basis`
    (`removeOccupiedComponent`, mirroring `RandomPhase.projectOntoUnoccupied`'s own, already-validated
    algorithm) -- a core electron cannot be excited into a filled orbital, and without this step a pseudostate
    whose kappa coincides with an occupied core shell would not be orthogonal to it. `checkFrozenOrbitalConsistency`
    then RAISES, rather than silently mixing inconsistent orbitals, if `initialLevel.basis` and `finalIonLevel.basis`
    disagree (overlap below `1.0 - 1.0e-6`) on any subshell they share -- the caller's responsibility is to have
    generated both multiplets against ONE consistent core, exactly as every frozen-orbital CI in JAC already
    assumes elsewhere.

    **THE ONE HONEST LIMITATION THIS BUYS, STATED RATHER THAN HIDDEN**: the pseudostate ladder is FROZEN, built
    once in a fixed local potential, never relaxed or EOL-optimized against the full combined CI. This is the
    same kind of approximation a DFS-based CI already makes elsewhere in JAC, not a new one invented here --
    but it IS an approximation, and a future round that compares against a relaxed treatment should expect some
    residual from it.

    THE NORMALIZATION: `crossSectionAt` turns the discrete "stick spectrum" of dipole strengths into a
    continuous cross section via the density-of-states recipe `sigma(E_k) = dipoleStrength_k / localSpacing_k`,
    MULTIPLIED by the SAME physical prefactor `PhotoIonization.computeCrossSection` applies to an energy-
    normalized continuum amplitude (`8*pi^3/(alpha*omega*gi)`) -- the prefactor is NOT optional and was missing
    from the first working version: Phase 0's toy validation only ever checked a RATIO (the Fano lineshape
    shape, and a bare box against its own defining formula), in which any constant prefactor cancels and so was
    never caught. `examples/example-Dz.jl` branch a's own PI-limit check is what found this, giving ratios of
    order 1e-4 before the fix and order 1 after it.

    **A SECOND LIMITATION, FOUND THE SAME WAY, NOT YET SOLVED: the pseudostate ladder shares the SAME grid the
    bound SCF used, with no mechanism to give it its OWN, independently-sized or independently-resolved box.**
    `buildPseudostateLadder` calls `Bsplines.generatePrimitives(grid)` on whatever `grid` the caller's bound
    orbitals were generated on. This is adequate near threshold but fails in two different, now-measured ways
    away from it (`examples/example-Dz.jl` branches b and c, 05-Nov-2026): a resonance far above threshold (Fe
    K-LL, ~171 Ha electron energy) sits where a box sized for Fe's own tiny bound orbitals (rbox ~ 4 a.u.) gives
    free-particle-like level spacing of order TEN Ha -- too coarse to hold even one pseudostate in a 5 Ha window
    there; and a genuinely NARROW resonance (He 2s2p, Gamma ~ 38 meV) sits inside a window whose LOCAL spacing,
    even on a more generous He-sized box (rbox = 25 a.u.), is still far coarser than Gamma itself, so an 8-Gamma-
    wide scan around it returns a single, perfectly flat stick with no visible lineshape at all. Both are the
    SAME root cause -- one shared grid asked to resolve a LOCAL feature far from where it was sized -- and
    fixing it properly means decoupling the pseudostate ladder's own box from the bound SCF's, which in turn
    means reconciling two different one-electron bases rather than the one this module currently assumes
    throughout. NOT attempted in this round; held at challenges-list S29.

    VALIDATED, before any `src/` code existed: two independent toy cases (see `work/photoionization-pseudostates-
    toy/`, Phase 0, 05-Nov-2026) -- a closed-form Fano two-level-plus-flat-continuum model (exact in the trivial
    uncoupled limit, converging toward the exact Fano lineshape as the box widens at fixed pseudostate density)
    and a direct check that `Bsplines.diagonalizeLocalMatrix`'s own positive-energy branch reproduces the exact
    free-particle box density of states to 1e-3..1e-4 relative accuracy, improving with box size. Both toys
    validated the SHAPE of the normalization, which is why neither caught the missing prefactor above.

    **A THIRD THING FOUND AND FIXED THE SAME WAY: `crossSectionAt` LINEARLY INTERPOLATES between the two
    sticks bracketing the requested energy, rather than reading off the single nearest one.** Measured against
    hydrogen's own EXACT Stobbe cross section (one-electron, no SCF or frozen-core approximation anywhere, so
    nothing else could be responsible) -- the first, nearest-neighbour version gave ratio = 1.0636 (6.4% high).
    Printing the individual stick-state sigma_k showed them SMOOTH and monotonically falling with energy, not
    noisy, so the error was never statistical scatter to be averaged away: it was evaluating at the nearest
    AVAILABLE stick's own energy rather than the one actually requested, on a cross section that falls off
    steeply here. Interpolating between the two bracketing sticks closed the gap to ratio = 1.0068 (0.7%).
    This is a DIFFERENT fix from the missing-prefactor one above, found the same way (a real example, not the
    toys) and independent of it. It does NOT touch the second limitation (the shared-grid one): a window with
    zero pseudostates in it, or one far coarser than Gamma, has nothing to interpolate between either.
"""
module PhotoIonizationPseudostates


using  Printf, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..Hamiltonian, ..ManyElectron, ..Nuclear,
       ..PhotoEmission, ..PhotoIonization, ..Radial, ..RadialIntegrals, ..TableStrings


"""
`struct  PhotoIonizationPseudostates.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of the diagonalize-once, pseudostate-continuum treatment
        of resonant photoionization.

    + multipoles               ::Array{EmMultipole,1}   ... Multipoles of the absorbed photon that are taken into account.
    + gauges                   ::Array{UseGauge,1}      ... Gauges to be included into the computations.
    + photonEnergies           ::Array{Float64,1}       ... Photon energies [a.u.] at which the SMOOTHED, density-of-states
                                                            normalized cross section is evaluated -- ONE diagonalization
                                                            per (initial level, final-ion level, kappa) serves the WHOLE
                                                            scan, unlike `PhotoIonizationInterference`'s per-energy solve.
    + lValues                  ::Array{Int64,1}         ... Orbital angular momenta of the outgoing electron that are
                                                            considered, exactly as in `PhotoIonization.Settings`.
    + pseudostateWindow        ::Tuple{Float64,Float64}  ... Energy window [a.u.], in the ejected electron's OWN one-
                                                            electron energy (i.e. directly comparable to photonEnergy minus
                                                            the ionization energy), of box-discretized pseudostates to
                                                            retain per kappa. Should comfortably bracket the electron-energy
                                                            range `photonEnergies` implies, with margin on both sides --
                                                            `crossSectionAt` is only trustworthy away from the window's own
                                                            edges, exactly as the Phase-0 toy validations were restricted.
    + continuumExchange        ::Union{Nothing,Bool}    ... True/false to add the frozen nonlocal Dirac-Fock exchange
                                                            operator (`Bsplines.setupExchangeMatrix`, item 99's mechanism)
                                                            to the pseudostate ladder's own generating potential;
                                                            `nothing` -- the default -- resolves to `false`. Unlike other
                                                            process modules this does NOT resolve automatically from a
                                                            `finalAsfSettings.scField`, since this entry point is handed a
                                                            `Level` rather than an `AsfSettings`.
    + printBefore              ::Bool                   ... True, if the combined-basis CI is printed before diagonalization.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                 ::Array{EmMultipole,1}
    gauges                      ::Array{UseGauge,1}
    photonEnergies               ::Array{Float64,1}
    lValues                       ::Array{Int64,1}
    pseudostateWindow              ::Tuple{Float64,Float64}
    continuumExchange                ::Union{Nothing,Bool}
    printBefore                        ::Bool
end


"""
`PhotoIonizationPseudostates.Settings()`
    ... constructor for an `empty` instance of PhotoIonizationPseudostates.Settings; a settings::Settings is returned.
"""
function Settings()
    Settings( EmMultipole[E1], UseGauge[UseCoulomb], Float64[], [0,1,2,3,4,5], (0., 1.), nothing, false )
end


"""
`PhotoIonizationPseudostates.Settings(set::PhotoIonizationPseudostates.Settings;`

        multipoles=..,          gauges=..,              photonEnergies=..,      lValues=..,
        pseudostateWindow=..,   continuumExchange=..,   printBefore=..)

    ... keyword copy-constructor for re-defining selected values of a settings::PhotoIonizationPseudostates.Settings; a
        settings::Settings is returned.
"""
function Settings(set::PhotoIonizationPseudostates.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,        gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,        lValues::Union{Nothing,Array{Int64,1}}=nothing,
    pseudostateWindow::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    continuumExchange::Union{Nothing,Bool}=nothing,                  printBefore::Union{Nothing,Bool}=nothing)

    if  isnothing(multipoles)          multipolesx         = set.multipoles         else   multipolesx         = multipoles         end
    if  isnothing(gauges)              gaugesx             = set.gauges             else   gaugesx             = gauges             end
    if  isnothing(photonEnergies)      photonEnergiesx     = set.photonEnergies     else   photonEnergiesx     = photonEnergies     end
    if  isnothing(lValues)             lValuesx            = set.lValues            else   lValuesx            = lValues            end
    if  isnothing(pseudostateWindow)   pseudostateWindowx  = set.pseudostateWindow  else   pseudostateWindowx  = pseudostateWindow  end
    if  isnothing(continuumExchange)   continuumExchangex  = set.continuumExchange  else   continuumExchangex  = continuumExchange  end
    if  isnothing(printBefore)         printBeforex        = set.printBefore        else   printBeforex        = printBefore        end

    Settings( multipolesx, gaugesx, photonEnergiesx, lValuesx, pseudostateWindowx, continuumExchangex, printBeforex )
end


# `Base.show(io::IO, settings::PhotoIonizationPseudostates.Settings)`  ... prepares a proper printout of settings.
function Base.show(io::IO, settings::PhotoIonizationPseudostates.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "lValues:                  $(settings.lValues)  ")
    println(io, "pseudostateWindow:        $(settings.pseudostateWindow)  ")
    println(io, "continuumExchange:        $(settings.continuumExchange)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
end


"""
`struct  PhotoIonizationPseudostates.StickState`
    ... defines a type for ONE eigenstate of the combined (bound + pseudostate-ladder) diagonalization.

    + energy                   ::Float64                ... Total energy of this eigenstate [a.u.].
    + dipoleStrength            ::EmProperty             ... |<this state|D|initialLevel>|^2, both gauges, summed
                                                            incoherently over the requested multipoles.
    + localSpacing               ::Float64                ... (E_{k+1}-E_{k-1})/2 -- the density-of-states normalization
                                                            divisor `crossSectionAt` uses.
    + weightOnReference            ::Float64                ... |c|^2 summed over the CSFs this state shares with
                                                            `initialLevel.basis` (i.e. everything that is NOT a newly
                                                            added pseudostate-ladder CSF) -- a resonant-vs-background
                                                            diagnostic: close to 1 marks a state dominated by the bound,
                                                            autoionizing character, close to 0 a smooth background state.
"""
struct StickState
    energy                      ::Float64
    dipoleStrength                ::EmProperty
    localSpacing                    ::Float64
    weightOnReference                 ::Float64
end


"""
`PhotoIonizationPseudostates.StickState()`
    ... constructor for an `empty` instance of PhotoIonizationPseudostates.StickState; a state::StickState is returned.
"""
function StickState()
    StickState( 0., EmProperty(0., 0.), 0., 0. )
end


# `Base.show(io::IO, state::PhotoIonizationPseudostates.StickState)`  ... prepares a proper printout of state.
function Base.show(io::IO, state::PhotoIonizationPseudostates.StickState)
    println(io, "energy:                   $(state.energy)  ")
    println(io, "dipoleStrength:           $(state.dipoleStrength)  ")
    println(io, "localSpacing:             $(state.localSpacing)  ")
    println(io, "weightOnReference:        $(state.weightOnReference)  ")
end


"""
`struct  PhotoIonizationPseudostates.Spectrum`
    ... defines a type for the full stick spectrum of the combined diagonalization, for ONE (initial level,
        final-ion level, kappa) combination.

    + kappa                     ::Int64                  ... Partial wave of the ejected electron this spectrum was
                                                            built for.
    + initialLevel               ::Level                  ... Initial (bound, possibly correlated) level.
    + finalIonLevel                ::Level                  ... Final-ion level, one electron fewer than initialLevel.
    + states                         ::Array{StickState,1}   ... Every eigenstate of the combined diagonalization, sorted
                                                            by energy.
"""
struct Spectrum
    kappa                        ::Int64
    initialLevel                   ::Level
    finalIonLevel                     ::Level
    states                               ::Array{StickState,1}
end


"""
`PhotoIonizationPseudostates.Spectrum()`
    ... constructor for an `empty` instance of PhotoIonizationPseudostates.Spectrum; a spectrum::Spectrum is returned.
"""
function Spectrum()
    Spectrum( 0, Level(), Level(), StickState[] )
end


# `Base.show(io::IO, spectrum::PhotoIonizationPseudostates.Spectrum)`  ... prepares a proper printout of spectrum.
function Base.show(io::IO, spectrum::PhotoIonizationPseudostates.Spectrum)
    println(io, "kappa:                    $(spectrum.kappa)  ")
    println(io, "initialLevel:             $(spectrum.initialLevel)  ")
    println(io, "finalIonLevel:            $(spectrum.finalIonLevel)  ")
    println(io, "states:                   $(length(spectrum.states)) stick state(s)  ")
end


"""
`PhotoIonizationPseudostates.crossSectionAt(spectrum::PhotoIonizationPseudostates.Spectrum, omega::Float64)`
    ... to evaluate the density-of-states-normalized cross section at photon energy omega [a.u.]. Each stick
        state k already carries `sigma_k = dipoleStrength_k / localSpacing_k`, the normalization this module's
        own docstring describes and Phase 0 validated; this function then LINEARLY INTERPOLATES sigma_k between
        the two sticks that bracket `spectrum.initialLevel.energy + omega`, rather than reading off the single
        nearest one. A value::EmProperty is returned.

        **WHY INTERPOLATION, NOT NEAREST-NEIGHBOUR -- FOUND 05-Nov-2026, read before reverting this.** The
        first working version used nearest-neighbour and gave a stable ~6-8% discrepancy against hydrogen's
        EXACT Stobbe cross section (`examples/example-Dz.jl` branch a). Printing the individual sigma_k around
        the evaluation point showed them to be SMOOTH and monotonic, not noisy -- so the discrepancy was never
        statistical scatter to be averaged away, it was simply evaluating at the nearest AVAILABLE stick's own
        energy instead of the requested one, on a function (hydrogen's photoionization cross section) that
        falls off steeply with energy. A plain linear interpolation between the two bracketing sticks closed
        the hydrogen residual from 6.4% to under 1%. At the array's own edges, where no upper or lower
        bracketing stick exists, this falls back to the nearest endpoint.
"""
function crossSectionAt(spectrum::PhotoIonizationPseudostates.Spectrum, omega::Float64)
    length(spectrum.states) == 0   &&   error("PhotoIonizationPseudostates.crossSectionAt(): spectrum carries no stick states.")
    totalEnergy = spectrum.initialLevel.energy + omega
    states = spectrum.states                                    ## already sorted by energy (computeSpectrum)
    sigmaOf(st) = EmProperty(st.dipoleStrength.Coulomb/st.localSpacing, st.dipoleStrength.Babushkin/st.localSpacing)

    if      totalEnergy <= states[1].energy       return( sigmaOf(states[1]) )
    elseif  totalEnergy >= states[end].energy     return( sigmaOf(states[end]) )
    end
    hi = findfirst(st -> st.energy >= totalEnergy, states);   lo = hi - 1
    frac = (totalEnergy - states[lo].energy) / (states[hi].energy - states[lo].energy)
    sigmaLo = sigmaOf(states[lo]);   sigmaHi = sigmaOf(states[hi])

    return( EmProperty( sigmaLo.Coulomb   + frac*(sigmaHi.Coulomb   - sigmaLo.Coulomb),
                        sigmaLo.Babushkin + frac*(sigmaHi.Babushkin - sigmaLo.Babushkin) ) )
end


"""
`PhotoIonizationPseudostates.removeOccupiedComponent(orb::Radial.Orbital, occupiedOrbitals::Dict{Subshell,Orbital},
                                                      grid::Radial.Grid)`
    ... removes from orb every component along an OCCUPIED, BOUND orbital of orb's own kappa (mirroring
        `RandomPhase.projectOntoUnoccupied`'s own, already-validated algorithm, written locally rather than
        reused since that function's input type is `DalgarnoLewis.PerturbedOrbital`, not a plain
        `Radial.Orbital`), and renormalizes the result to unit norm in this grid's quadrature. A core electron
        cannot be excited into a filled orbital, so without this step a pseudostate whose kappa coincides with
        an occupied ion-core shell would not be orthogonal to it. An orbital::Orbital is returned.
"""
function removeOccupiedComponent(orb::Radial.Orbital, occupiedOrbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)
    P = copy(orb.P);   Q = copy(orb.Q)
    for  (sh, occ)  in  occupiedOrbitals
        (sh.kappa == orb.subshell.kappa  &&  occ.energy < 0.0)   ||   continue
        ov = RadialIntegrals.overlap(occ.P, P, grid) + RadialIntegrals.overlap(occ.Q, Q, grid)
        nn = RadialIntegrals.overlap(occ.P, occ.P, grid) + RadialIntegrals.overlap(occ.Q, occ.Q, grid)
        nn <= 0.0   &&   continue
        for  i = 1:min(length(P), length(occ.P))   P[i] = P[i] - ov/nn * occ.P[i]   end
        for  i = 1:min(length(Q), length(occ.Q))   Q[i] = Q[i] - ov/nn * occ.Q[i]   end
    end
    norm = RadialIntegrals.overlap(P, P, grid) + RadialIntegrals.overlap(Q, Q, grid)
    norm <= 0.0   &&   error("PhotoIonizationPseudostates.removeOccupiedComponent(): a pseudostate was entirely " *
                              "consumed by the occupied-orbital projection -- the box or the window is too small.")
    P = P / sqrt(norm);   Q = Q / sqrt(norm)

    return( Orbital(orb.subshell, orb.isBound, orb.useStandardGrid, orb.energy, P, Q, orb.Pprime, orb.Qprime, orb.grid) )
end


"""
`PhotoIonizationPseudostates.buildPseudostateLadder(kappa::Int64, finalIonLevel::Level, nm::Nuclear.Model,
                                                      grid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings)`
    ... builds the box-discretized continuum pseudostate ladder of the given kappa, in the frozen local
        potential of finalIonLevel (nuclear + `Basics.DFSField(1.0)` mean field, the same combination
        `Continuum.generateOrbitalForLevel` already uses), keeping every eigenstate of
        `Bsplines.diagonalizeLocalMatrix`'s positive-energy branch that falls inside settings.pseudostateWindow.
        Each retained state is labelled `Subshell(1000+k, kappa)` (k = 1, 2, ..., purely a bookkeeping label,
        chosen far above any physical principal quantum number so it can never collide with an ordinary bound
        subshell) and Gram-Schmidt orthogonalized against finalIonLevel's own occupied orbitals of the same
        kappa via `removeOccupiedComponent`. A `ladder::Dict{Subshell,Orbital}` is returned.
"""
function buildPseudostateLadder(kappa::Int64, finalIonLevel::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::PhotoIonizationPseudostates.Settings)
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    meanPot    = Basics.computePotential(Basics.DFSField(1.0), grid, finalIonLevel)
    pot        = Basics.add(nuclearPot, meanPot)

    primitives = Bsplines.generatePrimitives(grid)
    storage    = Dict{String,Array{Float64,2}}()
    nsL        = primitives.grid.nsL;   nsS = primitives.grid.nsS
    overlap    = zeros(nsL+nsS, nsL+nsS)
    overlap[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    overlap[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    localMatrix = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
    if  settings.continuumExchange == true
        coreOccupations = Basics.extractMeanOccupation(finalIonLevel.basis)
        localMatrix      = localMatrix + Bsplines.setupExchangeMatrix(kappa, primitives, finalIonLevel.basis.orbitals,
                                                                       coreOccupations)
    end
    eig = Bsplines.diagonalizeLocalMatrix(kappa, localMatrix, overlap, primitives)

    windowLo, windowHi = settings.pseudostateWindow
    idxWindow = findall(E -> windowLo <= E <= windowHi, eig.values)
    length(idxWindow) < 3   &&   error("PhotoIonizationPseudostates.buildPseudostateLadder(): only " *
                                        "$(length(idxWindow)) pseudostate(s) fall inside pseudostateWindow = " *
                                        "$(settings.pseudostateWindow) for kappa = $kappa -- widen the window " *
                                        "or check that it is given in the ejected electron's OWN energy units.")

    ladder = Dict{Subshell,Orbital}()
    mtp    = primitives.grid.NoPoints
    for  (k, idx)  in  enumerate(idxWindow)
        sh  = Subshell(1000+k, kappa)
        raw = Bsplines.generateOrbitalFromPrimitives(sh, eig.values[idx], mtp, eig.vectors[idx], primitives)
        ladder[sh] = PhotoIonizationPseudostates.removeOccupiedComponent(raw, finalIonLevel.basis.orbitals, grid)
    end

    return( ladder )
end


"""
`PhotoIonizationPseudostates.checkFrozenOrbitalConsistency(initialLevel::Level, finalIonLevel::Level, grid::Radial.Grid)`
    ... RAISES, rather than silently mixing inconsistent orbitals, if initialLevel.basis and finalIonLevel.basis
        disagree -- overlap below 1.0 - 1.0e-6 -- on any subshell they share. The combined CI this module builds
        assumes ordinary (non-BiOrthogonal) Slater-Condon rules, which require every orbital shared by the two
        sides to be LITERALLY the same radial function; the caller's responsibility is to have generated both
        multiplets against one consistent frozen core. Nothing is returned on success.
"""
function checkFrozenOrbitalConsistency(initialLevel::Level, finalIonLevel::Level, grid::Radial.Grid)
    shared = intersect(keys(initialLevel.basis.orbitals), keys(finalIonLevel.basis.orbitals))
    for  sh  in  shared
        oi = initialLevel.basis.orbitals[sh];    of = finalIonLevel.basis.orbitals[sh]
        ov = RadialIntegrals.overlap(oi.P, of.P, grid) + RadialIntegrals.overlap(oi.Q, of.Q, grid)
        if  ov < 1.0 - 1.0e-6
            error("PhotoIonizationPseudostates.checkFrozenOrbitalConsistency(): initialLevel and finalIonLevel " *
                  "disagree on subshell $sh -- overlap = $ov, expected ~1.0.  The combined CI needs ONE " *
                  "consistent frozen core; re-generate both multiplets against the same orbitals.")
        end
    end

    return( nothing )
end


"""
`PhotoIonizationPseudostates.computeSpectrum(initialLevel::Level, finalIonLevel::Level, kappa::Int64,
        symmetries::Array{LevelSymmetry,1}, nm::Nuclear.Model, grid::Radial.Grid,
        settings::PhotoIonizationPseudostates.Settings)`
    ... builds the pseudostate ladder of the given kappa (`buildPseudostateLadder`), merges it with
        initialLevel's own CSFs into ONE combined basis (re-expanding initialLevel.basis onto the extended
        subshell list via `ManyElectron.Basis(basis, subshells)`'s existing zero-padding constructor, then
        appending the new pseudostate CSFs, both restricted to `symmetries`), diagonalizes it ONCE via
        `Hamiltonian.performCI`, and computes the dipole amplitude from initialLevel onto every resulting
        eigenstate via `PhotoEmission.amplitude(Absorption(), ...)` -- which already tolerates the two levels
        carrying different subshell lists/orbital dicts, so no further basis surgery is needed here. A
        spectrum::Spectrum is returned.
"""
function computeSpectrum(initialLevel::Level, finalIonLevel::Level, kappa::Int64, symmetries::Array{LevelSymmetry,1},
                         nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings)
    PhotoIonizationPseudostates.checkFrozenOrbitalConsistency(initialLevel, finalIonLevel, grid)
    ladder = PhotoIonizationPseudostates.buildPseudostateLadder(kappa, finalIonLevel, nm, grid, settings)

    # The final ion's own DOMINANT CSF supplies the core occupation every pseudostate configuration shares.
    dominantIdx = argmax( abs.(finalIonLevel.mc) )
    dominantW2  = finalIonLevel.mc[dominantIdx]^2
    dominantW2 < 0.9   &&
        @warn("PhotoIonizationPseudostates.computeSpectrum(): finalIonLevel's dominant CSF carries only " *
              "$dominantW2 of the weight -- the single-CSF ion-core picture the pseudostate ladder assumes " *
              "is a poorer approximation here than usual.")
    coreCsf          = finalIonLevel.basis.csfs[dominantIdx]
    coreOccupation   = Dict{Subshell,Int64}( finalIonLevel.basis.subshells[i] => coreCsf.occupation[i]
                                              for i = 1:length(finalIonLevel.basis.subshells)
                                              if  coreCsf.occupation[i] > 0 )

    # MUST be in Base.isless(::Subshell,...)'s own canonical (n, then |kappa|, then kappa) order -- not a bare
    # vcat/Dict-key order -- since PhotoEmission.amplitude's own internal Basics.merge(...) produces exactly
    # that sorted order and then requires THIS list to appear, in order, as a subsequence of it.
    extendedSubshells = Base.sort( vcat(initialLevel.basis.subshells, collect(keys(ladder))), lt=Base.isless )
    paddedBasis       = ManyElectron.Basis(initialLevel.basis, extendedSubshells)
    boundCsfs         = CsfR[];   for  csf in paddedBasis.csfs
                                      LevelSymmetry(csf.J, csf.parity) in symmetries   &&   push!(boundCsfs, csf)   end

    pseudoRelConfs = ConfigurationR[]
    for  sh  in  keys(ladder)
        subshells = merge( coreOccupation, Dict(sh => 1) )
        push!( pseudoRelConfs, ConfigurationR(subshells, finalIonLevel.basis.NoElectrons + 1) )
    end
    pseudoCsfs = CsfR[]
    for  relconf in pseudoRelConfs
        newCsfs = Basics.generateCsfRs(relconf, extendedSubshells)
        for  csf in newCsfs    LevelSymmetry(csf.J, csf.parity) in symmetries   &&   push!(pseudoCsfs, csf)    end
    end
    length(pseudoCsfs) == 0   &&   error("PhotoIonizationPseudostates.computeSpectrum(): no pseudostate CSF of " *
                                          "kappa = $kappa reaches any of $symmetries from this ion core.")

    combinedCsfs      = vcat(boundCsfs, pseudoCsfs);   nBound = length(boundCsfs)
    combinedOrbitals  = merge(initialLevel.basis.orbitals, ladder)
    coreSubshellList  = Subshell[]
    for  (k, sh)  in  enumerate(extendedSubshells)
        mocc = Basics.subshell_2j(sh) + 1;    isFilled = true
        for  csf in combinedCsfs    if  csf.occupation[k] != mocc    isFilled = false;   break   end   end
        isFilled   &&   push!(coreSubshellList, sh)
    end
    combinedBasis = Basis(true, initialLevel.basis.NoElectrons, extendedSubshells, combinedCsfs, coreSubshellList,
                          combinedOrbitals)

    mp = Hamiltonian.performCI(combinedBasis, nm, grid, AsfSettings(); printout=settings.printBefore, writeSummary=false)

    sortedLevels = sort(mp.levels, by = lev -> lev.energy)
    gi = Basics.twice(initialLevel.J) + 1.0
    states = StickState[]
    for  (n, lev)  in  enumerate(sortedLevels)
        weightOnReference = sum( lev.mc[1:nBound] .^ 2 )
        omega = lev.energy - initialLevel.energy
        dipoleStrength = EmProperty(0., 0.)
        for  mp_ in settings.multipoles
            if  string(mp_)[1] == 'E'
                aC = PhotoEmission.amplitude(Absorption(), mp_, Basics.Coulomb,   omega, lev, initialLevel, grid; printout=false)
                aB = PhotoEmission.amplitude(Absorption(), mp_, Basics.Babushkin, omega, lev, initialLevel, grid; printout=false)
                dipoleStrength = dipoleStrength + EmProperty(abs2(aC), abs2(aB))
            else
                aM = PhotoEmission.amplitude(Absorption(), mp_, Basics.Magnetic, omega, lev, initialLevel, grid; printout=false)
                dipoleStrength = dipoleStrength + EmProperty(abs2(aM), abs2(aM))
            end
        end
        # THE SAME PHYSICAL PREFACTOR PhotoIonization.computeCrossSection APPLIES TO AN ENERGY-NORMALIZED
        # continuum amplitude -- 8*pi^3/(alpha*omega*gi) -- converts a squared amplitude into a cross section;
        # it is MISSING here on purpose until this point, since nothing before it needs it, and omega varies
        # per stick state so it cannot be folded into a single constant. Found missing (05-Nov-2026) when
        # branch (a)'s own PI-limit check came out ~10^4 too small -- Phase 0's toy validated only the SHAPE
        # (a ratio, in which any missing constant factor cancels), never the absolute scale.
        if  omega > 0.
            csFactor = 8 * pi^3 / Defaults.getDefaults("alpha") / omega / gi
            dipoleStrength = EmProperty(csFactor * dipoleStrength.Coulomb, csFactor * dipoleStrength.Babushkin)
        end
        if      n == 1               localSpacing = sortedLevels[2].energy - sortedLevels[1].energy
        elseif  n == length(sortedLevels)   localSpacing = sortedLevels[end].energy - sortedLevels[end-1].energy
        else                           localSpacing = (sortedLevels[n+1].energy - sortedLevels[n-1].energy) / 2
        end
        push!( states, StickState(lev.energy, dipoleStrength, localSpacing, weightOnReference) )
    end

    return( Spectrum(kappa, initialLevel, finalIonLevel, states) )
end


"""
`PhotoIonizationPseudostates.computeSpectra(finalMultiplet::Multiplet, initialLevel::Level, nm::Nuclear.Model,
        grid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings; output::Bool=true)`
    ... the standard entry point of this module: for every final-ion level of finalMultiplet and every kappa
        that `PhotoIonization.determineChannels` finds reachable from initialLevel under settings.multipoles/
        lValues, builds and diagonalizes the combined CI (`computeSpectrum`) and collects the resulting stick
        spectra. A list of spectra::Array{PhotoIonizationPseudostates.Spectrum,1} is returned if output = true,
        and nothing otherwise.

        TAKES ONE initialLevel, NOT A WHOLE initialMultiplet -- deliberately, unlike `PhotoIonizationInterference
        .computePathways`. Looping over every level of a multiplet that may ALSO contain the autoionizing
        levels (see this module's own docstring) would have no way to tell a genuine below-threshold initial
        state apart from one of the higher-lying autoionizing ones; the caller already knows which level of
        their own multiplet is the one to ionize FROM, and passes it directly.
"""
function computeSpectra(finalMultiplet::Multiplet, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid,
                        settings::PhotoIonizationPseudostates.Settings; output::Bool=true)
    println("")
    printstyled("PhotoIonizationPseudostates.computeSpectra(): the diagonalize-once pseudostate-continuum " *
                "computation starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------- \n",
                color=:light_green)
    println("")
    piSettings = PhotoIonization.Settings(PhotoIonization.Settings(); multipoles=settings.multipoles,
                                          gauges=settings.gauges, lValues=settings.lValues)
    spectra = PhotoIonizationPseudostates.Spectrum[]
    for  finalIonLevel  in  finalMultiplet.levels
        partialWaves = PhotoIonization.determineChannels(finalIonLevel, initialLevel, piSettings)
        for  pw  in  partialWaves
            symmetries = [ch.symmetry for ch in pw.channels]
            length(symmetries) == 0   &&   continue
            push!( spectra, PhotoIonizationPseudostates.computeSpectrum(initialLevel, finalIonLevel, pw.kappa,
                                                                        symmetries, nm, grid, settings) )
        end
    end
    PhotoIonizationPseudostates.displayResults(stdout, spectra, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    PhotoIonizationPseudostates.displayResults(iostream, spectra, settings)    end

    if    output    return( spectra )
    else            return( nothing )
    end
end


"""
`PhotoIonizationPseudostates.displayResults(stream::IO, spectra::Array{PhotoIonizationPseudostates.Spectrum,1},
        settings::PhotoIonizationPseudostates.Settings)`
    ... to display, for every spectrum and every requested photon energy, the density-of-states normalized
        cross section from `crossSectionAt`. A neat table is printed but nothing is returned otherwise.
"""
function displayResults(stream::IO, spectra::Array{PhotoIonizationPseudostates.Spectrum,1},
                        settings::PhotoIonizationPseudostates.Settings)
    nx = 100
    println(stream, " ")
    println(stream, "  Diagonalize-once pseudostate-continuum photoionization:  cross sections [a.u.], Coulomb gauge:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(10, "i-level"; na=2)
    sa = sa * TableStrings.center(10, "f-level"; na=2)
    sa = sa * TableStrings.center(8,  "kappa"; na=2)
    sa = sa * TableStrings.center(16, "E_photon"; na=2)
    sa = sa * TableStrings.center(16, "sigma (Coulomb)"; na=2)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  spectrum in spectra,  omega in settings.photonEnergies
        omegaAu = Defaults.convertUnits("energy: to atomic", omega)
        sigma   = PhotoIonizationPseudostates.crossSectionAt(spectrum, omegaAu)
        sa = "  " * TableStrings.center(10, string(spectrum.initialLevel.index); na=2)
        sa = sa * TableStrings.center(10, string(spectrum.finalIonLevel.index); na=2)
        sa = sa * TableStrings.center(8,  string(spectrum.kappa); na=2)
        sa = sa * @sprintf("%14.6e", omega) * "  "
        sa = sa * @sprintf("%14.6e", sigma.Coulomb)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")

    return( nothing )
end

end # module
