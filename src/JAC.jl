#
#  Use git:                    git status ....;   git add <filenames>;   git commit -m "..";   git push;   git rm  <filenames>
#  Use Jupyter notebooks:      using IJulia;   notebook()
#  Activation:                 ];   pkg> up;   pkg> activate
#  Working with JAC:           using Revise;   using JAC;   include("../src/jac.jl");   pkg> test
#
"""
`module JAC`  
    ... Jena Atomic Calculator (JAC) provides tools for performing atomic (structure) calculations at various degrees of complexity 
        and sophistication. It has been designed to not only calculate atomic level structures and properties [such as g-factors or
        hyperfine and isotope-shift parameters] but also transition amplitudes between bound-state levels [for the anapole moment, dipole 
        operator, electron electric-dipole moment, parity non-conservation, etc.] and, in particular, (atomic) transition probabilities, 
        Auger rates, photoionization cross sections, radiative and dielectronic recombination rates as well as cross sections for many 
        other (elementary) processes. JAC also facilitates interactive computations, the simulation of atomic cascades, the time-evolution 
        of statistical tensors, a few semi-empirical estimates of atomic properties as well as the simplification of symbolic expressions
        from Racah's algebra. -- In addition, the JAC module supports the display of level energies, electron and photon spectra, 
        radial orbitals and and other atomic data.


**`Perform (atomic) computations of different complexity:`**  
    JAC will eventually support **ten kinds** of computations which can be summarized as follows:

   + Atomic computations, based on explicitly specified electron configurations.
   + Restricted active-space computations (RAS; partly implemented).
   + Interactive computations.
   + Atomic cascade computations (partly implemented).
   + Atomic representations (Green and close-coupling functions, complex rotation; not yet implemented).
   + Atomic responses (not yet implemented).
   + Atomic descriptors for machine learning algorithms (not yet implemented).
   + Time-evolution of statistical tensors in (intense) light pusles (not yet implemented).
   + Semi-empirical estimates of cross sections, etc. (partly implemented).
   + Symbolic evaluation of expressions from Racah's algebra, etc. (partly implemented).


**`Further details and information`**

        + Kinds of atomic implementation                                       [cf. ? Details.kindsOfComputation]
        + Atomic amplitudes (partly) implemented in JAC                        [cf. ? Details.amplitudes]
        + Atomic level properties (partly) implemented in JAC                  [cf. ? Details.properties]
        + Atomic processes (partly) implemented in JAC                         [cf. ? Details.processes]
        + Interactive use of JAC procedures                                    [cf. ? Details.interactive]
        + Design principles and limitations of the JAC program                 [cf. ? Details.design]
        + Data types, structs and name conventions of the JAC module           [cf. ? Details.datatypes]
        + Atomic cascade computations and approximations                       [cf. ? Details.decayCascades]
        + Use of (em) light pulses in the time evolution of statist. tensors   [cf. ? Details.pulses]
        + Why Julia ?                                                          [cf. ? Details.whyJulia]

"""
module JAC

using  Dates, Printf,  LinearAlgebra, IJulia, SpecialFunctions, FortranFiles, GaussQuadrature, QuadGK, GSL, JLD, SymEngine  ## , Interact

##x export @racahsum, 
export AbstractQedModel, add, analyze, AlphaX, AlphaVariation, AnapoleMoment, AngularJ64, AngularM64, AngularJ, AngularMomentum, 
       AsfSettings, Atomic, AtomicStructure, Auger, AutoIonization, 
       Basics, Basis, 
       CartesianVector, CiSettings, CiExpansion, CloseCoupling, compute, convertUnits, Cascade, Configuration, ConfigurationR, 
       Continuum, CsfR, CoulombExcitation, CoulombIonization,  
       diagonalize, Defaults, DecayYield, Details, Dielectronic, DoubleAuger,
       estimate, ElectricDipoleMoment, Einstein, EinsteinX, EmMultipole, evaluate,
       E1, M1, E2, M2, E3, M3, E4, M4,
       FormFactor, FormF, 
       generate, GreenSettings, GreenChannel, GreenExpansion, getDefaults, Green, Gui,
       Hfs, HighHarmonic, HFS, HydrogenicIon,
       interpolate, integrate, ImpactExcitation, ImpactExcitationAutoion, ImpactIonization, InternalConversion, Isotope, IsotopeShift, 
       Kronecker,
       LandeZeeman, LandeJ,  LandeF, Level, LevelSymmetry, LSjj, LSjjSettings,
       ManyElectron, MeanFieldSettings, MeanFieldBasis, Model, modify, MultiPhotonDeExcitation, MultiPhotonDoubleIon, 
       MultiPhotonIonization, MultipoleMoment, MultipolePolarizibility, Multiplet, 
       NoAmplitude, NoProcess, Nuclear, NoneQed, NoProperty, 
       Orbital, 
       perform, provide, PairAnnihilation1Photon, PairAnnihilation2Photon, PairProduction, ParityNonConservation, PeriodicTable,
       PhotoEmission, PhotoExcitation, PhotoExcitationAutoion, PhotoExcitationFluores, PhotoIonization, PhotoIonizationFluores, 
       PhotoIonizationAutoion, PhotoRecombination, PlasmaShift, Plasma, Polarity, 
       QedPetersburg, QedSydney,
       RacahAlgebra, RacahExpression, Radial, RadialIntegrals, Radiative, RadiativeAuger, RasSettings, RasStep, 
       RasExpansion, RayleighCompton, recast, REDA, READI, Representation,
       SchiffMoment, setDefaults, Shell, Spectroscopy, Subshell,
       tabulate, TestFrames, Triangle, tools,
       UseCoulomb, UseBabushkin, UseGauge,
       WeightedCartesian, W3j, W6j, W9j,
       Yields, 
       Zeeman
     
# Basic data and data structures
include("module-Basics.jl");        using ..Basics
include("module-Radial.jl");        using ..Radial
include("module-Math.jl");          using ..Math
include("module-Defaults.jl");      using ..Defaults
include("module-ManyElectron.jl");  using ..ManyElectron
include("module-Nuclear.jl");       using ..Nuclear
include("module-SpinAngular.jl");   using ..SpinAngular
include("module-BiOrthogonal.jl");  using ..BiOrthogonal

# Specialized functions/methods to manipulate these data
include("module-AngularMomentum.jl")
include("module-AngularCoefficients-Ratip2013.jl")
include("module-Bsplines.jl")
include("module-Continuum.jl")
include("module-ConfigurationSpace.jl")
include("module-Details.jl")
include("module-RadialIntegrals.jl")
include("module-HydrogenicIon.jl")
include("module-InteractionStrength.jl")
include("module-InteractionStrengthQED.jl")
include("module-PeriodicTable.jl")
include("module-TableStrings.jl")
include("module-Tools.jl")
include("module-LSjj.jl")

# Functions/methods for atomic amplitudes
include("module-MultipoleMoment.jl")
include("module-ParityNonConservation.jl") 

include("module-PhotoEmission.jl")

# Functions/methods for atomic properties
include("module-Einstein.jl")    
include("module-Hfs.jl")
include("module-IsotopeShift.jl")
include("module-LandeZeeman.jl")
include("module-AlphaVariation.jl")
include("module-FormFactor.jl")
include("module-DecayYield.jl")
##x include("module-GreenFunction.jl")
include("module-CloseCoupling.jl")
include("module-MultipolePolarizibility.jl")
include("module-PlasmaShift.jl")

# Functions/methods for atomic processes
include("module-PhotoExcitation.jl")
include("module-PhotoIonization.jl")
include("module-PhotoRecombination.jl")
include("module-AutoIonization.jl")
include("module-Dielectronic.jl")
include("module-PhotoExcitationFluores.jl")
include("module-PhotoExcitationAutoion.jl")
include("module-RayleighCompton.jl")
include("module-MultiPhotonDeExcitation.jl")
include("module-CoulombExcitation.jl")
include("module-PhotoIonizationFluores.jl")
include("module-PhotoIonizationAutoion.jl")
include("module-CoulombIonization.jl")
include("module-ImpactExcitation.jl")
include("module-ImpactExcitationAutoion.jl")
include("module-RadiativeAuger.jl")
include("module-MultiPhotonIonization.jl")
include("module-MultiPhotonDoubleIon.jl")
include("module-InternalConversion.jl")      ## up to here + up to jac-tools
#= Further processes, not yet included into the code
include("module-DoubleAuger.jl")
include("module-REDA.jl")
include("module-READI.jl")
include("module-PairProduction.jl")
include("module-PairAnnihilation1Photon.jl")
include("module-PairAnnihilation2Photon.jl")

# Functions/methods for the computation of atomic responses
include("module-HighHarmonic.jl")  =#

# Functions/methods for semi-empirical estimations
# include("module-ImpactIonization.jl")
include("module-Semiempirical.jl")

# Functions/methods for atomic and cascade computations
include("module-Atomic.jl");        using ..Atomic
include("module-Cascade.jl");       using ..Cascade

# Functions/methods for atomic time evolutions
# include("module-Pulse.jl")
# include("module-Statistical.jl")

# Functions/methods for symbolic computations
include("module-RacahAlgebra.jl");  using ..RacahAlgebra

# Basic functions/methods to manipulate these data
include("module-BasicsAG.jl")
include("module-BasicsCompute.jl")
include("module-BasicsGenerate.jl")
include("module-BasicsHP.jl")
include("module-BasicsPerform.jl")
include("module-BasicsQZ.jl")
include("module-ManyElectronAZ.jl")

# Specialized macros
##x include("macro-racahsum.jl")

# All test functions/methods stay with the JAC root module
include("module-TestFrames.jl");  using ..TestFrames
    
function __init__()
    # The following variables need to be initialized at runtime to enable precompilation
    global JAC_SUMMARY_IOSTREAM = stdout
    global JAC_TEST_IOSTREAM    = stdout
end

println("\nWelcome to JAC:  A community approach to the computation of atomic structures, cascades and time evolutions " * 
        "[(C) Copyright by Stephan Fritzsche, Jena (2018-2020)].")

end


