#
#  Use git:                    git status ....;   git add <filenames>;   git commit -m "..";   git push;   git rm  <filenames>
#  Use Jupyter notebooks:      using IJulia;   notebook()
#  Activation:                 ];   pkg> up;   pkg> activate
#  Working with JAC:           using Revise;   using JAC;   include("../src/jac.jl");   pkg> test
#  
#  Copy to desktop             scp -r JAC.jl/ fritzsch@10.140.119.236:~/fri/.
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
        radial orbitals and other atomic data.


**`Perform (atomic) computations of different complexity:`**  
    JAC will eventually support **ten kinds** of computations which can be summarized as follows:

   + Atomic computations, based on explicitly specified electron configurations.
   + Restricted active-space computations (RAS).
   + Interactive computations.
   + Atomic cascade computations (partly implemented).
   + Atomic representations (Green and close-coupling functions, complex rotation; partly implemented).
   + Atomic responses (partly implemented).
   + Atomic descriptors for machine learning algorithms (not yet implemented).
   + Time-evolution of statistical tensors in (intense) light pusles (not yet implemented).
   + Semi-empirical estimates of cross sections, etc. (partly implemented).
   + Symbolic evaluation of expressions from Racah's algebra, etc.


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

using  Dates, Printf,  LinearAlgebra, IJulia, SpecialFunctions, FortranFiles, QuadGK, GSL, JLD2, SymEngine, 
       HypergeometricFunctions  ## , Interact, GaussQuadrature

export AbstractConfigurationRestriction, AbstractEeInteraction, AbstractCImethod, AbstractPotential, AbstractQedModel, AbstractStartOrbitals,
       AbstractProcessSettings, AbstractEmpiricalSettings, AbstractPropertySettings, 
       add, analyze, AlphaX, AlphaVariation, AnapoleMoment, 
       AngularJ64, AngularM64, AngularJ, AngularMomentum, 
       AsfSettings, Atomic, AtomicState, AtomicStructure, Auger, AugerInPlasma, AutoIonization, AverageAtom,
       Basics, Basis, BreitInteraction, Bsplines,
       CartesianVector, CiSettings, CiExpansion, CloseCoupling, compute, convertUnits, Cascade, Compton, Configuration, ConfigurationR, 
       Continuum, CsfR, Coulex, CoulombExcitation, Coulion, CoulombIonization, CoulombBreit, CoulombInteraction, ClebschGordan, 
       CorePolarization,
       diagonalize, Defaults, DecayYield, Details, Dielectronic, Dierec, Djpq, DoubleAutoIonization, DoubleAuger,
       Eimex, ElectronCapture, ElecCapture, estimate, ElectricDipoleMoment, Einstein, EinsteinX, EmMultipole, evaluate, ExpStokes, Empirical,
       E1, M1, E2, M2, E3, M3, E4, M4,
       FormFactor, FormF, FullCIeigen,
       generate, GreenSettings, GreenChannel, GreenExpansion, getDefaults, Green, Gui,
       Hfs, HighHarmonic, HFS, HydrogenicIon,
       interpolate, integrate, Integral, ImpactExcAuto, ImpactExcitation, ImpactExcitationAutoion, ImpactIonization, 
       InteractionStrength, InternalConv, InternalConversion, Isotope, IsotopeShift, 
       Kronecker,
       LandeZeeman, LandeJ,  LandeF, Level, LevelSymmetry, LSjj, LSjjSettings, LevelSelection, LineSelection, 
       ManyElectron, MeanFieldSettings, MeanFieldBasis, minus, Model, modify, MultiPhotonDE, MultiPhotonDeExcitation, MultiPhotonDoubleIon, 
       MultiPI, MultiPDI, MultiPhotonIonization, MultipoleMoment, MultipolePolarizibility, Multiplet, 
       NoAmplitude, NoProcess, Nuclear, NoneQed, NoProperty, 
       OneElectronSettings, OneElectronSpectrum, Orbital, 
       PathwaySelection, perform, provide, PairA1P, PairAnnihilation1Photon, PairAnnihilation2Photon, 
       PairProduction, ParityNonConservation, PeriodicTable, Parity, plus,
       Photo, PhotoExc, PhotoExcAuto, PhotoExcFluor, PhotoEmission, PhotoExcitation, PhotoExcitationAutoion, PhotoExcitationFluores, 
       PhotoIonAuto, PhotoIonFluor, PhotoDouble, PhotoIonization, PhotoDoubleIonization, PhotoIonizationFluores, PhotoIonizationAutoion, 
       PhotoRecombination, PlasmaShift, Plasma, Polarity, Pulse,
       QedPetersburg, QedSydney,
       RacahAlgebra, RacahExpression, Radial, RadialIntegrals, Radiative, RadiativeAuger, RAuger, RasSettings, RasStep, 
       RasExpansion, RayleighCompton, recast, Rec, REDA, READI, Representation, ReducedDensityMatrix, RadiativeOpacity,
       RestrictNoElectronsTo, RestrictParity, RestrictToShellDoubles, RequestMinimumOccupation, RequestMaximumOccupation,
       SchiffMoment, Semiempirical, setDefaults, Shell, SolidAngle, Spectroscopy, SpinAngular, SphericalTensor, 
       StartFromHydrogenic, StartFromPrevious, StrongField, StrongField2, Subshell, ShellSelection,
       tabulate, TestFrames, Triangle, tools,
       UseCoulomb, UseBabushkin, UseGauge,
       WeightedCartesian, W3j, W6j, W9j,
       Yields, Ylm, 
       Zeeman
     
# Basic data and data structures
include("module-Basics.jl");            using ..Basics
include("module-Radial.jl");            using ..Radial
include("module-Math.jl");              using ..Math
include("module-Defaults.jl");          using ..Defaults
include("module-ManyElectron.jl");      using ..ManyElectron
include("module-Nuclear.jl");           using ..Nuclear
include("module-BiOrthogonal.jl");      using ..BiOrthogonal

# Specialized functions/methods to manipulate these data
include("module-AngularMomentum.jl")
## include("module-AngularCoefficients-Ratip2013.jl")  ## keep for internal test purposes only
include("module-SpinAngular.jl");       using ..SpinAngular
include("module-Bsplines.jl");          using ..Bsplines
include("module-Continuum.jl")
include("module-Details.jl")
include("module-RadialIntegrals.jl")
include("module-HydrogenicIon.jl")
include("module-InteractionStrength.jl")
include("module-InteractionStrengthQED.jl")
include("module-PeriodicTable.jl")
include("module-TableStrings.jl")
include("module-Tools.jl")
include("module-AtomicState.jl");       using ..AtomicState
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
include("module-FormFactor.jl")
include("module-ReducedDensityMatrix.jl")
include("module-PlasmaShift.jl")
include("module-AlphaVariation.jl")
include("module-DecayYield.jl")
include("module-RadiativeOpacity.jl")
##x include("module-CloseCoupling.jl")
include("module-MultipolePolarizibility.jl")

# Functions/methods for atomic processes
include("module-PhotoExcitation.jl")
include("module-PhotoIonization.jl")
include("module-PhotoRecombination.jl")
include("module-AutoIonization.jl")
include("module-ElectronCapture.jl")
include("module-Dielectronic.jl")
include("module-PhotoExcitationFluores.jl")
include("module-PhotoExcitationAutoion.jl")
include("module-RayleighCompton.jl")
include("module-MultiPhotonDeExcitation.jl")
include("module-DoubleAutoIonization.jl")
include("module-PhotoDoubleIonization.jl")
include("module-CoulombExcitation.jl")
include("module-PhotoIonizationFluores.jl")
include("module-PhotoIonizationAutoion.jl")
include("module-CoulombIonization.jl")
include("module-ImpactExcitation.jl")
include("module-ImpactExcitationAutoion.jl")
include("module-RadiativeAuger.jl")
include("module-MultiPhotonIonization.jl")
include("module-MultiPhotonDoubleIon.jl")
include("module-InternalConversion.jl") 
#= Further processes, not yet included into the code
include("module-REDA.jl")
include("module-READI.jl")
include("module-PairProduction.jl")
include("module-PairAnnihilation1Photon.jl")
include("module-PairAnnihilation2Photon.jl")  =#

# Functions/methods for atomic responses and time evolutions
include("module-Pulse.jl")
# include("module-Statistical.jl")

# Functions/methods for the computation of atomic responses
## include("module-HighHarmonic.jl")
include("module-StrongField.jl") 

# Functions/methods for semi-empirical estimations
include("module-ImpactIonization.jl")
include("module-Semiempirical.jl")

# Functions/methods for atomic and cascade computations
include("module-Atomic.jl");            using ..Atomic
include("module-Plasma.jl");            using ..Plasma
include("module-Empirical.jl");         using ..Empirical
include("module-Cascade.jl");           using ..Cascade

# Functions/methods for symbolic computations
include("module-RacahAlgebra.jl");      using ..RacahAlgebra
include("module-SphericalTensor.jl");   using ..SphericalTensor

# Basic functions/methods to manipulate these data
include("module-BasicsAZ.jl")
include("module-ManyElectronAZ.jl")

# Specialized macros

# All test functions/methods stay with the JAC root module
include("module-TestFrames.jl");        using ..TestFrames
    
function __init__()
    # The following variables need to be initialized at runtime to enable precompilation
    global JAC_SUMMARY_IOSTREAM = stdout
    global JAC_TEST_IOSTREAM    = stdout
end

println("\nWelcome to JAC:  A community approach to the computation of atomic structures, cascades and time evolutions " * 
        "[(C) Copyright by Stephan Fritzsche, Jena (2018-2023)].")

end


