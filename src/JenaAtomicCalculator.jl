#
#  Use git:                    git status ....;   git add <filenames>;   git commit -m "..";   git push;   git rm  <filenames>
#  Use Jupyter notebooks:      using IJulia;   notebook()
#  Activation:                 ];   pkg> up;   pkg> activate
#  Working with JAC:           using Revise;   using JAC;   include("../src/jac.jl");   pkg> test
#  
#  Copy to desktop             scp -r JAC.jl/ fritzsch@10.140.119.236:~/fri/.
"""
`module JenaAtomicCalculator`  
    ... Jena Atomic Calculator (JAC) provides tools for performing atomic (structure) calculations at various degrees of complexity 
        and sophistication. It has been designed to not only calculate atomic level structures and properties [such as g-factors or
        hyperfine and isotope-shift parameters] but also transition amplitudes between bound-state levels [for the dipole 
        operator, etc.] and, in particular, (atomic) transition probabilities, Auger rates, photoionization cross sections, 
        radiative and dielectronic recombination rates as well as cross sections for several other -- elementary or composed --
        processes. 

"""
module JenaAtomicCalculator

const JAC = JenaAtomicCalculator


# Restrict the size and functionality of code by just including certain modules, while others are not taken into account.
# Of course, this "exclusion" may cause later errors if missing functionality is invoked.  The selective use of code
# is introduced mainly for development purposes and for keeping the storage requirement moderate (if one wishes to focus
# on certain classes of applications)
incProperties           = true
incBasicProcesses       = true
incAdvancedProcesses    = true
incCascades             = true  ## Requires: incBasicProcesses
incPlasma               = true  ## Requires: incProperties
incStrongField          = true
incAtomicCompass        = true
incRacahAlgebra         = true

#==
incAdvancedProcesses    = false
incCascades             = false  ## Requires: incBasicProcesses
incPlasma               = false  ## Requires: incProperties
incStrongField          = false
incAtomicCompass        = false
incRacahAlgebra         = false  ==#



using  Dates,  Printf,  BSplineKit, LinearAlgebra, SpecialFunctions, QuadGK, Cubature, GSL, JLD2, SymEngine, 
       HypergeometricFunctions  ## , Interact, GaussQuadrature, IJulia, FortranFiles

export AbstractCImethod, AbstractConfigurationRestriction, AbstractEeInteraction, AbstractPotential, AbstractQedModel, AbstractStartOrbitals,
       AbstractProcessSettings, AbstractEmpiricalSettings, AbstractPlasmaModel, AbstractPropertySettings, AbstractLineShiftSettings,
       add, AlphaX, AlphaVariation, analyze, AnapoleMoment, 
       AngularJ64, AngularM64, AngularJ, AngularMomentum, 
       AsfSettings, Atomic, AtomicState, AtomicStructure, Auger, AugerInPlasma, AutoIonization, AverageAtom, AtomicCompass,
       Basics, Basis, Beam, BeamPhotoExcitation, BreitInteraction, Bsplines, BsplinesN,
       CartesianVector, Cartesian2DFieldVector, Cartesian3DFieldVector, CiSettings, CiExpansion, ClebschGordan, CloseCoupling, 
       compute, convertUnits, Compton, Configuration, ConfigurationR, 
       Cascade, Continuum, CorePolarization, Coulex, CoulombExcitation, Coulion, CoulombBreit, CoulombGaunt, 
       CoulombInteraction, CoulombIonization, CsfR,        
       diagonalize, Defaults, DecayYield, DielectronicRecombination, Dierec, Djpq, DoubleAutoIonization, DoubleAuger, 
       DiagonalCoulomb, DefaultQuantizationAxis, 
       Eimex, ElectronCapture, ElecCapture, estimate, ElectricDipoleMoment, Einstein, EinsteinX, EmMultipole, evaluate, ExpStokes, 
       Empirical,
       E1, M1, E2, M2, E3, M3, E4, M4,
       FormFactor, FormF, FullCIeigen,
       generate, GreenSettings, GreenChannel, GreenExpansion, getDefaults, Green, Gui,
       Hamiltonian, Hfs, HyperfineInduced, HighHarmonic, HFS, HydrogenicIon, HarmonicQuantizationAxis,
       interpolate, integrate, Integral, ImpactExcAuto, ImpactExcitation, ImpactExcitationAutoion, ImpactIonization, 
       InteractionStrength, InternalConv, InternalConversion, InternalRecombination, Isotope, IsotopeShift, IsotopicFraction, 
       Kronecker,
       LandeF, LandeJ, LandeZeeman, Level, LevelSelection, LevelSymmetry, LineSelection, LSjj, LSjjSettings, LeftCircular,
       ManyElectron, MeanFieldSettings, MeanFieldBasis, MeanFieldMultiplet, minus, Model, modify, 
       MultiPhotonDE, MultiPhotonDeExcitation, MultiPhotonDoubleIon, 
       MultiPI, MultiPDI, MultiPhotonIonization, MultipoleMoment, MultipolePolarizibility, Multiplet, 
       NoAmplitude, Nuclear, NoneQed, NoProcess, NoProperty, 
       OneElectronSettings, OneElectronSpectrum, Orbital, oplus,
       PairA1P, PairAnnihilation1Photon, PairAnnihilation2Photon, PairProduction, Parity, ParityNonConservation, ParticleScattering, 
       PathwaySelection, PeriodicTable, perform, 
       Photo, PhotoDouble, PhotoDoubleIonization, PhotoEmission, PhotoExc, PhotoExcAuto, PhotoExcFluor, 
       PhotoExcitation, PhotoExcitationAutoion, PhotoExcitationFluores, PhotoIonAuto, PhotoIonFluor, PhotoIonization, 
       PhotoIonizationAutoion, PhotoIonizationFluores, PhotoRecombination,   
       Plasma, plus, Polarity, provide, Pulse,
       QedPetersburg, QedSydney, 
       RacahAlgebra, RacahExpression, Radial, RadialIntegrals, Radiative, RadiativeAuger, RAuger, RasSettings, RasStep, 
       RasExpansion, RayleighCompton, recast, Rec, REDA, READI, Representation, ReducedDensityMatrix, RadiativeOpacity,
       RestrictMaximumDisplacements, RestrictNoElectronsTo, RestrictParity, RestrictToShellDoubles, RequestMinimumOccupation, RequestMaximumOccupation,
       ResonantInelastic,
       SchiffMoment, Semiempirical, setDefaults, Shell, ShellSelection, SolidAngle, Spectroscopy, SphericalTensor, SpinAngular, StarkShift,
       StartFromHydrogenic, StartFromPrevious, StrongField, StrongField2, Subshell, StaticQuantizationAxis, StaticField, SelfConsistent,
       tabulate, TestFrames, tools, Triangle, TwoElectronOnePhoton, TimeHarmonicField,
       UseBabushkin, UseCoulomb, UseGauge,
       WeightedCartesian, W3j, W6j, W9j,
       Yields, Ylm, 
       Zeeman
     
# Basic data and data structures
include("module-Basics.jl");            using ..Basics
include("module-Radial.jl");            using ..Radial
include("module-Math.jl");              using ..Math
include("module-Defaults.jl");          using ..Defaults
include("module-Distribution.jl");      using ..Defaults
include("module-ManyElectron.jl");      using ..ManyElectron
include("module-Nuclear.jl");           using ..Nuclear
# include("module-BiOrthogonal.jl");      using ..BiOrthogonal


# Specialized functions/methods to manipulate these data
include("module-AngularMomentum.jl")
## include("module-AngularCoefficients-Ratip2013.jl")  ## keep for internal test purposes only
include("module-SpinAngular.jl");       using ..SpinAngular
##x include("module-Bsplines.jl");          using ..Bsplines
include("module-BsplinesN.jl");         using ..BsplinesN
include("module-Pulse.jl")
include("module-Beam.jl")
include("module-Continuum.jl")
##x include("module-Details.jl")
include("module-RadialIntegrals.jl");   using ..RadialIntegrals
include("module-HydrogenicIon.jl")
include("module-InteractionStrength.jl")
include("module-InteractionStrengthQED.jl")
include("module-Hamiltonian.jl");       using ..Hamiltonian
include("module-SelfConsistent.jl");    using ..SelfConsistent
include("module-PeriodicTable.jl")
include("module-TableStrings.jl")
include("module-Tools.jl")
include("module-AtomicState.jl");       using ..AtomicState
include("module-LSjj.jl");              using ..LSjj

include("module-PhotoEmission.jl")

if  incProperties
# Functions/methods for atomic amplitudes
include("module-MultipoleMoment.jl")
include("module-ParityNonConservation.jl")
# Functions/methods for atomic properties
include("module-Einstein.jl")    
include("module-Hfs.jl")
include("module-IsotopeShift.jl")
include("module-LandeZeeman.jl")
include("module-StarkShift.jl")
include("module-FormFactor.jl")
include("module-ReducedDensityMatrix.jl")
include("module-AlphaVariation.jl")
include("module-RadiativeOpacity.jl")
include("module-MultipolePolarizibility.jl")
end

if  incBasicProcesses
# Functions/methods for atomic processes
include("module-PhotoExcitation.jl")
include("module-PhotoIonization.jl")
include("module-PhotoRecombination.jl")
include("module-AutoIonization.jl")
include("module-ElectronCapture.jl")
include("module-DielectronicRecombination.jl")
include("module-PhotoExcitationFluores.jl")
include("module-PhotoExcitationAutoion.jl")
include("module-RayleighCompton.jl")
include("module-ParticleScattering.jl")
include("module-BeamPhotoExcitation.jl") 
include("module-HyperfineInduced.jl") 
include("module-ResonantInelastic.jl") 
include("module-DecayYield.jl")
include("module-ImpactExcitation.jl")
end

if incAdvancedProcesses
# Functions/methods for more advanced atomic processes
include("module-MultiPhotonDeExcitation.jl")
include("module-CoulombExcitation.jl")
include("module-CoulombIonization.jl")
include("module-PhotoDoubleIonization.jl")
include("module-PhotoIonizationFluores.jl")
include("module-PhotoIonizationAutoion.jl")
include("module-ImpactExcitationAutoion.jl")
include("module-RadiativeAuger.jl")
include("module-MultiPhotonIonization.jl")
include("module-MultiPhotonDoubleIon.jl")
include("module-InternalConversion.jl") 
include("module-InternalRecombination.jl") 
include("module-TwoElectronOnePhoton.jl") 
include("module-DoubleAutoIonization.jl")
#= Further processes, not yet included into the code
include("module-REDA.jl")
include("module-READI.jl")
include("module-PairProduction.jl")
include("module-PairAnnihilation1Photon.jl")
include("module-PairAnnihilation2Photon.jl")  =#
end

# Functions/methods for atomic responses and time evolutions
# include("module-Statistical.jl")

if incStrongField
# Functions/methods for the computation of atomic responses
## include("module-HighHarmonic.jl")
include("module-StrongField.jl") 
end

if incAtomicCompass
# Functions/methods for the computation of atomic-compass simulations
include("module-AtomicCompass.jl") 
end

# Functions/methods for semi-empirical estimations
include("module-ImpactIonization.jl")
include("module-Semiempirical.jl")
include("module-Empirical.jl");         using ..Empirical

# Functions/methods for atomic computations
include("module-Atomic.jl");            using ..Atomic

if  incPlasma
# Functions/methods for plasma computations
include("module-Plasma.jl");            using ..Plasma
end


if  incCascades
# Functions/methods for cascade computations
include("module-Cascade.jl");           using ..Cascade
end

# Functions/methods for symbolic computations
if  incRacahAlgebra
include("module-RacahAlgebra.jl");      using ..RacahAlgebra
include("module-SphericalTensor.jl");   using ..SphericalTensor
end

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

println("\nWelcome to JenaAtomicCalculator (JAC):  A community approach to the computation of atomic structures, " *
        "cascades and time evolutions [(C) Copyright by Stephan Fritzsche, Jena (2018-2025)].")
        

end


