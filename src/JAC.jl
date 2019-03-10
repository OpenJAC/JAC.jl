#
#  Use git:                    git status ....;   git add <filenames>;   git commit -m "..";   git push;   git rm  <filenames>
#  Use Jupyter notebooks:      using IJulia;   notebook()
#  Activation:                 ];   pkg> up;   pkg> activate
#  Working with JAC:           using Revise;   using JAC;   include("../src/jac.jl");   pkg> test
#
"""
`module JAC`  ... Jena Atomic Calculator (JAC) provides tools for performing atomic (structure) calculations at various degrees of complexity 
                  and sophistication. It has been designed to not only calculate atomic level structures and properties [such as g-factors or
                  hyperfine and isotope-shift parameters] but also transition amplitudes between bound-state levels [for the anapole moment, dipole 
                  operator, electron electric-dipole moment, parity non-conservation, etc.] and, in particular, (atomic) transition probabilities, 
                  Auger rates, photoionization cross sections, radiative and dielectronic recombination rates as well as cross sections for many 
                  other (elementary) processes. In the future, JAC will also facilitate interactive computations, the simulation of atomic cascades, 
                  the time-evolution of statistical tensors as well as various semi-empirical estimates of atomic properties. -- 
                  In addition, the JAC module supports the display of level energies, electron and photon spectra, radial orbitals and 
                  and other atomic data.


**`Perform (atomic) computations of different complexity:`**  JAC will eventually support **seven kinds** of computations which can be 
                    summarized as follows:

   + Atomic computations, based on explicitly specified electron configurations.
   + Restricted active-space computations (RAS; not yet properly implemented).
   + Interactive computations.
   + Atomic cascade computations (not yet fully implemented).
   + Atomic responses (not yet implemented).
   + Time-evolution of statistical tensors in (intense) light pusles (not yet implemented).
   + Semi-empirical estimates of cross sections, etc. (not yet properly implemented).


**`Further details and information`**

      + Kinds of atomic implementation                                             [cf. ? JAC.kindsOfComputation]
      + Atomic amplitudes (partly) implemented in JAC                              [cf. ? JAC.amplitudes]
      + Atomic level properties (partly) implemented in JAC                        [cf. ? JAC.properties]
      + Atomic processes (partly) implemented in JAC                               [cf. ? JAC.processes]
      + Interactive use of JAC procedures                                          [cf. ? JAC.interactive]
      + Design principles and limitations of the JAC program                       [cf. ? JAC.design]
      + Data types, structs and name conventions of the JAC module                 [cf. ? JAC.datatypes]
      + Atomic cascade computations and approximations                             [cf. ? JAC.decayCascades]
      + Use of (em) light pulses in the time evolution of statistical tensors      [cf. ? JAC.pulses]
      + Why Julia ?                                                                [cf. ? JAC.whyJulia]

"""
module JAC

using  Dates, Printf,  LinearAlgebra, Interact, SpecialFunctions, FortranFiles, GaussQuadrature, QuadGK, GSL   #, JLD

export AlphaVariation, AnapoleMoment, AsfSettings, Atomic, Auger, Basis, Cascade, Configuration, ConfigurationR, Continuum, CsfR, Radial,
       CoulombExcitation, CoulombIonization,  
       DecayYield, Dielectronic, DoubleAuger,
       ElectricDipoleMoment, Einstein, 
       FormFactor,
       GreenFunction,
       Hfs, HydrogenicIon,
       ImpactExcitation, ImpactExcitationAutoion, ImpactIonization, InternalConversion, IsotopeShift, LandeZeeman, Level, 
       MultiPhotonDeExcitation, MultiPhotonDoubleIon, MultiPhotonIonization, MultipoleMoment,  MultipolePolarizibility, 
       Multiplet, 
       NoAmplitude, NoProcess, Nuclear, Model,
       Orbital, PairAnnihilation1Photon, PairAnnihilation2Photon, PairProduction, ParityNonConservation,
       PhotoExcitation, PhotoExcitationAutoion, PhotoExcitationFluores, PhotoIonization, PhotoIonizationFluores, 
       PhotoIonizationAutoion, PhotoRecombination, PlasmaShift, 
       Radiative, RadiativeAuger, RayleighCompton, REDA, READI,
       SchiffMoment, Shell, Subshell
    
global JAC_counter = 0

include("jac-basics.jl")
include("module-Radial.jl")
include("module-ManyElectron.jl")
include("module-Nuclear.jl")
include("module-Bsplines.jl")
include("module-Continuum.jl")

include("module-Einstein.jl")
include("module-Hfs.jl")
include("module-IsotopeShift.jl")
include("module-LandeZeeman.jl")
include("module-AlphaVariation.jl")
include("module-FormFactor.jl")
include("module-DecayYield.jl")
include("module-GreenFunction.jl")
include("module-MultipolePolarizibility.jl")
include("module-PlasmaShift.jl")

include("module-MultipoleMoment.jl")
include("module-ParityNonConservation.jl")

include("module-Radiative.jl")
include("module-PhotoExcitation.jl")
include("module-PhotoIonization.jl")
include("module-PhotoRecombination.jl")
include("module-Auger.jl")
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
include("module-InternalConversion.jl")
#=
include("module-DoubleAuger.jl")
include("module-REDA.jl")
include("module-READI.jl")
include("module-PairProduction.jl")
include("module-PairAnnihilation1Photon.jl")
include("module-PairAnnihilation2Photon.jl")  =#

include("module-Atomic.jl")
include("module-Cascade.jl")
include("module-HydrogenicIon.jl")
# include("module-Pulse.jl")
# include("module-Statistical.jl")

include("module-AngularMomentum.jl")
include("module-AngularCoefficients-Ratip2013.jl")
include("module-InteractionStrength.jl")
include("module-InteractionStrengthQED.jl")
include("module-Math.jl")
include("module-MessageHandling.jl")
include("module-RadialIntegrals.jl")
include("module-TableStrings.jl")
include("module-Tools.jl")
##x include("module-Tutorials.jl")

# include("module-ImpactIonization.jl")
# include("module-Semiempirical.jl")

using  JAC.Radial, JAC.ManyElectron, JAC.Nuclear, JAC.RadialIntegrals, JAC.InteractionStrength, JAC.AngularMomentum

include("jac-add.jl")
include("jac-analyze.jl")
include("jac-compute.jl")
include("jac-convert.jl")
include("jac-define.jl")
include("jac-determine.jl")
include("jac-diagonalize.jl")
include("jac-display.jl")
include("jac-estimate.jl")
include("jac-exclude.jl")
include("jac-evaluate.jl")
include("jac-generate.jl")
include("jac-give.jl")
include("jac-integrate.jl")
include("jac-interpolate.jl")
include("jac-merge.jl")
include("jac-modify.jl")
include("jac-perform.jl")
include("jac-provide.jl")
# include("jac-plot.jl")
include("jac-read.jl")
include("jac-recast.jl")
include("jac-sort.jl")
include("jac-store.jl")
include("jac-tabulate.jl")
include("jac-warn.jl")

include("jac-constants.jl")
include("jac-document.jl")

include("jac-test.jl")
include("jac-tools.jl")

println("\nWelcome to JAC:  A fresh computational approach to atomic structures, processes, casacdes and time evolutions [(C) Copyright by Stephan Fritzsche, Jena (2019)].")

end


