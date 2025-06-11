
"""
`module  JAC.BasicsAZ`  
 	   ... a submodel of JAC that contains methods from the Basics module but with reference to many other modules.
"""
module BasicsAZ

using  Printf,  LinearAlgebra, GaussQuadrature, JenaAtomicCalculator, ..Basics, ..TableStrings
       ## using JenaAtomicCalculator ... since otherwise almost all other modules must be included explicitly
       
#==    ..AlphaVariation, ..AngularMomentum, ..Atomic, ..AtomicState,  ..AutoIonization,  ..Basics, ..BsplinesN, ..Cascade, 
       ..DecayYield, ..Defaults, ..DielectronicRecombination, ..DoubleAutoIonization, ..Einstein, ..FormFactor,
       ..HydrogenicIon, ..Hfs, ..HyperfineInduced, ..ImpactExcitation, ..ImpactExcitationAutoion, ..InteractionStrength,  
       ..InteractionStrengthQED,  ..InternalRecombination, ..IsotopeShift, 
       ..LandeZeeman, ..LSjj, ..ManyElectron,  ..MultiPhotonDeExcitation, ..MultipolePolarizibility,  ..Nuclear, 
       ..PhotoDoubleIonization, ..PhotoEmission,  ..PhotoExcitation,  ..PhotoExcitationAutoion, ..PhotoExcitationFluores, 
       ..PhotoIonization, ..PhotoIonizationFluores, ..PhotoIonizationAutoion, ..PhotoRecombination,
       ..Radial, ..RadialIntegrals, ..RadiativeAuger, ..RayleighCompton, ..ResonantInelastic, ..SelfConsistent, ..SpinAngular, ..StarkShift, 
       ..TableStrings, ..TwoElectronOnePhoton
       
       ..AlphaVariation,  ..AutoIonization, ..Cascade, ..Continuum, ..DecayYield, ..Defaults, ..DielectronicRecombination, ..Einstein, 
       ..FormFactor, ..Hfs, ..IsotopeShift, ..LandeZeeman, ..LSjj, ..MultipolePolarizibility,  ..PhotoExcitation, 
       ..PhotoIonization  ==#


include("module-BascisAZ-inc-AG.jl")
include("module-BascisAZ-inc-compute.jl")
include("module-BascisAZ-inc-generate.jl")
include("module-BascisAZ-inc-HP.jl")
include("module-BascisAZ-inc-perform.jl")
include("module-BascisAZ-inc-QZ.jl")

end # module
