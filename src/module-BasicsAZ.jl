
"""
`module  JAC.BasicsAZ`  
 	   ... a submodel of JAC that contains methods from the Basics module but with reference to many other modules.
"""
module BasicsAZ

using  Printf,  LinearAlgebra, GaussQuadrature, JAC, 
       ..AngularMomentum, ..Atomic, ..AtomicState, ..Basics, ..Bsplines, ..HydrogenicIon, ..InteractionStrength, ..InteractionStrengthQED, 
       ..ManyElectron, ..Nuclear, ..PhotoEmission, ..Radial, ..RadialIntegrals, ..TableStrings 
       
#==    ..AlphaVariation,  ..AutoIonization, ..Cascade, ..Continuum, ..DecayYield, ..Defaults, ..Dielectronic, ..Einstein, 
       ..FormFactor, ..Hfs, ..IsotopeShift, ..LandeZeeman, ..LSjj, ..MultipolePolarizibility,  ..PhotoExcitation, 
       ..PhotoIonization  ==#


include("module-BascisAZ-inc-AG.jl")
include("module-BascisAZ-inc-compute.jl")
include("module-BascisAZ-inc-generate.jl")
include("module-BascisAZ-inc-HP.jl")
include("module-BascisAZ-inc-perform.jl")
include("module-BascisAZ-inc-QZ.jl")

end # module
