
"""
`module  JAC.BasicsAZ`  
 	   ... a submodel of JAC that contains methods from the Basics module but with reference to many other modules.
"""
module BasicsAZ

using  Printf,  LinearAlgebra, GaussQuadrature, JAC, ..AlphaVariation, ..AngularMomentum, ..Atomic, ..AtomicState, ..AutoIonization,
                ..Basics, ..Bsplines, ..Cascade, ..Continuum, ..DecayYield, ..Defaults, ..Dielectronic, 
                ..Einstein, ..FormFactor, ..Hfs, ..HydrogenicIon, ..InteractionStrength, ..InteractionStrengthQED, 
                ..IsotopeShift, ..LandeZeeman, ..LSjj, ..ManyElectron, ..MultipolePolarizibility, ..Nuclear, 
                ..PhotoEmission, ..PhotoExcitation, ..PhotoIonization, ..PlasmaShift, 
                ..Radial, ..RadialIntegrals, ..TableStrings


include("module-BascisAZ-inc-AG.jl")
include("module-BascisAZ-inc-compute.jl")
include("module-BascisAZ-inc-generate.jl")
include("module-BascisAZ-inc-HP.jl")
include("module-BascisAZ-inc-perform.jl")
include("module-BascisAZ-inc-QZ.jl")

end # module
