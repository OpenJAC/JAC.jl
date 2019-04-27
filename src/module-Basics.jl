
"""
`module JAC.Basics`  
    ... a submodel of JAC that contains various basic methods for dealing with general data structure.
"""
module Basics

    using  Dates, Printf, Interact,  LinearAlgebra 
    using  JAC, ..BasicTypes, ..Constants,  ..Radial, ..ManyElectron, ..Nuclear   ## , ..RadialIntegrals, ..InteractionStrength, ..AngularMomentum
    include("jac-add.jl")
    include("jac-analyze.jl")
    # include("jac-compute.jl")
    ##x include("jac-convert.jl")
    ##x include("jac-define.jl")
    include("jac-determine.jl")
    include("jac-diagonalize.jl")
    include("jac-display.jl")
    include("jac-estimate.jl")
    include("jac-exclude.jl")
    include("jac-evaluate.jl")
    include("jac-generate.jl")
    ##x include("jac-give.jl")
    include("jac-integrate.jl")
    include("jac-interpolate.jl")
    include("jac-merge.jl")
    include("jac-modify.jl")
    # include("jac-perform.jl")
    include("jac-provide.jl")
    # include("jac-plot.jl")
    include("jac-read.jl")
    # include("jac-recast.jl")
    # include("jac-sort.jl")
    # include("jac-tabulate.jl")
    ## include("jac-warn.jl")
    include("jac-tools.jl")

end # module


