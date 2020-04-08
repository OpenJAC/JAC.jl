
"""
`module  JAC.MultiPhotonDeExcitation`  
    ... a submodel of JAC that contains all methods for computing multi-photon excitation and decay rates.
"""
module MultiPhotonDeExcitation

    using Printf, QuadGK, ..AngularMomentum, ..AtomicState, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, 
                          ..PhotoEmission, ..Radial, ..TableStrings

    
    """
    `abstract type MultiPhotonDeExcitation.AbstractMultiPhotonProperty` 
        ... defines an abstract type to distinguish different multi-photon properties; however, by far not all of these properties
            are usually relevant for all the process::MultiPhotonDeExcitation.AbstractMultiPhotonProcess; cf. the JAC User Guide:
        
        + struct EnergyDiffCs           
                    ... Energy differential cross section for initially unpolarized atoms (TwoPhotonEmission, ...)
        + struct TotalCsLinear          
                    ... Total cross section for linerarly-polarized incident radiation (TwoPhotonAbsorptionMonochromatic, ...)
        + struct TotalCsRightCircular   
                    ... Total cross section for right-circularly polarized incident radiation (TwoPhotonAbsorptionMonochromatic, ...)
        + struct TotalCsUnpolarized     
                    ... Total cross section for unpolarized incident radiation (TwoPhotonAbsorptionMonochromatic, ...)
        + struct TotalCsDensityMatrix   
                    ... Total cross section for incident radiation with user-defined density matrix (not yet).
                        (TwoPhotonAbsorptionMonochromatic, ...)
    """
    abstract type  AbstractMultiPhotonProperty       end


    """
    `struct  MultiPhotonDeExcitation.EnergyDiffCs  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty`  
        ... a struct to request for the calculation of energy-diffrential cross sections in a two-photon emission process.
    """
    struct   EnergyDiffCs  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty      end

    function Base.show(property::EnergyDiffCs)
        sa = "energy-differential cross sections for multi-photon emission";  print(io, sa, "\n")
    end


    """
    `struct  MultiPhotonDeExcitation.TotalCsLinear  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty`  
        ... a struct to request for the calculation of the total (absorption) cross sections for linearly-polarized  
            incident radiation.
    """
    struct   TotalCsLinear  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty      end

    function Base.show(property::TotalCsLinear)
        sa = "total (absorption) cross sections for linearly-polarized radiation";  print(io, sa, "\n")
    end


    """
    `struct  MultiPhotonDeExcitation.TotalCsRightCircular  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty`  
        ... a struct to request for the calculation of the total (absorption) cross sections for right-circularly polarized  
            incident radiation.
    """
    struct   TotalCsRightCircular  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty      end

    function Base.show(property::TotalCsRightCircular)
        sa = "total (absorption) cross sections for right-circularly polarized radiation";  print(io, sa, "\n")
    end


    """
    `struct  MultiPhotonDeExcitation.TotalCsUnpolarized  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty`  
        ... a struct to request for the calculation of the total (absorption) cross sections for unpolarized  
            incident radiation.
    """
    struct   TotalCsUnpolarized  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProperty      end

    function Base.show(property::TotalCsUnpolarized)
        sa = "total (absorption) cross sections for unpolarized radiation";  print(io, sa, "\n")
    end
    

    
    """
    `abstract type MultiPhotonDeExcitation.AbstractMultiPhotonProcess` 
        ... defines an abstract type to distinguish different multi-photon processes; see also:
        
        + struct TwoPhotonAbsorptionMonochromatic
                    ... to specify a two-photon absorption process with monochromatic and equally-polarized photons,
                        usually from the same beam.
        + struct TwoPhotonAbsorptionBichromatic
                    ... to specify a two-photon absorption process by photons with well-defined frequency, propagation
                        direction and helicity, usually from different (oriented) beams.
        + struct TwoPhotonEmission
                    ... to specify a two-photon emission process, often for initially-unpolarized atoms. 
    """
    abstract type  AbstractMultiPhotonProcess       end


    """
    `struct  MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProcess`  
        ... a struct to specify a two-photon absorption process with monochromatic and equally-polarized photons,
            usually from the same beam.

        + properties        ::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1} 
            ... List of the multi-photon processes that are to be calculated for this absorption process.
    """
    struct   TwoPhotonAbsorptionMonochromatic  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProcess
        properties          ::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1} 
    end


    """
    `MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic()`  
        ... constructor for an 'default' instance of a MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic() process.
    """
    function TwoPhotonAbsorptionMonochromatic()
        TwoPhotonAbsorptionMonochromatic( [TotalCsLinear(), TotalCsUnpolarized()] )
    end


    function Base.string(process::MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic)
        sa = "Two-photon absorption process with mono-chromatic photons (from the same beam):"
        return( sa )
    end


    function Base.show(io::IO, process::MultiPhotonDeExcitation.TwoPhotonAbsorptionMonochromatic)
        sa = Base.string(process);               print(io, sa, "\n")
        println(io, "properties:                 $(process.properties)  ")
    end


    """
    `struct  MultiPhotonDeExcitation.TwoPhotonAbsorptionBichromatic  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProcess`  
        ... a struct to specify a two-photon absorption process by photons with well-defined frequency, propagation direction 
            and helicity, usually from different (oriented) beams.

        + properties        ::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1} 
            ... List of the multi-photon processes that are to be calculated for this absorption process.
        + omegaLess         ::Float64
            ... energy of the photon with smaller frequency; the larger frequency is derived from energy conservation.
    """
    struct   TwoPhotonAbsorptionBichromatic  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProcess
        properties          ::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1} 
        omegaLess           ::Float64
    end


    """
    `MultiPhotonDeExcitation.TwoPhotonAbsorptionBichromatic()`  
        ... constructor for an 'default' instance of a MultiPhotonDeExcitation.TwoPhotonAbsorptionBichromatic() process.
    """
    function TwoPhotonAbsorptionBichromatic()
        TwoPhotonAbsorptionBichromatic( [TotalCsUnpolarized()], 0.05 )
    end


    function Base.string(process::MultiPhotonDeExcitation.TwoPhotonAbsorptionBichromatic)
        sa = "Two-photon absorption process with bi-chromatic photons (from different beams):"
        return( sa )
    end


    function Base.show(io::IO, process::MultiPhotonDeExcitation.TwoPhotonAbsorptionBichromatic)
        sa = Base.string(process);               print(io, sa, "\n")
        println(io, "properties:                 $(process.properties)  ")
        println(io, "omegaLess:                  $(process.omegaLess)  ")
    end


    """
    `struct  MultiPhotonDeExcitation.TwoPhotonEmission  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProcess`  
        ... a struct to specify a two-photon emission process.

        + properties        ::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1} 
            ... List of the multi-photon processes that are to be calculated for this emission process.
        + noSharings  ::Int64   ... number of energy sharings for atomic transition, using the zeros of a Gauß-Legendre integration
    """
    struct   TwoPhotonEmission  <:  MultiPhotonDeExcitation.AbstractMultiPhotonProcess
        properties          ::Array{MultiPhotonDeExcitation.AbstractMultiPhotonProperty,1} 
        noSharings          ::Int64
    end


    """
    `MultiPhotonDeExcitation.TwoPhotonEmission()`  
        ... constructor for an 'default' instance of a MultiPhotonDeExcitation.TwoPhotonEmission() process.
    """
    function TwoPhotonEmission()
        TwoPhotonEmission( [EnergyDiffCs()], 4 )
    end


    function Base.string(process::MultiPhotonDeExcitation.TwoPhotonEmission)
        sa = "Two-photon emission process:"
        return( sa )
    end


    function Base.show(io::IO, process::MultiPhotonDeExcitation.TwoPhotonEmission)
        sa = Base.string(process);               print(io, sa, "\n")
        println(io, "properties:                 $(process.properties)  ")
    end

    
    
    """
    `struct  MultiPhotonDeExcitation.Settings`  ... defines a type for the settings in estimating multi-photon excitation and decay rates.

        + process                 ::MultiPhotonDeExcitation.AbstractMultiPhotonProcess      ... considered multi-photon process.
        + multipoles              ::Array{EmMultipole}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge}              ... Specifies the gauges to be included into the computations.
        + greenChannels           ::Array{AtomicState.GreenChannel,1}    ... Green function channels for these multi-photon calculations 
                                                                             which need to be generated independently.
        + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).

    """
    struct Settings 
        process                   ::MultiPhotonDeExcitation.AbstractMultiPhotonProcess
        multipoles                ::Array{EmMultipole}
        gauges                    ::Array{UseGauge}
        greenChannels             ::Array{AtomicState.GreenChannel,1}
        printBefore               ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1} 
    end 


    """
    `MultiPhotonDeExcitation.Settings()`  ... constructor for the default values of multi-photon excitation and decay rates.
    """
    function Settings()
        Settings( TwoPhotonEmission(), EmMultipole[], UseGauge[], AtomicState.GreenExpansion(), false, false, Tuple{Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::MultiPhotonDeExcitation.Settings)`  
    #		... prepares a proper printout of the variable settings::MultiPhotonDeExcitation.Settings.
     function Base.show(io::IO, settings::MultiPhotonDeExcitation.Settings) 
        println(io, "process:                  $(settings.process)  ")
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "gauges:                   $(settings.gauges)  ")
        print(io,   "greenChannels:             (settings.greenChannels)  ... with symmetries:   ")
        for channel in settings.greenChannels  print(io, "   $(channel.symmetry)")  end;   println(io, " ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
    end
    
    include("module-MultiPhotonDeExc-inc-2pAbsMono.jl")
    include("module-MultiPhotonDeExc-inc-2pEmission.jl")

end # module
