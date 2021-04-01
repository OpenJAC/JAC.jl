
"""
`module  JAC.Cascade
    ... a submodel of JAC that contains all methods to set-up and process (simple) atomic cascade computations.
"""
module Cascade

    using Dates, JLD2, Printf, ..AutoIonization, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..DecayYield, ..Dielectronic,
                               ..ElectronCapture, ..Radial, ..ManyElectron, ..Nuclear, 
                               ..PhotoEmission, ..PhotoExcitation, ..PhotoIonization, ..Semiempirical, ..TableStrings

    
    """
    `abstract type Cascade.AbstractCascadeScheme` 
        ... defines an abstract type to distinguish different excitation, ionization and decay schemes of an atomic cascade; see also:
        
        + struct StepwiseDecayScheme       
            ... to represent a standard decay scheme, starting from the levels of one or several initial multiplets.
        + struct PhotonIonizationScheme    
            ... to represent a (prior) photoionization part of a cascade for a list of given photon energies.
        + struct PhotonExcitationScheme    
            ... to represent a (prior) photo-excitation part of a cascade for a specified set of excitations.
        + struct DielectronicCaptureScheme  
            ... to represent the (dielectronic) capture of an electron for a maximum excitation energy and a given list of
                subshells; the excitation energy refers to the levels of the reference configurations, and cascade blocks are built
                only with the given subshells.
        + struct HollowIonScheme    
            ... to represent the capture of one or several electrons into a list of subshells; various distributions of the
                electron among these shells are supported. For the subsequent decay, the list of decay shells need to be specified
                as well.
        + struct ElectronExcitationScheme  
            ... to represent a (prior) electron-impact excitation part of a cascade for given electron energies (not yet).
    """
    abstract type  AbstractCascadeScheme       end


    """
    `struct  Cascade.StepwiseDecayScheme  <:  Cascade.AbstractCascadeScheme`  
        ... a struct to represent (and generate) a mean-field orbital basis.

        + processes             ::Array{Basics.AbstractProcess,1} 
            ... List of the atomic processes that are supported and should be included into the cascade.
        + maxElectronLoss       ::Int64             
            ... (Maximum) Number of electrons in which the initial- and final-state configurations can differ from each other; 
                this also determines the maximal steps of any particular decay path.
        + chargeStateShifts     ::Dict{Int64,Float64} 
            ... (N => en) total energy shifts of all levels with N electrons; these shifts [in a.u.] help open/close decay 
                channels by simply shifting the total energies of all levels.
        + NoShakeDisplacements  ::Int64             
            ... Maximum number of electron displacements due to shake-up  or shake-down processes in any individual step of cascade.
        + shakeFromShells       ::Array{Shell,1}        ... List of shells from which shake transitions may occur.
        + shakeToShells         ::Array{Shell,1}        ... List of shells into which shake transitions may occur.
    """
    struct   StepwiseDecayScheme  <:  Cascade.AbstractCascadeScheme
        processes               ::Array{Basics.AbstractProcess,1}
        maxElectronLoss         ::Int64
        chargeStateShifts       ::Dict{Int64,Float64}
        NoShakeDisplacements    ::Int64
        shakeFromShells         ::Array{Shell,1}
        shakeToShells           ::Array{Shell,1}
    end


    """
    `Cascade.StepwiseDecayScheme()`  ... constructor for an 'default' instance of a Cascade.StepwiseDecayScheme.
    """
    function StepwiseDecayScheme()
        StepwiseDecayScheme([Radiative()], 0, Dict{Int64,Float64}(), 0, Shell[], Shell[] )
    end


    # `Base.string(scheme::StepwiseDecayScheme)`  ... provides a String notation for the variable scheme::StepwiseDecayScheme.
    function Base.string(scheme::StepwiseDecayScheme)
        sa = "Stepwise decay (scheme) of an atomic cascade with:"
        return( sa )
    end


    # `Base.show(io::IO, scheme::StepwiseDecayScheme)`  ... prepares a proper printout of the scheme::StepwiseDecayScheme.
    function Base.show(io::IO, scheme::StepwiseDecayScheme)
        sa = Base.string(scheme);                print(io, sa, "\n")
        println(io, "processes:                  $(scheme.processes)  ")
        println(io, "maxElectronLoss:            $(scheme.maxElectronLoss)  ")
        println(io, "chargeStateShifts:          $(scheme.chargeStateShifts)  ")
        println(io, "NoShakeDisplacements:       $(scheme.NoShakeDisplacements)  ")
        println(io, "shakeFromShells:            $(scheme.shakeFromShells)  ")
        println(io, "shakeToShells:              $(scheme.shakeToShells)  ")
    end


    """
    `struct  Cascade.PhotonIonizationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... a struct to represent (and generate) a mean-field orbital basis.

        + processes             ::Array{Basics.AbstractProcess,1} 
            ... List of the atomic processes that are supported and should be included into the cascade.
        + maxPhotoElectrons     ::Int64                 
            ... (Maximum) Number of photo-electrons in which the initial-and (photo-ionized) final-state configurations can differ 
                from each other; this affects the number of ionized multiplets in the cascade.
        + photonEnergies        ::Array{Float64,1}      
            ... List of photon energies for which this photo-ionization scheme is to be calculated.
    """
    struct   PhotonIonizationScheme  <:  Cascade.AbstractCascadeScheme
        processes               ::Array{Basics.AbstractProcess,1}
        maxPhotoElectrons       ::Int64
        photonEnergies          ::Array{Float64,1}
    end


    """
    `Cascade.PhotonIonizationScheme()`  ... constructor for an 'default' instance of a Cascade.PhotonIonizationScheme.
    """
    function PhotonIonizationScheme()
        PhotonIonizationScheme([PhotoIonization], 0, Float64[] )
    end


    # `Base.string(scheme::PhotonIonizationScheme)`  ... provides a String notation for the variable scheme::PhotonIonizationScheme.
    function Base.string(scheme::PhotonIonizationScheme)
        sa = "Photo-ionization (scheme) as prior part of an atomic cascade:"
        return( sa )
    end


    # `Base.show(io::IO, scheme::PhotonIonizationScheme)`  ... prepares a proper printout of the scheme::PhotonIonizationScheme.
    function Base.show(io::IO, scheme::PhotonIonizationScheme)
        sa = Base.string(scheme);                print(io, sa, "\n")
        println(io, "processes:                  $(scheme.processes)  ")
        println(io, "maxPhotoElectrons:          $(scheme.maxPhotoElectrons)  ")
        println(io, "photonEnergies:             $(scheme.photonEnergies)  ")
    end


    """
    `struct  Cascade.PhotonExcitationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... a struct to define and describe a photo-excitation calculation for an atom in some initial state/configuration
            and for a given set of shell-excitations

        + processes             ::Array{Basics.AbstractProcess,1} 
            ... List of the atomic processes that are supported and should be included into the cascade.  
        + multipoles            ::Array{EmMultipole}           
            ... Multipoles of the radiation field that are to be included into the excitation processes.
        + minPhotonEnergy       ::Float64                 
            ... Minimum photon energy [in a.u.] that restrict the number of excited configurations to be taken into accout.
        + maxPhotonEnergy       ::Float64                 
            ... Maximum photon energy [in a.u.] that restrict the number of excited configurations to be taken into accout.
        + NoExcitations         ::Int64                 
            ... (Maximum) Number of electron replacements with regard to the initial configurations/multiplets.
        + excitationFromShells  ::Array{Shell,1}    
            ... List of shells from which photo-excitations are to be considered.
        + excitationToShells    ::Array{Shell,1}    
            ... List of shells into which photo-excitations are to be considered, including possibly already occupied shells.
    """
    struct   PhotonExcitationScheme  <:  Cascade.AbstractCascadeScheme
        processes               ::Array{Basics.AbstractProcess,1}    
        multipoles              ::Array{EmMultipole}  
        minPhotonEnergy         ::Float64                 
        maxPhotonEnergy         ::Float64                 
        NoExcitations           ::Int64
        excitationFromShells    ::Array{Shell,1}
        excitationToShells      ::Array{Shell,1}
    end


    """
    `Cascade.PhotonExcitationScheme()`  ... constructor for an 'default' instance of a Cascade.PhotonExcitationScheme.
    """
    function PhotonExcitationScheme()
        PhotonExcitationScheme([PhotoExc()], [E1], 1.0, 100., 0, Shell[], Shell[] )
    end


    # `Base.string(scheme::PhotonExcitationScheme)`  ... provides a String notation for the variable scheme::PhotoAbsorptionScheme.
    function Base.string(scheme::PhotonExcitationScheme)
        sa = "Photon-excitation (scheme) due to processes:"
        return( sa )
    end


    # `Base.show(io::IO, scheme::PhotonExcitationScheme)`  ... prepares a proper printout of the scheme::PhotonExcitationScheme.
    function Base.show(io::IO, scheme::PhotonExcitationScheme)
        sa = Base.string(scheme);                print(io, sa, "\n")
        println(io, "processes:                  $(scheme.processes)  ")
        println(io, "multipoles:                 $(scheme.multipoles)  ")
        println(io, "minPhotonEnergy:            $(scheme.minPhotonEnergy)  ")
        println(io, "maxPhotonEnergy:            $(scheme.maxPhotonEnergy)  ")
        println(io, "NoExcitations :             $(scheme.NoExcitations )  ")
        println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
        println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
    end


    """
    `struct  Cascade.DielectronicCaptureScheme  <:  Cascade.AbstractCascadeScheme`  
        ... a struct to define and describe the (dielectronic) capture of electrons for an atom in some initial state/configuration;
            for such a scheme, the doubly-excited configurations due to the electron capture are generated automatically due to
            given maximal numbers of the (into-) shells (nl) as well as the maximum displacement with regard to the initial configuration.
            An additional maximum excitation energy need to be provided due to the maximum temperatures for which DR plasma rate coefficients
            are to be determined, cf. Basics.convert().

        + multipoles            ::Array{EmMultipole}           
            ... Multipoles of the radiation field that are to be included into the radiative stabilization processes.
        + maxExcitationEnergy   ::Float64                 
            ... Maximum excitation energy [in a.u.] with regard to the reference configurations/levels that restrict the number 
                of excited configurations to be taken into accout. This maximum excitation energy has to be derived from the maximum 
                temperature for which DR coefficients need to be derived and is typically set to 5x T_e,max.
        + minPhotonEnergy        ::Float64                 
            ... Minimum (mean) photon energy [in a.u.] in the radiative stabilization of the doubly-excited configurations; 
                If cascade blocks are separated be less than this energy, the radiative stabilization is neglected.
        + NoExcitations         ::Int64                 
            ... (Maximum) Number of electron replacements in the doubly-excited configuration with regard to the initial 
                configurations/multiplets, apart from one additional electron due to the electron capture itself.
        + excitationFromShells  ::Array{Shell,1}    
            ... List of shells from which excitations are to be considered.
        + excitationToShells  ::Array{Shell,1}    
            ... List of shells to which (core-shell) excitations are to be considered.
        + intoShells            ::Array{Shell,1}
            ... List of shells into which electrons are initially placed (captured).
        + decayShells           ::Array{Shell,1}
            ... List of shells into which electrons the electrons can decay (apart from the core shells).
    """
    struct   DielectronicCaptureScheme  <:  Cascade.AbstractCascadeScheme
        multipoles              ::Array{EmMultipole}  
        maxExcitationEnergy     ::Float64   
        minPhotonEnergy         ::Float64                 
        NoExcitations           ::Int64
        excitationFromShells    ::Array{Shell,1}
        excitationToShells      ::Array{Shell,1}
        intoShells              ::Array{Shell,1}
        decayShells             ::Array{Shell,1}
    end


    """
    `Cascade.DielectronicCaptureScheme()`  ... constructor for an 'default' instance of a Cascade.DielectronicCaptureScheme.
    """
    function DielectronicCaptureScheme()
        DielectronicCaptureScheme([E1], 1.0, 0., 1, Shell[], Shell[], Shell[], Shell[] )
    end


    # `Base.string(scheme::DielectronicCaptureScheme)`  ... provides a String notation for the variable scheme::DielectronicCaptureScheme.
    function Base.string(scheme::DielectronicCaptureScheme)
        sa = "Dielectronic capture & stabilization (scheme):"
        return( sa )
    end


    # `Base.show(io::IO, scheme::DielectronicCaptureScheme)`  ... prepares a proper printout of the scheme::DielectronicCaptureScheme.
    function Base.show(io::IO, scheme::DielectronicCaptureScheme)
        sa = Base.string(scheme);                print(io, sa, "\n")
        println(io, "multipoles:                 $(scheme.multipoles)  ")
        println(io, "maxExcitationEnergy:        $(scheme.maxExcitationEnergy)  ")
        println(io, "minPhotonEnergy:            $(scheme.minPhotonEnergy)  ")
        println(io, "NoExcitations:              $(scheme.NoExcitations)  ")
        println(io, "excitationFromShells:       $(scheme.excitationFromShells)  ")
        println(io, "excitationToShells:         $(scheme.excitationToShells)  ")
        println(io, "intoShells:                 $(scheme.intoShells)  ")
        println(io, "decayShells:                $(scheme.decayShells)  ")
    end


    """
    `struct  Cascade.HollowIonScheme  <:  Cascade.AbstractCascadeScheme`  
        ... a struct to define and describe the formation and decay of a hollow ion, e.g. an electronic core configuration 
            into which one or several additional electrons are captured into (high) nl shells. Both the shell (lists) for 
            the initial capture (intoShells) and the subsequent decay (decayShells) need to be specified explicitly to 
            readily control the size of the computations.

        + processes             ::Array{Basics.AbstractProcess,1} 
            ... List of the atomic processes that are supported and should be included into the decay scheme.  
        + multipoles            ::Array{EmMultipole,1}           
            ... Multipoles of the radiation field that are to be included into the radiative stabilization processes.
        + NoCapturedElectrons   ::Int64   
            ... Number of captured electrons, e.g. placed in the intoShells.
        + intoShells            ::Array{Shell,1}
            ... List of shells into which electrons are initially placed (captured).
        + decayShells           ::Array{Shell,1}
            ... List of shells into which electrons the electrons can decay (apart from the core shells).
    """
    struct   HollowIonScheme  <:  Cascade.AbstractCascadeScheme
        processes               ::Array{Basics.AbstractProcess,1}    
        multipoles              ::Array{EmMultipole}  
        NoCapturedElectrons     ::Int64
        intoShells              ::Array{Shell,1}
        decayShells             ::Array{Shell,1}
    end


    """
    `Cascade.HollowIonScheme()`  ... constructor for an 'default' instance of a Cascade.HollowIonScheme.
    """
    function HollowIonScheme()
        HollowIonScheme([Auger(), Radiative()], [E1], 0, Shell[], Shell[] )
    end


    # `Base.string(scheme::HollowIonScheme)`  ... provides a String notation for the variable scheme::HollowIonScheme.
    function Base.string(scheme::HollowIonScheme)
        sa = "Hollow ion (scheme) due to processes:"
        return( sa )
    end


    # `Base.show(io::IO, scheme::HollowIonScheme)`  ... prepares a proper printout of the scheme::HollowIonScheme.
    function Base.show(io::IO, scheme::HollowIonScheme)
        sa = Base.string(scheme);                print(io, sa, "\n")
        println(io, "processes:                  $(scheme.processes)  ")
        println(io, "multipoles:                 $(scheme.multipoles)  ")
        println(io, "NoCapturedElectrons:        $(scheme.NoCapturedElectrons)  ")
        println(io, "intoShells:                 $(scheme.intoShells)  ")
        println(io, "decayShells:                $(scheme.decayShells)  ")
    end


    """
    `struct  Cascade.ElectronExcitationScheme  <:  Cascade.AbstractCascadeScheme`  
        ... a struct to represent (and generate) a mean-field orbital basis.

        + processes             ::Array{Basics.AbstractProcess,1} 
            ... List of the atomic processes that are supported and should be included into the cascade.
        + electronEnergies      ::Array{Float64,1}                
            ... List of electron energies for which this electron-impact excitation scheme is to be calculated.
    """
    struct   ElectronExcitationScheme  <:  Cascade.AbstractCascadeScheme
        processes               ::Array{Basics.AbstractProcess,1}
        electronEnergies        ::Array{Float64,1}
    end


    """
    `Cascade.ElectronExcitationScheme()`  ... constructor for an 'default' instance of a Cascade.ElectronExcitationScheme.
    """
    function ElectronExcitationScheme()
        ElectronExcitationScheme([Eimex], Float64[] )
    end


    # `Base.string(scheme::ElectronExcitationScheme)`  ... provides a String notation for the variable scheme::ElectronExcitationScheme.
    function Base.string(scheme::ElectronExcitationScheme)
        sa = "Electron-impact excitation (scheme) as prior part of an atomic cascade:"
        return( sa )
    end


    # `Base.show(io::IO, scheme::ElectronExcitationScheme)`  ... prepares a proper printout of the scheme::ElectronExcitationScheme.
    function Base.show(io::IO, scheme::ElectronExcitationScheme)
        sa = Base.string(scheme);                 print(io, sa, "\n")
        println(io, "processes:                   $(scheme.processes)  ")
        println(io, "electronEnergies:            $(scheme.electronEnergies)  ")
    end
    

    """
    `abstract type Cascade.AbstractCascadeApproach` 
        ... defines an abstract and a number of singleton types for the computational approach/model that is applied in order to 
            generate and evaluate all many-electron amplitudes of a given cascade.

      + struct AverageSCA         
        ... all levels in the cascade are described in single-configuration and single-CSF approximation; this (rather crude) approach 
            neglects all configuration-interactions and also applies just a single set of one-electron orbitals (from the least-ionized charge
            state) for all considered charge states.
            
      + struct SCA                
        ... all levels in the cascade are described in single-configuration approximation but with 'mixtures' within the configuration;
            an individual mean-field is generated for each charge state and all continuum orbitals are generated for the correct transition
            energy in the field of the remaining ion. Moreover, all the fine-structure transitions are calculated individually.
    """
    abstract type  AbstractCascadeApproach                   end
    struct         AverageSCA  <:  AbstractCascadeApproach   end
    struct         SCA         <:  AbstractCascadeApproach   end
    struct         UserMCA     <:  AbstractCascadeApproach   end
    


    """
    `struct  Cascade.Block`  
        ... defines a type for an individual block of configurations that are treatet together within the cascade. Such an block is given 
            by a list of configurations that may occur as initial- and/or final-state configurations in some step of the canscade and that 
            give rise to a common multiplet in order to allow for configuration interactions but to avoid 'double counting' of individual 
            levels in the cascade.

        + NoElectrons     ::Int64                     ... Number of electrons in this block.
        + confs           ::Array{Configuration,1}    ... List of one or several configurations that define the multiplet.
        + hasMultiplet    ::Bool                      
            ... true if the (level representation in the) multiplet has already been computed and false otherwise.
        + multiplet       ::Multiplet                 ... Multiplet of the this block.
    """
    struct  Block
        NoElectrons       ::Int64
        confs             ::Array{Configuration,1} 
        hasMultiplet      ::Bool
        multiplet         ::Multiplet  
    end 


    """
    `Cascade.Block()`  ... constructor for an 'empty' instance of a Cascade.Block.
    """
    function Block()
        Block( )
    end


    # `Base.show(io::IO, block::Cascade.Block)`  ... prepares a proper printout of the variable block::Cascade.Block.
    function Base.show(io::IO, block::Cascade.Block) 
        println(io, "NoElectrons:        $(block.NoElectrons)  ")
        println(io, "confs :             $(block.confs )  ")
        println(io, "hasMultiplet:       $(block.hasMultiplet)  ")
        println(io, "multiplet:          $(block.multiplet)  ")
    end


    """
    `struct  Cascade.Step`  
        ... defines a type for an individual step of an excitation and/or decay cascade. Such a step is determined by the two lists 
            of initial- and final-state configuration as well as by the atomic process, such as Auger, PhotoEmission, or others, which related
            the initial- and final-state levels to each other. Since the (lists of) initial- and final-state configurations treated (each)
            by a single multiplet (for parities and total angular momenta), a cascade step supports full configuration interaction within
            the multiplet but also help  avoid 'double counting' of individual levels. Indeed, each electron configuration may occur only in
            one cascade block. In contrast, each list of initial- and final-state multiplets (cascade blocks) can occur in quite different 
            steps due to the considered processes and parallel decay pathes in a cascade.

        + process          ::JBasics.AbstractProcess   ... Atomic process that 'acts' in this step of the cascade.
        + settings         ::Union{PhotoEmission.Settings, AutoIonization.Settings, PhotoIonization.Settings, PhotoExcitation.Settings,
                                   DielectronicCapture.Settings}        
                                                       ... Settings for this step of the cascade.
        + initialConfigs   ::Array{Configuration,1}    ... List of one or several configurations that define the initial-state multiplet.
        + finalConfigs     ::Array{Configuration,1}    ... List of one or several configurations that define the final-state multiplet.
        + initialMultiplet ::Multiplet                 ... Multiplet of the initial-state levels of this step of the cascade.
        + finalMultiplet   ::Multiplet                 ... Multiplet of the final-state levels of this step of the cascade.
    """
    struct  Step
        process            ::Basics.AbstractProcess
        settings           ::Union{PhotoEmission.Settings, AutoIonization.Settings, PhotoIonization.Settings, PhotoExcitation.Settings}
        initialConfigs     ::Array{Configuration,1}
        finalConfigs       ::Array{Configuration,1}
        initialMultiplet   ::Multiplet
        finalMultiplet     ::Multiplet
    end 


    """
    `Cascade.Step()`  ... constructor for an 'empty' instance of a Cascade.Step.
    """
    function Step()
        Step( Basics.NoProcess, PhotoEmission.Settings, Configuration[], Configuration[], Multiplet(), Multiplet())
    end


    # `Base.show(io::IO, step::Cascade.Step)`  ... prepares a proper printout of the variable step::Cascade.Step.
    function Base.show(io::IO, step::Cascade.Step) 
        println(io, "process:                $(step.process)  ")
        println(io, "settings:               $(step.settings)  ")
        println(io, "initialConfigs:         $(step.initialConfigs)  ")
        println(io, "finalConfigs:           $(step.finalConfigs)  ")
        println(io, "initialMultiplet :      $(step.initialMultiplet )  ")
        println(io, "finalMultiplet:         $(step.finalMultiplet)  ")
    end


    """
    `struct  Cascade.Computation`  
        ... defines a type for a cascade computation, i.e. for the computation of a whole photon excitation, photon ionization and/or 
            decay cascade. The -- input and control -- data from this computation can be modified, adapted and refined to the practical needs, 
            and before the actual computations are carried out explictly. Initially, this struct just contains the physical meta-data about the 
            cascade, that is to be calculated, but a new instance of the same Cascade.Computation gets later enlarged in course of the 
            computation in order to keep also wave functions, level multiplets, etc.

        + name               ::String                          ... A name for the cascade
        + nuclearModel       ::Nuclear.Model                   ... Model, charge and parameters of the nucleus.
        + grid               ::Radial.Grid                     ... The radial grid to be used for the computation.
        + asfSettings        ::AsfSettings                     ... Provides the settings for the SCF process.
        + scheme             ::Cascade.AbstractCascadeScheme   ... Scheme of the atomic cascade (photoionization, decay, ...)
        + approach           ::Cascade.AbstractCascadeApproach 
            ... Computational approach/model that is applied to generate and evaluate the cascade; possible approaches are: 
                {AverageSCA(), SCA(), ...}
        + initialConfs       ::Array{Configuration,1}          
            ... List of one or several configurations that contain the level(s) from which the cascade starts.
        + initialMultiplets  ::Array{Multiplet,1}              
            ... List of one or several (initial) multiplets; either initialConfs 'xor' initialMultiplets can  be specified 
                for a given cascade computation.
    """
    struct  Computation
        name                 ::String
        nuclearModel         ::Nuclear.Model
        grid                 ::Radial.Grid
        asfSettings          ::AsfSettings
        scheme               ::Cascade.AbstractCascadeScheme 
        approach             ::Cascade.AbstractCascadeApproach
        initialConfigs       ::Array{Configuration,1}
        initialMultiplets    ::Array{Multiplet,1} 
    end 


    """
    `Cascade.Computation()`  ... constructor for an 'default' instance of a Cascade.Computation.
    """
    function Computation()
        Computation("Default cascade computation",  Nuclear.Model(10.), Radial.Grid(), AsfSettings(), Cascade.StepwiseDecayScheme(), 
                    Cascade.AverageSCA(), Configuration[], Multiplet[] )
    end


    """
    `Cascade.Computation(comp::Cascade.Computation;`
        
                name=..,               nuclearModel=..,             grid=..,              asfSettings=..,     
                scheme=..,             approach=..,                 initialConfigs=..,    initialMultiplets=..)
                
        ... constructor for re-defining the computation::Cascade.Computation.
    """
    function Computation(comp::Cascade.Computation;                              
        name::Union{Nothing,String}=nothing,                                  nuclearModel::Union{Nothing,Nuclear.Model}=nothing,
        grid::Union{Nothing,Radial.Grid}=nothing,                             asfSettings::Union{Nothing,AsfSettings}=nothing,  
        scheme::Union{Nothing,Cascade.AbstractCascadeScheme}=nothing,         approach::Union{Nothing,Cascade.AbstractCascadeApproach}=nothing, 
        initialConfigs::Union{Nothing,Array{Configuration,1}}=nothing,        initialMultiplets::Union{Nothing,Array{Multiplet,1}}=nothing)
        
        if  initialConfigs != nothing   &&   initialMultiplets != nothing   
            error("Only initialConfigs=..  'xor'  initialMultiplets=..  can be specified for a Cascade.Computation().")
        end

        if  name                 == nothing   namex                 = comp.name                    else  namex = name                                   end 
        if  nuclearModel         == nothing   nuclearModelx         = comp.nuclearModel            else  nuclearModelx = nuclearModel                   end 
        if  grid                 == nothing   gridx                 = comp.grid                    else  gridx = grid                                   end 
        if  asfSettings          == nothing   asfSettingsx          = comp.asfSettings             else  asfSettingsx = asfSettings                     end 
        if  scheme               == nothing   schemex               = comp.scheme                  else  schemex = scheme                               end 
        if  approach             == nothing   approachx             = comp.approach                else  approachx = approach                           end 
        if  initialConfigs       == nothing   initialConfigsx       = comp.initialConfigs          else  initialConfigsx = initialConfigs               end 
        if  initialMultiplets    == nothing   initialMultipletsx    = comp.initialMultiplets       else  initialMultipletsx = initialMultiplets         end
    	
    	Computation(namex, nuclearModelx, gridx, asfSettingsx, schemex, approachx, initialConfigsx, initialMultipletsx)
    end


    # `Base.string(comp::Cascade.Computation)`  ... provides a String notation for the variable comp::Cascade.Computation.
    function Base.string(comp::Cascade.Computation)
        if        typeof(comp.scheme) == Cascade.StepwiseDecayScheme       sb = "stepwise decay scheme"
        elseif    typeof(comp.scheme) == Cascade.PhotonIonizationScheme    sb = "(prior) photo-ionization scheme"
        elseif    typeof(comp.scheme) == Cascade.PhotonExcitationScheme    sb = "photo-excitation scheme"
        elseif    typeof(comp.scheme) == Cascade.DielectronicCaptureScheme     sb = "(di-) electronic capture scheme"
        elseif    typeof(comp.scheme) == Cascade.HollowIonScheme           sb = "hollow ion scheme"
        else      error("unknown typeof(comp.scheme)")
        end
        
        sa = "Cascade computation   $(comp.name)  for a $sb,  in $(comp.approach) approach as well as "
        sa = sa * "for Z = $(comp.nuclearModel.Z) and initial configurations: \n "
        if  comp.initialConfigs    != Configuration[]  for  config  in  comp.initialConfigs     sa = sa * string(config) * ",  "     end     end
        if  comp.initialMultiplets != Multiplet[]      
            for  mp      in  comp.initialMultiplets  sa = sa * mp.name * " with " * string(length(mp.levels)) * " levels,  "         end     end
        return( sa )
    end


    # `Base.show(io::IO, comp::Cascade.Computation)`  ... prepares a proper printout comp::Cascade.Computation.
    function Base.show(io::IO, comp::Cascade.Computation)
        sa = Base.string(comp)
        sa = sa * "\n ... in addition, the following parameters/settings are defined: ";       print(io, sa, "\n")
        println(io, "> cascade scheme:           $(comp.scheme)  ")
        println(io, "> nuclearModel:             $(comp.nuclearModel)  ")
        println(io, "> grid:                     $(comp.grid)  ")
        println(io, "> asfSettings:            \n$(comp.asfSettings) ")
    end
    
    #######################################################################################################################################
    #######################################################################################################################################
    #######################################################################################################################################


    """
    abstract type  Cascade.AbstractData`  
        ... defines an abstract type to distinguish different output data of a cascade computation; see also:

        + struct DecayData       ... to comprise the amplitudes/rates from a stepswise decay cascade.
        + struct PhotoIonData    ... to comprise the amplitudes/rates from a photo-ionization part of a cascade.
        + struct ExcitationData  ... to comprise the amplitudes/rates from a photo-excitation part of a cascade.
        + struct CaptureData     ... to comprise the amplitudes/rates from a electron-capture part of a cascade.
    """
    abstract type  AbstractData      end


    """
    `struct  Cascade.DecayData  <:  Cascade.AbstractData`  ... defines a type for an atomic cascade, i.e. lists of radiative and Auger lines.

        + linesR         ::Array{PhotoEmission.Line,1}          ... List of radiative lines.
        + linesA         ::Array{AutoIonization.Line,1}         ... List of Auger lines.
    """  
    struct  DecayData  <:  Cascade.AbstractData
        linesR           ::Array{PhotoEmission.Line,1}
        linesA           ::Array{AutoIonization.Line,1}
    end 


    """
    `Cascade.DecayData()`  ... (simple) constructor for cascade  decay data.
    """
    function DecayData()
        DecayData(Array{PhotoEmission.Line,1}[], Array{AutoIonization.Line,1}[] )
    end


    # `Base.show(io::IO, data::Cascade.DecayData)`  ... prepares a proper printout of the variable data::Cascade.DecayData.
    function Base.show(io::IO, data::Cascade.DecayData) 
        println(io, "linesR:                  $(data.linesR)  ")
        println(io, "linesA:                  $(data.linesA)  ")
    end


    """
    `struct  Cascade.PhotoIonData  <:  Cascade.AbstractData`  ... defines a type for an atomic cascade, i.e. lists of photoionization lines.

        + photonEnergy   ::Float64                              ... Photon energy for this (part of the) cascade.
        + linesP         ::Array{PhotoIonization.Line,1}        ... List of photoionization lines.
    """  
    struct  PhotoIonData  <:  Cascade.AbstractData
        photonEnergy     ::Float64
        linesP           ::Array{PhotoIonization.Line,1}
    end 


    """
    `Cascade.PhotoIonData()`  ... (simple) constructor for cascade photo-ionization data.
    """
    function PhotoIonData()
        PhotoIonData(0., Array{PhotoIonization.Line,1}[] )
    end


    # `Base.show(io::IO, data::Cascade.PhotoIonData)`  ... prepares a proper printout of the variable data::Cascade.PhotoIonData.
    function Base.show(io::IO, data::Cascade.PhotoIonData) 
        println(io, "photonEnergy:            $(data.photonEnergy)  ")
        println(io, "linesP:                  $(data.linesP)  ")
    end


    """
    `struct  Cascade.ExcitationData  <:  Cascade.AbstractData`  ... defines a type for an atomic cascade, i.e. lists of excitation lines.

        + linesE         ::Array{PhotoExcitation.Line,1}        ... List of photoionization lines.
    """  
    struct  ExcitationData  <:  Cascade.AbstractData
        linesE           ::Array{PhotoExcitation.Line,1}
    end 


    """
    `Cascade.ExcitationData()`  ... (simple) constructor for cascade excitation data.
    """
    function ExcitationData()
        ExcitationData(Array{PhotoExcitation.Line,1}[] )
    end


    # `Base.show(io::IO, data::Cascade.ExcitationData)`  ... prepares a proper printout of the variable data::Cascade.ExcitationData.
    function Base.show(io::IO, data::Cascade.ExcitationData) 
        println(io, "linesE:                  $(data.linesE)  ")
    end
    

    #######################################################################################################################################
    #######################################################################################################################################
    #######################################################################################################################################

    
    """
    `abstract type  Cascade.AbstractSimulationProperty`  
        ... defines an abstract and various singleton types for the different properties that can be obtained from the simulation of 
            cascade data.

        + struct IonDistribution         
            ... simulate the 'ion distribution' as it is found after all cascade processes are completed.
        + struct FinalLevelDistribution  
            ... simulate the 'final-level distribution' as it is found after all cascade processes are completed.
        + struct PhotoAbsorption  
            ... simulate the photoabsorption cross sections for a given set of photo-excitation and ionization processes.
        + struct DecayPathes             ... determine the major 'decay pathes' of the cascade.
        + struct ElectronIntensities     ... simulate the electron-line intensities as function of electron energy.
        + struct PhotonIntensities       ... simulate  the photon-line intensities as function of electron energy. 
        + struct DrRateCoefficients      ... simulate  the DR (plasma) rate coefficients for given plasma temperatures. 
        + struct ElectronCoincidence     ... simulate electron-coincidence spectra.
    """
    abstract type  AbstractSimulationProperty                              end


    """
    `struct  Cascade.IonDistribution   <:  Cascade.AbstractSimulationProperty`  
        ... defines a type for simulating the 'ion distribution' as it is found after all cascade processes are completed.

        + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
            ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    """  
    struct  IonDistribution   <:  Cascade.AbstractSimulationProperty
        initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    end 


    """
    `Cascade.IonDistribution()`  ... (simple) constructor for cascade IonDistribution data.
    """
    function IonDistribution()
        IonDistribution([(1, 1.0)])
    end


    # `Base.show(io::IO, data::Cascade.IonDistribution)`  ... prepares a proper printout of the variable data::Cascade.IonDistribution.
    function Base.show(io::IO, dist::Cascade.IonDistribution) 
        println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    end


    """
    `struct  Cascade.FinalLevelDistribution   <:  Cascade.AbstractSimulationProperty`  
        ... defines a type for simulating the 'final-level distribution' as it is found after all cascade processes are completed.

        + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
            ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    """  
    struct  FinalLevelDistribution   <:  Cascade.AbstractSimulationProperty
        initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    end 


    """
    `Cascade.FinalLevelDistribution()`  ... (simple) constructor for cascade FinalLevelDistribution.
    """
    function FinalLevelDistribution()
        FinalLevelDistribution([(1, 1.0)])
    end


    # `Base.show(io::IO, data::Cascade.FinalLevelDistribution)`  
    #       ... prepares a proper printout of the variable data::Cascade.FinalLevelDistribution.
    function Base.show(io::IO, dist::Cascade.FinalLevelDistribution) 
        println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    end


    """
    `struct  Cascade.ElectronIntensities   <:  Cascade.AbstractSimulationProperty`  
        ... defines a type for simulating the electron-line intensities as function of electron energy.

        + minElectronEnergy   ::Float64     ... Minimum electron energy for the simulation of electron spectra.
        + maxElectronEnergy   ::Float64     ... Maximum electron energy for the simulation of electron spectra.
        + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
            ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    """  
    struct  ElectronIntensities   <:  Cascade.AbstractSimulationProperty
        minElectronEnergy     ::Float64
        maxElectronEnergy     ::Float64
        initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    end 


    """
    `Cascade.ElectronIntensities()`  ... (simple) constructor for cascade ElectronIntensities.
    """
    function ElectronIntensities()
        ElectronIntensities(0., 1.0e6,  [(1, 1.0)])
    end


    # `Base.show(io::IO, data::Cascade.ElectronIntensities)`  ... prepares a proper printout of the variable data::Cascade.ElectronIntensities.
    function Base.show(io::IO, dist::Cascade.ElectronIntensities) 
        println(io, "minElectronEnergy:        $(dist.minElectronEnergy)  ")
        println(io, "maxElectronEnergy:        $(dist.maxElectronEnergy)  ")
        println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    end


    """
    `struct  Cascade.PhotonIntensities   <:  Cascade.AbstractSimulationProperty`  
        ... defines a type for simulating the photon-line intensities as function of photon energy.

        + minPhotonEnergy     ::Float64     ... Minimum photon energy for the simulation of photon spectra.
        + maxPhotonEnergy     ::Float64     ... Maximum photon energy for the simulation of photon spectra.
        + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
            ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
        + leadingConfigs      ::Array{Configuration,1}   
            ... List of leading configurations whose levels are equally populated, either initially or ....
    """  
    struct  PhotonIntensities   <:  Cascade.AbstractSimulationProperty
        minPhotonEnergy       ::Float64
        maxPhotonEnergy       ::Float64
        initialOccupations    ::Array{Tuple{Int64,Float64},1} 
        leadingConfigs        ::Array{Configuration,1}
    end 


    """
    `Cascade.PhotonIntensities()`  ... (simple) constructor for cascade PhotonIntensities.
    """
    function PhotonIntensities()
        PhotonIntensities(0., 1.0e6,  [(1, 1.0)], Configuration[])
    end


    # `Base.show(io::IO, dist::Cascade.PhotonIntensities)`  ... prepares a proper printout of the variable data::Cascade.PhotonIntensities.
    function Base.show(io::IO, dist::Cascade.PhotonIntensities) 
        println(io, "minPhotonEnergy:          $(dist.minPhotonEnergy)  ")
        println(io, "maxPhotonEnergy:          $(dist.maxPhotonEnergy)  ")
        println(io, "initialOccupations:       $(dist.initialOccupations)  ")
        println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
    end


    """
    `struct  Cascade.DrRateCoefficients   <:  Cascade.AbstractSimulationProperty`  
        ... defines a type for simulating the DR plasma rate coefficients as function of the (free) electron energy and
            plasma temperature.

        + minPhotonEnergy     ::Float64     ... Minimum photon energy to be included in the radiative stabilization.
        + maxPhotonEnergy     ::Float64     ... Maximum photon energy to be included in the radiative stabilization.
        + temperatures        ::Array{Float64,1}
            ... temperatures [K] for which the DR plasma rate coefficieints to be calculated.
        + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
            ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
        + leadingConfigs      ::Array{Configuration,1}   
            ... List of leading configurations whose levels are equally populated, either initially or ....
    """  
    struct  DrRateCoefficients   <:  Cascade.AbstractSimulationProperty
        minPhotonEnergy       ::Float64
        maxPhotonEnergy       ::Float64
        temperatures          ::Array{Float64,1}
        initialOccupations    ::Array{Tuple{Int64,Float64},1} 
        leadingConfigs        ::Array{Configuration,1}
    end 


    """
    `Cascade.DrRateCoefficients()`  ... (simple) constructor for cascade DrRateCoefficients.
    """
    function DrRateCoefficients()
        DrRateCoefficients(0., 1.0e6,  Float64[], [(1, 1.0)], Configuration[])
    end


    # `Base.show(io::IO, dist::Cascade.DrRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.DrRateCoefficients.
    function Base.show(io::IO, dist::Cascade.DrRateCoefficients) 
        println(io, "minPhotonEnergy:          $(dist.minPhotonEnergy)  ")
        println(io, "maxPhotonEnergy:          $(dist.maxPhotonEnergy)  ")
        println(io, "temperatures:             $(dist.temperatures)  ")
        println(io, "initialOccupations:       $(dist.initialOccupations)  ")
        println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
    end
    
    
    struct   PhotoAbsorption              <:  AbstractSimulationProperty   end
    struct   DecayPathes                  <:  AbstractSimulationProperty   end
    struct   ElectronCoincidence          <:  AbstractSimulationProperty   end


    """
    abstract type  Cascade.AbstractSimulationMethod`  
        ... defines a abstract and a list of singleton data types for the properties that can be 'simulated' from a given
            set of a set Cascade.AbstractData.

        + struct ProbPropagation     ... to propagate the (occupation) probabilites of the levels until no further changes occur.
        + struct MonteCarlo          ... to simulate the cascade decay by a Monte-Carlo approach of possible pathes (not yet considered).
        + struct RateEquations       ... to solve the cascade by a set of rate equations (not yet considered).
    """
    abstract type  AbstractSimulationMethod                  end
    struct   ProbPropagation  <:  AbstractSimulationMethod   end
    struct   MonteCarlo       <:  AbstractSimulationMethod   end
    struct   RateEquations    <:  AbstractSimulationMethod   end


    """
    `struct  Cascade.SimulationSettings`  ... defines settings for performing the simulation of some cascade (data).

        + printTree           ::Bool        ... Print the cascade tree in a short form
        + printLongTree       ::Bool        ... Print the cascade tree in a long form
        + initialPhotonEnergy ::Float64     ... Photon energy for which photoionization data are considered. 
    """
    struct  SimulationSettings
        printTree             ::Bool
        printLongTree         ::Bool 
        initialPhotonEnergy   ::Float64
    end 


    """
    `Cascade.SimulationSettings()`  ... constructor for an 'empty' instance of a Cascade.Block.
    """
    function SimulationSettings()
        SimulationSettings(false, false, 0.)
    end


    # `Base.show(io::IO, settings::SimulationSettings)`  ... prepares a proper printout of the variable settings::SimulationSettings.
    function Base.show(io::IO, settings::SimulationSettings) 
        println(io, "printTree:                $(settings.printTree)  ")
        println(io, "printLongTree:            $(settings.printLongTree)  ")
        println(io, "initialPhotonEnergy:      $(settings.initialPhotonEnergy)  ")
    end


    """
    `struct  Cascade.Simulation`  ... defines a simulation on some given cascade (data).

        + name            ::String                              ... Name of the simulation
        + property        ::Cascade.AbstractSimulationProperty 
            ... Property that is to be considered in this simulation of the cascade (data).
        + method          ::Cascade.AbstractSimulationMethod    
            ... Method that is used in the cascade simulation; cf. Cascade.SimulationMethod.
        + settings        ::Cascade.SimulationSettings          ... Settings for performing these simulations.
        + computationData ::Array{Dict{String,Any},1}           ... Date on which the simulations are performed
    """
    struct  Simulation
        name              ::String
        property          ::Cascade.AbstractSimulationProperty
        method            ::Cascade.AbstractSimulationMethod
        settings          ::Cascade.SimulationSettings 
        computationData   ::Array{Dict{String,Any},1}
    end 


    """
    `Cascade.Simulation()`  ... constructor for an 'default' instance of a Cascade.Simulation.
    """
    function Simulation()
        Simulation("Default cascade simulation", Cascade.PhotoAbsorption(), Cascade.ProbPropagation(), 
                   Cascade.SimulationSettings(), Array{Dict{String,Any},1}[] )
    end


    """
    `Cascade.Simulation(sim::Cascade.Simulation;`
        
                name=..,               property=..,             method=..,              settings=..,     computationData=.. )
                
        ... constructor for re-defining the computation::Cascade.Simulation.
    """
    function Simulation(sim::Cascade.Simulation;                              
        name::Union{Nothing,String}=nothing,                                  property::Union{Nothing,Cascade.AbstractSimulationProperty}=nothing,
        method::Union{Nothing,Cascade.AbstractSimulationMethod}=nothing,      settings::Union{Nothing,Cascade.SimulationSettings}=nothing,    
        computationData::Union{Nothing,Array{Dict{String,Any},1}}=nothing )
 
        if  name            == nothing   namex            = sim.name              else  namex            = name                end 
        if  property        == nothing   propertyx        = sim.property          else  propertyx        = property            end 
        if  method          == nothing   methodx          = sim.method            else  methodx          = method              end 
        if  settings        == nothing   settingsx        = sim.settings          else  settingsx        = settings            end 
        if  computationData == nothing   computationDatax = sim.computationData   else  computationDatax = computationData     end 
    	
    	Simulation(namex, propertyx, methodx, settingsx, computationDatax)
    end


    # `Base.show(io::IO, simulation::Cascade.Simulation)`  ... prepares a proper printout of the variable simulation::Cascade.Simulation.
    function Base.show(io::IO, simulation::Cascade.Simulation) 
        println(io, "name:              $(simulation.name)  ")
        println(io, "property:          $(simulation.property)  ")
        println(io, "method:            $(simulation.method)  ")
        println(io, "computationData:   ... based on $(length(simulation.computationData)) cascade data sets ")
        println(io, "> settings:      \n$(simulation.settings)  ")
    end


    """
    `struct  Cascade.LineIndex`  ... defines a line index with regard to the various lineLists of data::Cascade.LineIndex.

        + lineSet      ::Cascade.AbstractData    ... refers to the data set for which this index is defined.
        + process      ::Basics.AbstractProcess  ... refers to the particular lineList of cascade (data).
        + index        ::Int64                   ... index of the corresponding line.
    """
    struct  LineIndex
        lineSet        ::Cascade.AbstractData 
        process        ::Basics.AbstractProcess
        index          ::Int64 
    end 


    # `Base.show(io::IO, index::Cascade.LineIndex)`  ... prepares a proper printout of the variable index::Cascade.LineIndex.
    function Base.show(io::IO, index::Cascade.LineIndex) 
        println(io, "lineSet (typeof):      $(typeof(index.lineSet))  ")
        println(io, "process:               $(index.process)  ")
        println(io, "index:                 $(index.index)  ")
    end


    """
    `mutable struct  Cascade.Level`  ... defines a level specification for dealing with cascade transitions.

        + energy       ::Float64                     ... energy of the level.
        + J            ::AngularJ64                  ... total angular momentum of the level
        + parity       ::Basics.Parity               ... total parity of the level
        + NoElectrons  ::Int64                       ... total number of electrons of the ion to which this level belongs.
        + relativeOcc  ::Float64                     ... relative occupation  
        + parents      ::Array{Cascade.LineIndex,1}  ... list of parent lines that (may) populate the level.     
        + daugthers    ::Array{Cascade.LineIndex,1}  ... list of daugther lines that (may) de-populate the level.     
    """
    mutable struct  Level
        energy         ::Float64 
        J              ::AngularJ64 
        parity         ::Basics.Parity 
        NoElectrons    ::Int64 
        relativeOcc    ::Float64 
        parents        ::Array{Cascade.LineIndex,1} 
        daugthers      ::Array{Cascade.LineIndex,1} 
    end 


    # `Base.show(io::IO, level::Cascade.Level)`  ... prepares a proper printout of the variable level::Cascade.Level.
    function Base.show(io::IO, level::Cascade.Level) 
        println(io, "energy:        $(level.energy)  ")
        println(io, "J:             $(level.J)  ")
        println(io, "parity:        $(level.parity)  ")
        println(io, "NoElectrons:   $(level.NoElectrons)  ")
        println(io, "relativeOcc:   $(level.relativeOcc)  ")
        println(io, "parents:       $(level.parents)  ")
        println(io, "daugthers:     $(level.daugthers)  ")
    end
    
    
    """
    `Base.:(==)(leva::Cascade.Level, levb::Cascade.Level)`  ... returns true if both levels are equal and false otherwise.
    """
    function  Base.:(==)(leva::Cascade.Level, levb::Cascade.Level)
        if  leva.energy == levb.energy   &&   leva.J == levb.J  && leva.parity == levb.parity  && leva.NoElectrons == levb.NoElectrons 
                return( true )
        else    return( false )
        end
    end


    """
    `struct  Cascade.AbsorptionCrossSection`  
        ... defines the absorption cross section for a particular (incident) photon energy in terms of its discrete and 
            (direct photoionization) contributions. Of course, this absorption cross section depends on the relative population
            of the initial levels.

        + photonEnergy ::Float64              ... incident photon energy/photon-energy dependence of the absorption spectrum.
        + excitationCS ::Basics.EmProperty    ... contribution due to discrete excitation processes.
        + ionizationCS ::Basics.EmProperty    ... contribution due to contineous ionization processes.
    """
    struct  AbsorptionCrossSection
        photonEnergy   ::Float64
        excitationCS   ::Basics.EmProperty
        ionizationCS   ::Basics.EmProperty
    end 


    # `Base.show(io::IO, cs::Cascade.AbsorptionCrossSection)`  ... prepares a proper printout of the variable cs::Cascade.AbsorptionCrossSection.
    function Base.show(io::IO, cs::Cascade.AbsorptionCrossSection) 
        println(io, "photonEnergy:      $(cs.photonEnergy)  ")
        println(io, "excitationCS:      $(cs.excitationCS)  ")
        println(io, "ionizationCS:      $(cs.ionizationCS)  ")
    end
    
    
    include("module-Cascade-inc-computations.jl")
    include("module-Cascade-inc-dielectronic-capture.jl")
    include("module-Cascade-inc-hollow-ion.jl")
    include("module-Cascade-inc-photoabsorption.jl")
    include("module-Cascade-inc-photoexcitation.jl")
    include("module-Cascade-inc-photoionization.jl")
    include("module-Cascade-inc-stepwise-decay.jl")
    include("module-Cascade-inc-simulations.jl")
    

end # module


