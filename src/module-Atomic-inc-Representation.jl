    # This file is 'included' into the module Atomic

    
    """
    `abstract type Atomic.AbstractRepresentationType` 
        ... defines an abstract type and a number of data types to work with and distinguish different atomic representations; see also:
        
        + struct MeanFieldBasis       ... to represent (and generate) a mean-field basis and, especially, a set of (mean-field) orbitals.
        + struct CiExpansion          ... to represent (and generate) a configuration-interaction representation.
        + struct RasExpansion         ... to represent (and generate) a restricted active-space representation.
        + struct GreenFunction        ... to represent (and generate) an approximate (many-electron) Green functions.
    """
    abstract type  AbstractRepresentationType                           end

    
    """
    `struct  Atomic.MeanFieldSettings`  
        ... a struct for defining the settings for a mean-field basis (orbital) representation.

        + methodScf            ::String             ... Specify the SCF method: ["meanDFS", "meanHS"].
    """
    struct  MeanFieldSettings
        methodScf              ::String
    end

    """
    `Atomic.MeanFieldSettings()`  ... constructor for setting the default values.
    """
    function MeanFieldSettings()
    	MeanFieldSettings("meanDFS")
    end
    
    
    # `Base.show(io::IO, settings::MeanFieldSettings)`  ... prepares a proper printout of the settings::MeanFieldSettings.
    function Base.show(io::IO, settings::MeanFieldSettings)
    	  println(io, "methodScf:                $(settings.methodScf)  ")
    end


    """
    `struct  Atomic.MeanFieldBasis  <:  AbstractRepresentationType`  
        ... a struct to represent (and generate) a mean-field orbital basis.

        + settings         ::Atomic.MeanFieldSettings      ... Settings for the given mean-field orbital basis
    """
    struct   MeanFieldBasis  <:  AbstractRepresentationType
        settings           ::Atomic.MeanFieldSettings
    end


    # `Base.string(basis::MeanFieldBasis)`  ... provides a String notation for the variable basis::MeanFieldBasi.
    function Base.string(basis::MeanFieldBasis)
        sa = "Mean-field orbital basis for a $(basis.settings.methodScf) SCF field:"
        return( sa )
    end


    # `Base.show(io::IO, basis::MeanFieldBasis)`  ... prepares a proper printout of the basis::MeanFieldBasis.
    function Base.show(io::IO, basis::MeanFieldBasis)
        sa = Base.string(basis);       print(io, sa, "\n")
    end

    
    
    ### RAS: Restricted-Active-Space expansions ##################################################################################
    """
    `struct  Atomic.RasSettings`  
        ... a struct for defining the settings for a restricted active-space computations.

        + levelsScf            ::Array{Int64,1}     ... Levels on which the optimization need to be carried out.
        + maxIterationsScf     ::Int64              ... maximum number of SCF iterations in each RAS step.
        + accuracyScf          ::Float64            ... convergence criterion for the SCF field.
        
    	+ breitCI              ::Bool               ... logical flag to include Breit interactions.
    	+ selectLevelsCI       ::Bool               ... true, if specific level (number)s have been selected.
    	+ selectedLevelsCI     ::Array{Int64,1}     ... Level number that have been selected.
    """
    struct  RasSettings
        levelsScf              ::Array{Int64,1}
        maxIterationsScf       ::Int64  
        accuracyScf            ::Float64 
    	breitCI                ::Bool 
    	selectLevelsCI         ::Bool 
    	selectedLevelsCI       ::Array{Int64,1} 
    end

    """
    `Atomic.RasSettings()`  ... constructor for setting the default values.
    """
    function RasSettings()
    	RasSettings(Int64[1], 24, 1.0e-6, false, false, Int64[])
    end
    
    
    # `Base.show(io::IO, settings::RasSettings)`  ... prepares a proper printout of the settings::RasSettings.
    function Base.show(io::IO, settings::RasSettings)
    	  println(io, "levelsScf:            $(settings.levelsScf)  ")
    	  println(io, "maxIterationsScf:     $(settings.maxIterationsScf)  ")
    	  println(io, "accuracyScf:          $(settings.accuracyScf)  ")
    	  println(io, "breitCI:              $(settings.breitCI)  ")
    	  println(io, "selectLevelsCI:       $(settings.selectLevelsCI)  ")
    	  println(io, "selectedLevelsCI:     $(settings.selectedLevelsCI)  ")
    end


    """
    `struct  Atomic.RasStep`  
        ... specifies an individual step of a (relativistic) restricted active space computation for a set of levels. This struct 
            comprises all information to generate the orbital basis and to perform the associated SCF and multiplet computations for a 
            selected number of levels.

        + seFrom            ::Array{Shell,1}        ... Single-excitations from shells   [sh_1, sh_2, ...]
        + seTo              ::Array{Shell,1}        ... Single-excitations to shells  [sh_1, sh_2, ...]
        + deFrom            ::Array{Shell,1}        ... Double-excitations from shells   [sh_1, sh_2, ...]
        + deTo              ::Array{Shell,1}        ... Double-excitations to shells  [sh_1, sh_2, ...]
        + teFrom            ::Array{Shell,1}        ... Triple-excitations from shells   [sh_1, sh_2, ...]
        + teTo              ::Array{Shell,1}        ... Triple-excitations to shells  [sh_1, sh_2, ...]
        + qeFrom            ::Array{Shell,1}        ... Quadrupole-excitations from shells   [sh_1, sh_2, ...]
        + qeTo              ::Array{Shell,1}        ... Quadrupole-excitations to shells  [sh_1, sh_2, ...]
        + frozenShells      ::Array{Shell,1}        ... List of shells that are kept 'frozen' in this step.
        + constraints       ::Array{String,1}       ... List of Strings to define 'constraints/restrictions' to the generated CSF basis.
    """
    struct  RasStep
        seFrom              ::Array{Shell,1}
        seTo                ::Array{Shell,1}
        deFrom              ::Array{Shell,1}
        deTo                ::Array{Shell,1}
        teFrom              ::Array{Shell,1}
        teTo                ::Array{Shell,1}
        qeFrom              ::Array{Shell,1}
        qeTo                ::Array{Shell,1}
        frozenShells        ::Array{Shell,1}
        constraints         ::Array{String,1}
    end

    """
    `Atomic.RasStep()`  ... constructor for an 'empty' instance of a variable::Atomic.RasStep
    """
    function RasStep()
        RasStep(Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[],    Shell[], String[])
    end


    """
    `Atomic.RasStep(rasStep::Atomic.RasStep;`
    
                        seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[], 
                        deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[], 
                        teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[], 
                        qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[], 
                        frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[]  
                        
        ... constructor for modifying the given rasStep by specifying all excitations, frozen shells and constraints optionally.
    """
    function RasStep(rasStep::Atomic.RasStep;
                     seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[], 
                     deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[], 
                     teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[], 
                     qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[], 
                     frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[])
        if  seFrom == Shell[]   sxFrom = rasStep.seFrom   else      sxFrom = seFrom      end
        if  seTo   == Shell[]   sxTo   = rasStep.seTo     else      sxTo   = seTo        end
        if  deFrom == Shell[]   dxFrom = rasStep.deFrom   else      dxFrom = deFrom      end
        if  deTo   == Shell[]   dxTo   = rasStep.deTo     else      dxTo   = deTo        end
        if  teFrom == Shell[]   txFrom = rasStep.teFrom   else      txFrom = teFrom      end
        if  teTo   == Shell[]   txTo   = rasStep.teTo     else      txTo   = teTo        end
        if  qeFrom == Shell[]   qxFrom = rasStep.qeFrom   else      qxFrom = qeFrom      end
        if  qeTo   == Shell[]   qxTo   = rasStep.qeTo     else      qxTo   = qeTo        end
        if  frozen == Shell[]        frozx  = rasStep.frozenShells   else      frozx  = frozen             end
        if  constraints ==String[]   consx  = rasStep.constraints    else      consx  = constraints        end
        
        RasStep( sxFrom, sxTo, dxFrom, dxTo, txFrom, txTo, qxFrom, qxTo, frozx, consx)
    end


    # `Base.string(step::Atomic.RasStep)`  ... provides a String notation for the variable step::Atomic.RasStep.
    function Base.string(step::Atomic.RasStep)
        sa = "\nCI or RAS step with $(length(step.frozenShells)) (explicitly) frozen shell(s): $(step.frozenShells)  ... and virtual excitations"
        return( sa )
    end


    # `Base.show(io::IO, step::Atomic.RasStep)`  ... prepares a proper printout of the (individual step of computations) step::Atomic.RasStep.
    function Base.show(io::IO, step::Atomic.RasStep)
        sa = Base.string(step);                 print(io, sa, "\n")
        if  length(step.seFrom) > 0
           sa = "   Singles from:          { ";      for  sh in step.seFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.seTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
        if  length(step.deFrom) > 0
           sa = "   Doubles from:          { ";      for  sh in step.deFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.deTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
        if  length(step.teFrom) > 0
           sa = "   Triples from:          { ";      for  sh in step.teFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.teTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
        if  length(step.qeFrom) > 0
           sa = "   Quadruples from:       { ";      for  sh in step.qeFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.qeTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
    end
    
    
    """
    `struct  Atomic.RasExpansion    <:  AbstractRepresentationType`  
        ... a struct to represent (and generate) a restricted active-space representation.

        + symmetry         ::LevelSymmetry             ... Symmetry of the levels/CSF in the many-electron basis.
        + NoElectrons      ::Int64                     ... Number of electrons.
        + steps            ::Array{Atomic.RasStep,1}   ... List of SCF steps that are to be done in this model computation.
        + settings         ::Atomic.RasSettings        ... Settings for the given RAS computation
    """
    struct         RasExpansion    <:  AbstractRepresentationType
        symmetry           ::LevelSymmetry
        NoElectrons        ::Int64
        steps              ::Array{Atomic.RasStep,1}
        settings           ::Atomic.RasSettings 
     end


    """
    `Atomic.RasExpansion()`  ... constructor for an 'empty' instance of the a variable::Atomic.RasExpansion
    """
    function RasExpansion()
        RasExpansion(Basics.LevelSymmetry(0, Basics.plus), 0, 0, Atomic.RasStep[], Atomic.RasSettings())
    end


    # `Base.string(expansion::RasExpansion)`  ... provides a String notation for the variable expansion::RasExpansion.
    function Base.string(expansion::RasExpansion)
        sa = "RAS expansion for symmétry $(expansion.symmetry) and with $(length(expansion.steps)) steps:"
        return( sa )
    end


    # `Base.show(io::IO, expansion::RasExpansion)`  ... prepares a proper printout of the (individual step of computations) expansion::RasExpansion.
    function Base.show(io::IO, expansion::RasExpansion)
        sa = Base.string(expansion);       print(io, sa, "\n")
        println(io, "$(expansion.steps)  ")
        println(io, "... and the current settings:")
        println(io, "$(expansion.settings)  ")
    end

    
    
    ### CI: Configuration Interaction expansions ##################################################################################
    """
    `struct  Atomic.CiSettings`  
        ... a struct for defining the settings for a configuration-interaction (CI) expansion.

    	+ breitCI              ::Bool               ... logical flag to include Breit interactions.
    	+ selectLevelsCI       ::Bool               ... true, if specific level (number)s have been selected.
    	+ selectedLevelsCI     ::Array{Int64,1}     ... Level number that have been selected.
    	+ selectSymmetriesCI   ::Bool                          ... true, if specific level symmetries have been selected.
    	+ selectedSymmetriesCI ::Array{LevelSymmetry,1}        ... Level symmetries that have been selected.
    """
    struct  CiSettings
    	breitCI                ::Bool 
    	selectLevelsCI         ::Bool 
    	selectedLevelsCI       ::Array{Int64,1} 
    	selectSymmetriesCI     ::Bool 	
    	selectedSymmetriesCI   ::Array{LevelSymmetry,1} 
    end

    """
    `Atomic.CiSettings()`  ... constructor for setting the default values.
    """
    function CiSettings()
    	CiSettings(false, false, Int64[], false, LevelSymmetry[])
    end


    """
    `Atomic.CiSettings(settings::Atomic.CiSettings;`
    
            breitCI::Union{Nothing,Bool}=nothing, 
            selectLevelsCI::Union{Nothing,Bool}=nothing,                    selectedLevelsCI::Union{Nothing,Array{Int64,1}}=nothing, 
            selectSymmetriesCI::Union{Nothing,Bool}=nothing,                selectedSymmetriesCI::Union{Nothing,Array{LevelSymmetry,1}}=nothing)
                        
        ... constructor for modifying the given CiSettings by 'overwriting' the explicitly selected parameters.
    """
    function CiSettings(settings::Atomic.CiSettings;
        breitCI::Union{Nothing,Bool}=nothing, 
        selectLevelsCI::Union{Nothing,Bool}=nothing,                        selectedLevelsCI::Union{Nothing,Array{Int64,1}}=nothing, 
        selectSymmetriesCI::Union{Nothing,Bool}=nothing,                    selectedSymmetriesCI::Union{Nothing,Array{LevelSymmetry,1}}=nothing)
        
        if  breitCI             == nothing   breitCIx              = settings.breitCI               else   breitCIx              = breitCI              end 
        if  selectLevelsCI      == nothing   selectLevelsCIx       = settings.selectLevelsCI        else   selectLevelsCIx       = selectLevelsCI       end 
        if  selectedLevelsCI    == nothing   selectedLevelsCIx     = settings.selectedLevelsCI      else   selectedLevelsCIx     = selectedLevelsCI     end 
        if  selectSymmetriesCI  == nothing   selectSymmetriesCIx   = settings.selectSymmetriesCI    else   selectSymmetriesCIx   = selectSymmetriesCI   end 
        if  selectedSymmetriesCI== nothing   selectedSymmetriesCIx = settings.selectedSymmetriesCI  else   selectedSymmetriesCIx = selectedSymmetriesCI end 
        
        CiSettings( breitCIx, selectLevelsCIx, selectedLevelsCIx, selectSymmetriesCIx, selectedSymmetriesCIx)
    end
    
    
    # `Base.show(io::IO, settings::CiSettings)`  ... prepares a proper printout of the settings::CiSettings.
    function Base.show(io::IO, settings::CiSettings)
    	  println(io, "breitCI:                  $(settings.breitCI)  ")
    	  println(io, "selectLevelsCI:           $(settings.selectLevelsCI)  ")
    	  println(io, "selectedLevelsCI:         $(settings.selectedLevelsCI)  ")
    	  println(io, "selectSymmetriesCI:       $(settings.selectSymmetriesCI)  ")
    	  println(io, "selectedSymmetriesCI:     $(settings.selectedSymmetriesCI)  ")
    end


    """
    `struct  Atomic.CiExpansion  <:  AbstractRepresentationType`  
        ... a struct to represent (and generate) a configuration-interaction representation.

        + applyOrbitals    ::Dict{Subshell, Orbital}
        + excitations      ::Atomic.RasStep            ... Excitations beyond refConfigs.
        + settings         ::Atomic.CiSettings         ... Settings for the given CI expansion
    """
    struct         CiExpansion     <:  AbstractRepresentationType
        applyOrbitals      ::Dict{Subshell, Orbital}
        excitations        ::Atomic.RasStep
        settings           ::Atomic.CiSettings
    end


    # `Base.string(expansion::CiExpansion)`  ... provides a String notation for the variable expansion::CiExpansion.
    function Base.string(expansion::CiExpansion)
        sa = "CI expansion with (additional) excitations:"
        return( sa )
    end


    # `Base.show(io::IO, expansion::CiExpansion)`  ... prepares a proper printout of the (individual step of computations) expansion::CiExpansion.
    function Base.show(io::IO, expansion::CiExpansion)
        sa = Base.string(expansion);       print(io, sa, "\n")
        println(io, "$(expansion.excitations)  ")
        println(io, "... and the current settings:")
        println(io, "$(expansion.settings)  ")
    end


     
    ### Green function representation ##################################################################################

    """
    `abstract type Atomic.AbstractGreenApproach` 
        ... defines an abstract and a number of singleton types for approximating a many-electron Green
            function expansion for calculating second-order processes.

      + struct SingleCSFwithoutCI        
        ... to approximate the many-electron multiplets (gMultiplet) for every chosen level symmetry by dealing with each CSF 
            independently, and without any configuration interaction. This is a fast but also very rough approximation.
            
      + struct CoreSpaceCI                
        ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction between the 
            bound-state orbitals into account.
            
      + struct DampedSpaceCI                
        ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction for all, the 
            bound and free-electron, orbitals into account but by including a damping factor e^{tau*r}, tau > 0 into the 
            electron densities rho_ab (r) --> rho_ab (r) * e^{tau*r}
    """
    abstract type  AbstractGreenApproach                                  end
    struct         SingleCSFwithoutCI  <:  Atomic.AbstractGreenApproach   end
    struct         CoreSpaceCI         <:  Atomic.AbstractGreenApproach   end
    struct         DampedSpaceCI       <:  Atomic.AbstractGreenApproach   end

   
    """
    `struct  Atomic.GreenSettings`  
        ... defines a type for defining the details and parameters of the approximate Green (function) expansion.

        + nMax                     ::Int64            ... maximum principal quantum numbers of (single-electron) 
                                                          excitations that are to be included into the representation.
        + lValues                  ::Array{Int64,1}   ... List of (non-relativistic) orbital angular momenta for which
                                                          (single-electron) excitations are to be included.
        + printBeforeComputation   ::Bool             ... True if a short overview is to be printed before. 
        + selectLevels             ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels           ::Array{Int64,1}   ... List of selected levels.
    """
    struct GreenSettings 
        nMax                       ::Int64
        lValues                    ::Array{Int64,1}
        printBeforeComputation     ::Bool 
        selectLevels               ::Bool
        selectedLevels             ::Array{Int64,1}
    end 


    """
    `Atomic.GreenSettings()`  ... constructor for an `empty` instance of Atomic.GreenSettings.
    """
    function GreenSettings()
        Settings( 0, Int64[], LevelSymmetry[], false, false, Int64[])
    end


    # `Base.show(io::IO, settings::Atomic.GreenSettings)`  ... prepares a proper printout of the variable settings::Atomic.GreenSettings.
    function Base.show(io::IO, settings::Atomic.GreenSettings) 
        println(io, "nMax:                     $(settings.nMax)  ")
        println(io, "lValues:                  $(settings.lValues)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    """
    `struct  Atomic.GreenChannel`  
        ... defines a type for a single symmetry channel of an (approximate) Green (function) expansion.

        + symmetry          ::LevelSymmetry    ... Level symmetry of this part of the representation.
        + gMultiplet        ::Multiplet        ... Multiplet of (scattering) levels of this symmetry.
    """
    struct GreenChannel 
        symmetry            ::LevelSymmetry
        gMultiplet          ::Multiplet
    end   


    """
    `Atomic.GreenChannel()`  ... constructor for an `empty` instance of Atomic.GreenChannel.
    """
    function GreenChannel()
        GreenChannel( LevelSymmetry(0, Basics.plus), ManyElectron.Multiplet)
    end


    # `Base.show(io::IO, channel::Atomic.GreenChannel)`  ... prepares a proper printout of the variable channel::Atomic.GreenChannel.
    function Base.show(io::IO, channel::Atomic.GreenChannel) 
        println(io, "symmetry:                $(channel.symmetry)  ")
        println(io, "gMultiplet:              $(channel.gMultiplet)  ")
    end


    """
    `struct  Atomic.GreenExpansion  <:  AbstractRepresentationType`  
        ... defines a type to keep an (approximate) Green (function) expansion that is associated with a given set of reference
            configurations.

        + approach          ::Atomic.AbstractGreenApproach    ... Approach used to approximate the representation.
        + excitationScheme  ::Basics.AbstractExcitationScheme ... Applied excitation scheme w.r.t. refConfigs. 
        + levelSymmetries   ::Array{LevelSymmetry,1}          ... Total symmetries J^P to be included into Green expansion.
        + NoElectrons       ::Int64                           ... Number of electrons.
        + settings          ::Atomic.GreenSettings            ... settings for the Green (function) expansion.
    """
    struct GreenExpansion  <:  AbstractRepresentationType
        approach            ::Atomic.AbstractGreenApproach
        excitationScheme    ::Basics.AbstractExcitationScheme 
        levelSymmetries     ::Array{LevelSymmetry,1}
        NoElectrons         ::Int64 
        settings            ::Atomic.GreenSettings
    end   


    """
    `Atomic.GreenExpansion()`  ... constructor for an `empty` instance of Atomic.GreenExpansion.
    """
    function GreenExpansion()
        GreenExpansion( Atomic.SingleCSFwithoutCI(), Basics.NoExcitationScheme(), LevelSymmetry[], 0, Atomic.GreenSettings())
    end


    # `Base.string(expansion::GreenExpansion)`  ... provides a String notation for the variable expansion::GreenExpansion.
    function Base.string(expansion::GreenExpansion)
        sa = "Green (function) expansion in $(expansion.approach) approach and for excitation scheme  $(expansion.excitationScheme)," *
             "\nincluding (Green function channels with) symmetries $(expansion.levelSymmetries):"
        return( sa )
    end


    # `Base.show(io::IO, expansion::GreenExpansion)`  ... prepares a proper printout of the (individual step of computations) expansion::GreenExpansion.
    function Base.show(io::IO, expansion::GreenExpansion)
        sa = Base.string(expansion);       print(io, sa, "\n")
        println(io, "... and the current settings:")
        println(io, "$(expansion.settings)  ")
    end



    ### Representation ##################################################################################
    """
    `struct  Atomic.Representation`  
        ... a struct for defining an atomic representation. Such representations often refer to approximate wave function approximations of
            one or several levels but may concern also a mean-field basis (for some multiplet of some given configurations) or Green functions,
            etc.

        + name             ::String                      ... to assign a name to the given model.
        + nuclearModel     ::Nuclear.Model               ... Model, charge and parameters of the nucleus.
        + grid             ::Radial.Grid                 ... The radial grid to be used for the computation.
        + refConfigs       ::Array{Configuration,1}      ... List of references configurations, at least 1.
        + repType          ::AbstractRepresentationType  ... Specifies the particular representation.
    """
    struct  Representation
        name               ::String  
        nuclearModel       ::Nuclear.Model
        grid               ::Radial.Grid
        refConfigs         ::Array{Configuration,1} 
        repType            ::AbstractRepresentationType
    end


    """
    `Atomic.Representation()`  ... constructor for an 'empty' instance of the a variable::Atomic.Representation
    """
    function Representation()
        Representation("", Nuclear.Model(1.0), Radial.Grid(), ManyElectron.Configuration[], CiExpansion())
    end


    # `Base.string(rep::Representation)`  ... provides a String notation for the variable rep::Representation.
    function Base.string(rep::Representation)
        sa = "Atomic representation:   $(rep.name) for Z = $(rep.nuclearModel.Z) and with reference configurations: \n   "
        for  refConfig  in  rep.refConfigs     sa = sa * string(refConfig) * ",  "     end
        return( sa )
    end


    # `Base.show(io::IO, rep::Representation)`  ... prepares a printout of rep::Representation.
    function Base.show(io::IO, rep::Representation)
        sa = Base.string(rep);            print(io, sa, "\n")
        println(io, "representation type:   $(rep.repType)  ")
        println(io, "nuclearModel:          $(rep.nuclearModel)  ")
        println(io, "grid:                  $(rep.grid)  ")
    end
    
