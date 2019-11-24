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
    `struct  Atomic.MeanFieldBasis`  
        ... a struct to represent (and generate) a mean-field basis and, especially, a set of (mean-field) orbitals.

        + levelsScf            ::Array{Int64,1}     ... Levels on which the optimization need to be carried out.
    """
    struct         MeanFieldBasis  <:  AbstractRepresentationType
        levelsScf              ::Array{Int64,1}
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
    `struct  Atomic.CiExpansion`  
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


     

    #==
    ### Green function representation ##################################################################################
    """
    `struct  Atomic.GreenFunction`  
        ... a struct to represent (and generate) an approximate (many-electron) Green functions; this struct enables one to specify further
            about the excitation scheme and the particular approximation used.

        + levelsScf            ::Array{Int64,1}     ... Levels on which the optimization need to be carried out.
    """
    struct         GreenFunction   <:  AbstractRepresentationType
        levelsScf              ::Array{Int64,1}
     end
    ==#



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
    
