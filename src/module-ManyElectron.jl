
"""
`module  JAC.ManyElectron`  
    ... a submodel of JAC that contains all structs and methods for many-electron CSF, basis and levels.
"""
module ManyElectron

    using  ..Basics, ..Defaults,  ..Radial
    export Configuration, ConfigurationR, AsfSettings, CsfR, Basis, Level, Multiplet, 
           AbstractQedModel, NoneQed, QedPetersburg, QedSydney, LSjjSettings,
           AbstractConfigurationRestriction, RestrictNoElectronsTo, RestrictParity, RestrictToShellDoubles, RequestMinimumOccupation, RequestMaximumOccupation

    
    
    """
    `abstract type ManyElectron.AbstractQedModel` 
        ... defines an abstract and a number of singleton types for dealing with QED estimates in the many-electron
            computations.

      + struct QedPetersburg         
        ... to estimate the QED corrections due to the model Hamiltonian by Shabaev and coworkers (2013).
            
      + struct QedSydney                
        ... to estimate the QED corrections by means of the radiative potential by Flambaum and Ginges (2004).
            
      + struct NoneQed               
        ... No QED estimates are included into the computations.
    """
    abstract type  AbstractQedModel                      end
    struct         QedPetersburg  <:  AbstractQedModel   end
    struct         QedSydney      <:  AbstractQedModel   end
    struct         NoneQed        <:  AbstractQedModel   end

    
    
    """
    `abstract type ManyElectron.AbstractConfigurationRestriction` 
        ... defines an abstract types for dealing with restrictions that need to be applied to a list of configurations.
            Typically, a loop through is made through all given restrictions and all configurations are tested to obey all
            these restrictions. Two contradicting restrictions, for instance RestrictParity(plus) & RestrictParity(minus),
            therefore leads zero configurations in all cases. It remains the reponsibility of the user to make sure that
            the given restrictions are consistent with what is to be achieved. The given set of restrictions can be easily
            extended if this need arises by the users.

        + RestrictNoElectronsTo(..)    ... to restrict the total number of electron in high subshells.
        + RestrictParity(..)           ... to restrict to configurations with a given parity.
        + RestrictToShellDoubles(..)   ... to allow only double occupations in high subshells.
        + RequestMinimumOccupation(..) ... to request a minimum occupation in a given set of shells.
        + RequestMaximumOccupation(..) ... to request a maximum occupation in a given set of shells.
    """
    abstract type  AbstractConfigurationRestriction      end


    """
    `struct  ManyElectron.RestrictNoElectronsTo  <: AbstractConfigurationRestriction`   
        ... restrict the number of electron in all shells with principal quantum number n >= nmin or orbital angular momentum l >= lmin 
            to a total of ne electrons.

        + ne            ::Int64     ... maximum number of (allowed) electrons in the specicied higher subshells.
        + nmin          ::Int64     ... principal quantum number nmin.
        + lmin          ::Int64     ... orbital angular momentum lmin.
    """
    struct   RestrictNoElectronsTo  <: AbstractConfigurationRestriction
        ne              ::Int64
        nmin            ::Int64
        lmin            ::Int64
    end


    function Base.string(res::RestrictNoElectronsTo)
        sa = "Restrict to configurations with a maximum of $(res.ne) electrons in shells with n >= $(res.nmin) & l >= $(res.lmin)."
        return( sa )
    end

    function Base.show(io::IO, res::RestrictNoElectronsTo)
        sa = string(res);       print(io, sa)
    end


    """
    `struct  ManyElectron.RestrictParity  <: AbstractConfigurationRestriction`   
        ... restrict to configurations with a given parity.

        + parity        ::Basics.Parity   ... given parity.
    """
    struct   RestrictParity  <: AbstractConfigurationRestriction
        parity          ::Basics.Parity
    end


    function Base.string(res::RestrictParity)
        sa = "Restrict to configurations with parity $(res.parity)."
        return( sa )
    end

    function Base.show(io::IO, res::RestrictParity)
        sa = string(res);       print(io, sa)
    end


    """
    `struct  ManyElectron.RestrictToShellDoubles  <: AbstractConfigurationRestriction`   
        ... restrict to a double electron occupation in all shells with principal quantum number n >= nmin or orbital angular momentum l >= lmin. 

        + nmin          ::Int64     ... principal quantum number nmin.
        + lmin          ::Int64     ... orbital angular momentum lmin.
    """
    struct   RestrictToShellDoubles  <: AbstractConfigurationRestriction
        nmin            ::Int64
        lmin            ::Int64
    end


    function Base.string(res::RestrictToShellDoubles)
        sa = "Restrict to configurations with a double electron occupation in shells with n >= $(res.nmin) & l >= $(res.lmin)."
        return( sa )
    end

    function Base.show(io::IO, res::RestrictToShellDoubles)
        sa = string(res);       print(io, sa)
    end


    """
    `struct  ManyElectron.RequestMinimumOccupation  <: AbstractConfigurationRestriction`   
        ... request a minimum occupation ne in the given (list of) shells. 

        + ne            ::Int64           ... minimum electron occupation.
        + shells        ::Array{Shell,1}  ... list of shells.
    """
    struct   RequestMinimumOccupation  <: AbstractConfigurationRestriction
        ne              ::Int64 
        shells          ::Array{Shell,1}
    end


    function Base.string(res::RequestMinimumOccupation)
        sa = "Request a minimum occupation of $(res.ne) electron in the (list of) shells $(res.shells)."
        return( sa )
    end

    function Base.show(io::IO, res::RequestMinimumOccupation)
        sa = string(res);       print(io, sa)
    end


    """
    `struct  ManyElectron.RequestMaximumOccupation  <: AbstractConfigurationRestriction`   
        ... request a maximum occupation ne in the given (list of) shells. 

        + ne            ::Int64           ... maximum electron occupation.
        + shells        ::Array{Shell,1}  ... list of shells.
    """
    struct   RequestMaximumOccupation  <: AbstractConfigurationRestriction
        ne              ::Int64 
        shells          ::Array{Shell,1}
    end


    function Base.string(res::RequestMaximumOccupation)
        sa = "Request a maximum occupation of $(res.ne) electron in the (list of) shells $(res.shells)."
        return( sa )
    end

    function Base.show(io::IO, res::RequestMaximumOccupation)
        sa = string(res);       print(io, sa)
    end

    

    """
    `struct  ManyElectron.Configuration`  
        ... defines a struct for a non-relativistic electron configuration that is fully speficied by its shell notations 
            (such as "1s", "2s", "2p", ....) and the corresponding occupation numbers (>= 0). An electron configuration 
            is independent of any order of the shells, and a zero occupation is assumed for all shells that do not appear as a 
            key in the shell (dictionary). Therefore, the number of keys do not allow any conclusion about the underlying orbital 
            space of any considered computation that include more than just a single configuration.

        + shells         ::Dict{Shell,Int64}   ... Dictionary that maps shells to their occupation.
        + NoElectrons    ::Int64               ... No. of electrons.
    """
    struct  Configuration
        shells           ::Dict{Shell,Int64}  
        NoElectrons      ::Int64               
    end


    """
    `ManyElectron.Configuration(sa::String)`  
        ... constructor for a given configuration string, such as "[He]", "[Ne]", "[Ne] 3s 3p^6"  or "1s 2p^6 3s^2 3p".
    """
    function Configuration(sa::String)
        sax = strip(sa)
        shellOccList = split(sax, " ")
        ##x println("shellOccList = $shellOccList")
        NoElectrons = 0;    shellDict    = Dict{Shell,Int64}()
        for shocc in shellOccList
            if       shocc == " "    continue
            elseif   shocc  in [ "[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]" ]
                shellDict = Base.merge( shellDict, Dict( Shell("1s") => 2));                                       NoElectrons = NoElectrons + 2
                shocc == "[He]"  &&  continue
                shellDict = Base.merge( shellDict, Dict(Shell("2s") => 2,  Shell("2p") => 6));                     NoElectrons = NoElectrons + 8  
                shocc == "[Ne]"  &&  continue
                shellDict = Base.merge( shellDict, Dict(Shell("3s") => 2,  Shell("3p") => 6));                     NoElectrons = NoElectrons + 8  
                shocc == "[Ar]"  &&  continue 
                shellDict = Base.merge( shellDict, Dict(Shell("3d") => 10, Shell("4s") => 2, Shell("4p") => 6));   NoElectrons = NoElectrons + 18   
                shocc == "[Kr]"  &&  continue   
                shellDict = Base.merge( shellDict, Dict(Shell("5s") => 10, Shell("5s") => 2, Shell("5p") => 6));   NoElectrons = NoElectrons + 18  
                shocc == "[Xe]"  &&  continue 
            else
                ia = findall(in("^"), shocc);    
                if     length(ia) == 0     shellDict = Base.merge( shellDict, Dict( Shell(shocc) => 1));           NoElectrons = NoElectrons + 1
                else   occ = parse(Int64, shocc[ia[1]+1:end] );    sh = Shell( shocc[1:ia[1]-1] )
                       shellDict = Base.merge( shellDict, Dict( sh => occ));                                       NoElectrons = NoElectrons + occ
                end
            end
        end
        return( Configuration(shellDict, NoElectrons) )
    end


    # `Base.show(io::IO, conf::Configuration)`  ... prepares a proper printout of the non-relativistic configuration conf::Configuration.
    function Base.show(io::IO, conf::Configuration)
        sa = "Configuration: " * string(conf)
        print(io, sa)
    end


    # `Base.unique(confs::Array{Configuration,1})`  ... return a unique list of configurations.
    function Base.unique(confs::Array{Configuration,1})
        confList       = Configuration[]  
        configurations = deepcopy(confs)
        for  confa in configurations
            addTo = true
            for confb in confList    if   confa == confb    addTo = false;    break     end     end
            if  addTo    push!(confList, confa)     end
        end
    
        return( confList )
    end


    # `Base.string(conf::Configuration)`  ... provides a String notation for the variable conf::Configuration.
    function Base.string(conf::Configuration)
        wa = keys(conf.shells);   va = values(conf.shells)
        wb = Defaults.getDefaults("ordered shell list: non-relativistic", 7)

        sa = ""
        for  k  in  wb
            if   k  in  wa
              occ = conf.shells[k]
              sa = sa * string(k) * "^$occ "
            end
        end

        return( sa )
    end


    """
    `Base.:(==)(confa::Configuration, confb::Configuration)`  
        ... compares (recursively) two non-relativistic configurations and return true if all subfields are equal, 
            and false otherwise
    """
    function  Base.:(==)(confa::Configuration, confb::Configuration)
        if   confa.NoElectrons  != confb.NoElectrons    return( false )    end
        wk = keys(confa.shells)
        for  k in wk
            if  !haskey(confb.shells, k)  ||   confa.shells[k] != confb.shells[k]    return( false )    end
        end
        return( true )
    end


    """
    `struct  ManyElectron.ConfigurationR`  
        ... defines a type for a relativistic electron configuration that is fully speficied by its subshell notations 
            (such as "1s_1/2", "2s_1/2", "2p_3/2", ....) and the corresponding occupation numbers (>= 0). An electron 
            configuration is independent of any order of the subshells, and a zero occupation is assumed for all shells that 
            do not appears as a key. Therefore, the number of keys do not allow any conclusion about the underlying orbital 
            space of any considered computation that include more than a single configuration.

        + subshells      ::Dict{Subshell,Int64}  ... Dictionary that maps subshells to their occupation.
        + NoElectrons    ::Int64                 ... No. of electrons.
    """
    struct  ConfigurationR
        subshells        ::Dict{Subshell,Int64}  
        NoElectrons      ::Int64               
    end


    # `Base.show(io::IO, conf::ConfigurationR)`  ... prepares a proper printout of the relativistic configuration conf::ConfigurationR.
    function Base.show(io::IO, conf::ConfigurationR)
        sa = "Configuration: " * string(conf)
        print(io, sa)
    end


    # `Base.string(conf::ConfigurationR)`  ... provides a String notation for the variable conf::ConfigurationR.
    function Base.string(conf::ConfigurationR)
        wa = keys(conf.subshells);   va = values(conf.subshells)
        wb = Defaults.getDefaults("ordered subshell list: relativistic", 7)

        sa = ""
        for  k  in  wb
           if   k  in  wa
               occ = conf.subshells[k]
               sa = sa * string(k) * "^$occ "
           end
        end

        return( sa )
    end


    """
    `Base.:(==)(confa::ConfigurationR, confb::ConfigurationR)`  
        ... compares (recursively) two relativistic configurations and return true if all subfields are equal, 
            and false otherwise
    """
    function  Base.:(==)(confa::ConfigurationR, confb::ConfigurationR)
        if   confa.NoElectrons  != confb.NoElectrons    return( false )    end
        wk = keys(confa.subshells)
        for  k in wk
            if  !haskey(confb.subshells, k)  ||   confa.subshells[k] != confb.subshells[k]    return( false )    end
        end
        return( true )
    end

    
    """
    `abstract type ManyElectron.AbstractStartOrbitals` 
        ... defines an abstract and a number of singleton types to determine the orbitals from which the SCF computations start.

        + struct StartFromHydrogenic   ... to start the SCF procedure from hydrogenic orbitals.        
        + struct StartFromPrevious     ... to start the SCF procedure from previously calculated orbitals.        
    """
    abstract type  AbstractStartOrbitals                          end
    struct     StartFromHydrogenic    <:  AbstractStartOrbitals   end
    
    export  AbstractStartOrbitals, StartFromHydrogenic, StartFromPrevious


    """
    `struct  ManyElectron.StartFromPrevious  <: ManyElectron.AbstractStartOrbitals`   ... to start from a set of previously calculated orbitals.

        + orbitals      ::Dict{Subshell, Orbital}  ... set of previously-calculated orbitals. 
    """
    struct   StartFromPrevious  <: ManyElectron.AbstractStartOrbitals
        orbitals        ::Dict{Subshell, Orbital}
    end


    function Base.string(from::ManyElectron.StartFromPrevious)
        sa = " "
        for  (sh,v) in from.orbitals  sa = sa * "  " * string(sh)   end
        return( sa )
    end

    function Base.show(io::IO, from::ManyElectron.StartFromPrevious)
        sa = "Start from previously-calculated orbitals for shells:\n" * string(from);       print(io, sa, "\n")
    end


    """
    `struct  ManyElectron.LSjjSettings`  ... defines a type for the details and parameters for performing jj-LS expansions.

        + makeIt         ::Bool            ... True, if the jj-LS expansion is to be made and false otherwise.
        + minWeight      ::Float64         ... minimum weight with which a (relativistic) CSF must contribute to (at least one) level.
        + printWeight    ::Float64         ... minimum weight of a nonrelativistic CSF to be printed out in the final expansion.
        + selectLevels   ::Bool            ... True, if levels are selected individually for the computations.
        + selectedLevels ::Array{Int64,1}  ... List of levels.
    """
    struct LSjjSettings 
        makeIt          ::Bool
        minWeight       ::Float64 
        printWeight     ::Float64
        selectLines     ::Bool
        selectedLines   ::Array{Int64,1}    
    end 


    """
    `LSjjSettings(makeIt::Bool)`  ... constructor for the default values of jj-LS transformations.
    """
    function LSjjSettings(makeIt::Bool)
        LSjjSettings(makeIt, 0.05, 0.1, false, Int64[])
    end


    # `Base.show(io::IO, settings::LSjjSettings)`  ... prepares a proper printout of the variable settings::LSjjSettings.
    function Base.show(io::IO, settings::LSjjSettings) 
        println(io, "makeIt:            $(settings.makeIt)  ")
        println(io, "minWeight:         $(settings.minWeight)  ")
        println(io, "printWeight:       $(settings.printWeight)  ")
        println(io, "selectLines:       $(settings.selectLines)  ")
        println(io, "selectedLines:     $(settings.selectedLines)  ")
    end


    """
    `struct  ManyElectron.AsfSettings`  
        ... a struct for defining the settings for the atomic state functions, i.e. the self-consistent-field (SCF) 
            and CI computations

        + generateScf          ::Bool                   ... True, if a SCF need to be generated, and false otherwise (frozen orbitals).
        + eeInteraction        ::AbstractEeInteraction  ... Specify the e-e interaction to be included into the SCF computations.
        + scField              ::AbstractScField        ... Specify the self-consistent field, for instance, Basics.ALField(), etc.
        + startScfFrom         ::AbstractStartOrbitals  ... Specify the orbitals to start the SCF computations
        + maxIterationsScf     ::Int64                  ... maximum number of SCF iterations
        + accuracyScf          ::Float64                ... convergence criterion for the SCF field.
        + shellSequenceScf     ::Array{Subshell,1}      ... Sequence of subshells to be optimized.
        + frozenSubshells      ::Array{Subshell,1}      ... Sequence of subshells to be kept frozen.
        
    	+ eeInteractionCI      ::AbstractEeInteraction  ... Specify the e-e interaction to be included into the CI computations.
    	+ qedModel             ::AbstractQedModel       ... model for estimating QED corrections {NoneQed(), QedPetersburg(), QedSydney()}.
    	+ methodCI             ::AbstractCImethod       ... method for diagonalizing the matrix.
    	+ jjLS                 ::LSjjSettings           ... settings to control a jj-LS transformation of atomic level, if requested.
    	
    	+ selectLevelsCI       ::Bool                   ... true, if specific level (number)s have been selected.
    	+ selectedLevelsCI     ::Array{Int64,1}         ... Level number that have been selected.
    	+ selectSymmetriesCI   ::Bool                   ... true, if specific level symmetries have been selected.
    	+ selectedSymmetriesCI ::Array{LevelSymmetry,1} ... Level symmetries that have been selected.
    """
    struct  AsfSettings
        generateScf            ::Bool 
        eeInteraction          ::AbstractEeInteraction 
        scField                ::AbstractScField 
        startScfFrom           ::AbstractStartOrbitals
        maxIterationsScf       ::Int64 
        accuracyScf            ::Float64   
        shellSequenceScf       ::Array{Subshell,1}
        frozenSubshells        ::Array{Subshell,1}
        #
    	eeInteractionCI        ::AbstractEeInteraction
    	qedModel               ::AbstractQedModel 	
    	methodCI               ::AbstractCImethod
    	jjLS                   ::LSjjSettings
    	selectLevelsCI         ::Bool 
    	selectedLevelsCI       ::Array{Int64,1}
    	selectSymmetriesCI     ::Bool 	
    	selectedSymmetriesCI   ::Array{LevelSymmetry,1} 
     end


    """
    `ManyElectron.AsfSettings()`  ... constructor for setting the default values.
    """
    function AsfSettings()
    	AsfSettings(false, CoulombInteraction(), Basics.DFSField(), StartFromHydrogenic(), 24, 1.0e-6, Subshell[], Subshell[],  
    	            CoulombInteraction(), NoneQed(), FullCIeigen(), LSjjSettings(false), false, Int64[], false, LevelSymmetry[] )
    end


    """
    `ManyElectron.AsfSettings(settings::AsfSettings;`
        
                generateScf=..,       eeInteraction=..,       scField=..,            startScfFrom=..,           maxIterationsScf=..,    
                accuracyScf=..,       shellSequenceScf=..,    frozenSubshells=..,    eeInteractionCI=..,        qedModel=..,           
                methodCI=..,          jjLS=..,   
                selectLevelsCI=..,    selectedLevelsCI=..,    selectSymmetriesCI=.., selectedSymmetriesCI=..,   printout::Bool=false)
        ... constructor for re-defining a settings::AsfSettings.
    """
    function AsfSettings(settings::AsfSettings; 
        generateScf::Union{Nothing,Bool}=nothing,                       eeInteraction::Union{Nothing,AbstractEeInteraction}=nothing,         
        scField::Union{Nothing,AbstractScField}=nothing,                startScfFrom::Union{Nothing,AbstractStartOrbitals}=nothing,
        maxIterationsScf::Union{Nothing,Int64}=nothing,                 accuracyScf::Union{Nothing,Float64}=nothing,     
        shellSequenceScf::Union{Nothing,Array{Subshell,1}}=nothing,     frozenSubshells::Union{Nothing,Array{Subshell,1}}=nothing, 
        eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,  qedModel::Union{Nothing,AbstractQedModel}=nothing,              
        methodCI::Union{Nothing,AbstractCImethod}=nothing,              jjLS::Union{Nothing,LSjjSettings}=nothing,  
        selectLevelsCI::Union{Nothing,Bool}=nothing,                    selectedLevelsCI::Union{Nothing,Array{Int64,1}}=nothing,         
        selectSymmetriesCI::Union{Nothing,Bool}=nothing,                selectedSymmetriesCI::Union{Nothing,Array{LevelSymmetry,1}}=nothing,         
        printout::Bool=false)

        if  generateScf         == nothing   generateScfx          = settings.generateScf           else   generateScfx          = generateScf          end 
        if  eeInteraction       == nothing   eeInteractionx        = settings.eeInteraction         else   eeInteractionx        = eeInteraction        end 
        if  scField             == nothing   scFieldx              = settings.scField               else   scFieldx              = scField              end 
        if  startScfFrom        == nothing   startScfFromx         = settings.startScfFrom          else   startScfFromx         = startScfFrom         end 
        if  maxIterationsScf    == nothing   maxIterationsScfx     = settings.maxIterationsScf      else   maxIterationsScfx     = maxIterationsScf     end 
        if  accuracyScf         == nothing   accuracyScfx          = settings.accuracyScf           else   accuracyScf           = accuracyScf          end 
        if  shellSequenceScf    == nothing   shellSequenceScfx     = settings.shellSequenceScf      else   shellSequenceScfx     = shellSequenceScf     end 
        if  frozenSubshells     == nothing   frozenSubshellsx      = settings.frozenSubshells       else   frozenSubshellsx      = frozenSubshells      end 
        if  eeInteractionCI     == nothing   eeInteractionCIx      = settings.eeInteractionCI       else   eeInteractionCIx      = eeInteractionCI      end 
        if  qedModel            == nothing   qedModelx             = settings.qedModel              else   qedModelx             = qedModel             end 
        if  methodCI            == nothing   methodCIx             = settings.methodCI              else   methodCI              = methodCI             end 
        if  jjLS                == nothing   jjLSx                 = settings.jjLS                  else   jjLSx                 = jjLS                 end 
        if  selectLevelsCI      == nothing   selectLevelsCIx       = settings.selectLevelsCI        else   selectLevelsCIx       = selectLevelsCI       end 
        if  selectedLevelsCI    == nothing   selectedLevelsCIx     = settings.selectedLevelsCI      else   selectedLevelsCIx     = selectedLevelsCI     end 
        if  selectSymmetriesCI  == nothing   selectSymmetriesCIx   = settings.selectSymmetriesCI    else   selectSymmetriesCIx   = selectSymmetriesCI   end 
        if  selectedSymmetriesCI== nothing   selectedSymmetriesCIx = settings.selectedSymmetriesCI  else   selectedSymmetriesCIx = selectedSymmetriesCI end 
        
    	AsfSettings(generateScfx, eeInteractionx, scFieldx, startScfFromx, maxIterationsScfx, accuracyScfx, 
    	            shellSequenceScfx, frozenSubshellsx, eeInteractionCIx, qedModelx, methodCIx, jjLSx, 
    	            selectLevelsCIx, selectedLevelsCIx, selectSymmetriesCIx, selectedSymmetriesCIx)
    end
    
    
    # `Base.show(io::IO, settings::AsfSettings)`  ... prepares a proper printout of the settings::AsfSettings.
    function Base.show(io::IO, settings::AsfSettings)
    	  println(io, "generateScf:          $(settings.generateScf)  ")
    	  println(io, "eeInteraction:        $(settings.eeInteraction)  ")
    	  println(io, "scFieldf:             $(settings.scFieldf)  ")
    	  println(io, "startScfFrom:         $(settings.startScfFrom)  ")
    	  println(io, "maxIterationsScf:     $(settings.maxIterationsScf)  ")
    	  println(io, "accuracyScf:          $(settings.accuracyScf)  ")
    	  println(io, "shellSequenceScf:     $(settings.shellSequenceScf)  ")
    	  println(io, "frozenSubshells:      $(settings.frozenSubshells)  ")
    	  #
    	  println(io, "eeInteractionCI:      $(settings.eeInteractionCI)  ")
    	  println(io, "qedModel :            $(settings.qedModel)  ")
    	  println(io, "methodCI:             $(settings.methodCI)  ")
    	  println(io, "jjLS :                $(settings.jjLS.makeIt)  ")
    	  println(io, "selectLevelsCI:       $(settings.selectLevelsCI)  ")
    	  println(io, "selectedLevelsCI:     $(settings.selectedLevelsCI)  ")
    	  println(io, "selectSymmetriesCI:   $(settings.selectSymmetriesCI)  ")
    	  println(io, "selectedSymmetriesCI: $(settings.selectedSymmetriesCI)  ")
    end
    
    
    # `Base.string(settings::AsfSettings)`  ... provides a String notation for the variable settings::AsfSettings.
    function Base.string(settings::AsfSettings)
    	  error("Not yet implemented.")
    	  sa = "Asf settings: maximum No. of iterations = $(settings.maxIterationsScf), accuracy = (settings.accuracyScf)"
    	  return( sa )
    end
    

    """
    `Base.:(==)(seta::AsfSettings, setb::AsfSettings)`  
        ... compares (recursively) two instances of AsfSettings and return true if all subfields are equal, 
            and false otherwise
    """
    function  Base.:(==)(seta::AsfSettings, setb::AsfSettings)
    
    	if  seta.generateScf      !=  setb.generateScf                return( false )    end
    	if  seta.eeInteraction    !=  setb.eeInteraction              return( false )    end
    	if  seta.scField          !=  setb.scField                    return( false )    end
    	if  seta.startScfFrom     !=  setb.startScfFrom               return( false )    end
    	#   startOrbitals
    	if  seta.maxIterationsScf !=  setb.maxIterationsScf           return( false )    end
    	if  seta.accuracyScf      !=  setb.accuracyScf                return( false )    end
    	if  seta.shellSequenceScf !=  setb.shellSequenceScf           return( false )    end
    	if  seta.frozenSubshells  !=  setb.frozenSubshells            return( false )    end
    	if  seta.eeInteractionCI  !=  setb.eeInteractionCI            return( false )    end
    	#   qedModel
    	#   methodCI
    	#   jjLS
    	
    	  return( true )
    end

    
    
    """
    `struct  ManyElectron.CsfR`  
        ... defines a type for a relativistic configuration state function (CSF) in terms of its orbitals (sequence of orbitals), 
            their occupation as well as the angular momenta, seniority and the coupling of the corresponding antisymmetric subshell states.
    
        + useStandardSubshells ::Bool                 ... Determines whether the subshell list has to be obtained from some outer structure
                                                          given by the subshells field below. 
        + J                    ::AngularJ64           ... Total angular momentum J.
        + parity               ::Parity               ... Total parity.
        + occupation           ::Array{Int64,1}       ... occupation of the orbitals with regard to the specified orbital list.
        + seniority            ::Array{Int64,1}       ... list of seniority values of the antisymmetric subshell states
        + subshellJ            ::Array{AngularJ64,1}  ... list of J-values of the antisymmetric subshell states.
        + subshellX            ::Array{AngularJ64,1}  ... intermediate X-values in the coupling of the antisymmetric subshell states.
        + subshells            ::Array{Subshell,1}    ... Explicitly given subshell list if useStandardSubshells == false.
    """
    struct  CsfR
        useStandardSubshells   ::Bool 
        J		               ::AngularJ64
        parity  	           ::Parity
        occupation	           ::Array{Int64,1}
        seniority	           ::Array{Int64,1}
        subshellJ	           ::Array{AngularJ64,1}
        subshellX	           ::Array{AngularJ64,1}
        subshells	           ::Array{Subshell,1}
    end
    
    
    """
    `ManyElectron.CsfR(J::AngularJ64, parity::Parity)`  
        ... simple constructor for given J and parity, and where standardOrbitals is set to false.
    """
    function CsfR(J::AngularJ64, parity::String)
    	  CsfR(false, J, parity, Int64[], Int64[], AngularJ64[], AngularJ64[], Subshell[] )
    end
    
    
    """
    `ManyElectron.CsfR(sa::String)`  
        ... simple constructor for given closed-shell occupation (such as "[Ne]"), and where standardOrbitals is set to true.
    """
    function CsfR(sa::String)
    	  wa = subshellsFromClosedShellConfiguration(sa)
    
    	  occupation = Int64[];	  subshellJ = AngularJ64[];	subshellX  = AngularJ64[]
    	  subshells  = Subshell[];  seniority = Int64[]
    	
        for  sb  in  wa
    	      #x wb = Basics.SubshellQuantumNumbers(string(sb));   kappa = wb[2];   occ = wb[4] + 1
    	      occ = Basics.subshell_2j(sb) + 1
    	      occupation = push!(occupation, occ)
    	      seniority  = push!(seniority,    0) 
    	      subshellJ  = push!(subshellJ, AngularJ64(0) ) 
    	      subshellX  = push!(subshellX, AngularJ64(0) ) 
    	      subshells  = push!(subshells, sb ) 
    	  end
    	  J = AngularJ64(0);    parity = plus
    
    	  CsfR(false, J, parity, occupation, seniority, subshellJ, subshellX, subshells)
    end
    
    
    # `Base.show(io::IO, csf::CsfR)`  ... prepares a proper printout of the variable orbital::CsfR.
    function Base.show(io::IO, csf::CsfR) 
    	  print(io, string(csf) ) 
    end
    
    
    # `Base.string(csf::CsfR)`  ... provides a String notation for csf::CsfR.
    function Base.string(csf::CsfR) 
    	  sa = "\n   CSF: "
        if    csf.useStandardSubshells

    	      for  i in 1:length(csf.occupation)
    		    sa = sa * Basics.subshellStateString(string(Defaults.GBL_STANDARD_SUBSHELL_LIST[i]), csf.occupation[i], csf.seniority[i], 
                                                  csf.subshellJ[i], csf.subshellX[i])
    		    sa = sa * ", "
    	      end
    	      sa = sa * ": J=" * string(csf.J) * string(csf.parity)

            else
    	      for  i in 1:length(csf.subshells)
    		    sa = sa * Basics.subshellStateString(string(csf.subshells[i]), csf.occupation[i], csf.seniority[i], 
                                                  csf.subshellJ[i], csf.subshellX[i])
    		    sa = sa * ", "
    	      end
    	      sa = sa * ": J=" * string(csf.J) * string(csf.parity)
    	  end
    	  return( sa )
    end
    

    """
    `Base.:(==)(csfa::CsfR, csfb::CsfR)`  
        ... compares (recursively) two relativistic CSFs and return true if all subfields are equal, 
            and false otherwise
    """
    function  Base.:(==)(csfa::CsfR, csfb::CsfR)
        # Stop the comparison with an error message if the CSF are not defined w.r.t a common subshell list.
        if   !(csfa.useStandardSubshells)   ||   !(csfb.useStandardSubshells)	## ||   length(csfa.occupation) != length(csfb.occupation)
            error("To determine the equivalence of two CSF requires that both are defined on a common subshell list.")
        end
    
        if  length(csfa.occupation) != length(csfb.occupation)      return( false )    end
        if  csfa.J  !=  csfb.J  ||  csfa.parity  !=  csfb.parity	return( false )    end
    	if  csfa.occupation      !=  csfb.occupation			    return( false )    end
        if  csfa.seniority       !=  csfb.seniority			        return( false )    end
        if  csfa.subshellJ       !=  csfb.subshellJ			        return( false )    end
        if  csfa.subshellX       !=  csfb.subshellX			        return( false )    end
    	
    	  return( true )
    end
    
    
    """
    `ManyElectron.CsfRGrasp92(subshells::Array{Subshell,1}, coreSubshells::Array{Subshell,1}, sa::String, sb::String, sc::String)`  
        ... to construct a CsfR from 3 Grasp-Strings and the given list of subshells and core subshells.
    """
    function CsfRGrasp92(subshells::Array{Subshell,1}, coreSubshells::Array{Subshell,1}, sa::String, sb::String, sc::String)
    	occupation = Int64[];	 seniority = Int64[];	 subshellJ = AngularJ64[];    subshellX = AngularJ64[];    subshellsx = Subshell[]

    	# Define the occupation and quantum numbers of the closed shells
    	for i in 1:length(coreSubshells)
    	    push!(occupation, Basics.subshell_2j( coreSubshells[i] ) +1);    
    	    push!(seniority, 0);     push!(subshellJ, AngularJ64(0));	 push!(subshellX, AngularJ64(0)) 
    	end
    
    	scc = sc * "	  ";   scx = " "
    	i  = -1;   ishell = length(coreSubshells)
    	while true
    	    i = i + 1;    sax = sa[9i+1:9i+9];     sbx = sb[9i+1:9i+9];    scx = strip( scc[9i+6:9i+14] )
    
    	    sh = strip( sax[1:5] );    sh = Basics.subshellGrasp(sh);    occ = parse( sax[7:8] )
    	    search(sbx, ',') > 0    &&    error("stop a: missing seniority; sb = $sb")   # Include seniority more properly if it occurs
    	    subJx = parse( sbx )
    	    if      typeof(subJx) == Void     subJ = AngularJ64(0)
    	    elseif  typeof(subJx) == Int64    subJ = AngularJ64(subJx)
    	    elseif  typeof(subJx) == Expr     subJ = AngularJ64( subJx.args[2]// subJx.args[3] )
    	    else    error("stop b; type = $(typeof(subJx))")
    	    end
    
    	    if   length(scx) > 0  &&  scx[end:end] in ["+", "-"]    scx = scx[1:end-1]   end
    	    subXx = parse( scx )
    	    if      ishell == 0  &&  typeof(subXx)     == Void                                       subX = subJ
    	    elseif  ishell >= 1  &&  subshellX[ishell] == AngularJ64(0)  && typeof(subXx) == Void    subX = subJ
    	    elseif  ishell >= 1  &&  typeof(subXx)     == Void                                       subX = subshellX[ishell]
    	    elseif  typeof(subXx) == Void                     subX = subshellX[ishell-1]
    	    elseif  typeof(subXx) == Int64                    subX = AngularJ64(subXx)
    	    elseif  typeof(subXx) == Expr                     subX = AngularJ64( subXx.args[2]// subXx.args[3] )
    	    else    error("stop c")
    	    end
    
    	    # Fill the arrays due to the pre-defined sequence of subshells
    	    while  true
    		ishell = ishell + 1
    		if  subshells[ishell] == sh
    		    wa = ManyElectron.provideSubshellStates(sh, occ);    nu = -1
    		    for  a in wa
    			if   AngularJ64( a.Jsub2//2 ) == subJ	 nu = a.nu;    break   end
    		    end
    		    nu < 0    &&    error("stop d")
    		    push!(occupation, occ);    push!(seniority, nu);	 push!(subshellJ, subJ);    push!(subshellX, subX) 
    		    break
    		else
    		    push!(occupation, 0);      push!(seniority, 0);	 push!(subshellJ, AngularJ64(0));    push!(subshellX, subshellX[ishell-1]) 
    		end
    	    end
    
    	    ishell > length(subshells)    &&	error("stop e")
    
    	    # Terminate if all subshells are read
    	    if  9i + 10 >  length(sa)	 break    end
    	end

      while  ishell < length(subshells)
    	    ishell = ishell + 1
          push!(occupation, 0);    push!(seniority, 0);	 push!(subshellJ, AngularJ64(0));    push!(subshellX, subshellX[ishell-1]) 
      end

    	# Fill the remaining subshells

    	J = subshellX[end];    scx = strip(sc)
    	if	  scx[end:end] == "+"  parity = Parity("+")
    	elseif  scx[end:end] == "-"  parity = Parity("-")
    	else	  error("stop f")
    	end
    
    	wa = CsfR(true, J, parity, occupation, seniority, subshellJ, subshellX, subshellsx)
    end
    

    """
    `ManyElectron.CsfRTransformToStandardSubshells(csf::CsfR, subshells::Array{Subshell,1})`  
        ... constructor to re-define a CSF with regard to the given list of subshells. It incorporates 'empty' subshell state 
            and their coupling (if necessary) but does not support a re-arrangement of the subshell order. This constructor 
            is used to 'combine' individual CsfR into a common list (such as in a Basis) and, thus, use StandardSubshells = true 
            is set. This method terminates with an error message if the given subshell list requires a re-coupling of the CSF.
    """
    function CsfRTransformToStandardSubshells(csf::CsfR, subshells::Array{Subshell,1})
    	J = csf.J;    parity = csf.parity; 
    	occupation = Int64[];	 seniority = Int64[];	  subshellJ = AngularJ64[];	subshellX = AngularJ64[]
    	
    	ish = 0
    	for  subsh in subshells
    	    push!(subshells, subsh)
    	    found = false
    	    for  i in 1:length(csf.subshells)
    		if     subsh == csf.subshells[i]  &&   i >   ish    i = ish;  found = true;   break    
    		elseif subsh == csf.subshells[i]  &&   i <=  ish    error("CSF requires recoupling; this is not supported.")
    		end
    	    end
    	    if  found
    		push!(occupation, csf.occupation[ish]);    push!(seniority, csf.seniority[ish])   
    		push!(subshellJ,  csf.subshellJ[ish]);     push!(subshellX, csf.subshellX[ish])   
    	    else
    		push!(occupation, 0);	 push!(seniority, 0);	 push!(subshellJ,  0)
    		if    ish == 0    push!(subshellX, AngularJ64(0) )   
    		else		  push!(subshellX, subshellX[end])  
    		end 
    	    end
    	end
    
    	CsfR(true, J, parity, occupation, seniority, subshellJ, subshellX, Subshell[] )
    end
    
    
    """
    `struct  ManyElectron.Basis`  
        ... defines a type for a relativistic atomic basis, including the (full) specification of the configuration space 
            and radial orbitals.

        + isDefined      ::Bool                      ... Determines whether this atomic basis 'points' to some a well-defined basis
                                                         or to an empty instance, otherwise.
        + NoElectrons    ::Int64                     ... No. of electrons.
        + subshells      ::Array{Subshell,1}         ... Explicitly given subshell list for this basis.
        + csfs           ::Array{CsfR,1}             ... List of CSF.	    
        + coreSubshells  ::Array{Subshell,1}         ... List of (one-electron) core orbitals.
        + orbitals       ::Dict{Subshell, Orbital}   ... Dictionary of (one-electron) orbitals.
    """
    struct  Basis
        isDefined	       ::Bool
    	  NoElectrons	 ::Int64 
        subshells	       ::Array{Subshell,1}
        csfs             ::Array{CsfR,1}	      
    	  coreSubshells	 ::Array{Subshell,1}
    	  orbitals	       ::Dict{Subshell, Orbital}
    end 
    
    
    """
    `ManyElectron.Basis()`  ... constructor for an 'empty' instance::Basis with isDefined == false.
    """
    function Basis()
    	  Basis(false, 0, Subshell[], CsfR[], Subshell[], Dict{Subshell, Orbital}() )
    end
    
    
    """
    `ManyElectron.Basis("from Grasp2013", cslFilename::String, rwfFilename::String)`  
        ... to construct an instance::Basis from the Grasp92/Grasp2013 .csl and .rwf files.
    """
    function Basis(sa::String, cslFilename::String, rwfFilename::String)
        if  sa == "from Grasp2013"
           isDefined = true 
           #
           println("ManyElectron.Basis-aa: warning ... The standard grid is set; make an interactive request if not appropriate.") 
           grid     = Radial.Grid(true)
           basis    = Basics.read("CSF list: Grasp92", cslFilename)
           orbitals = Basics.readOrbitalFileGrasp92(rwfFilename, grid)
        else  error("stop a")
        end

    	  Basis( isDefined, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals) 
    end
    
    
    """
    `Base.show(io::IO, basis::Basis; full=false)`  ... prepares a proper printout of the variable basis::Basis.
    """
    function Base.show(io::IO, basis::Basis) 
        print(io, string(basis), "\n")
    	  if  basis.isDefined
            println(io, "isDefined:          $(basis.isDefined)  ")
            println(io, "NoElectrons:        $(basis.NoElectrons)  ")
            println(io, "subshells:          $(basis.subshells)  ")
            println(io, "csfs:               $(basis.csfs)  ")
            println(io, "coreSubshells:      $(basis.coreSubshells)  ")
            println(io, "orbitals:           $(basis.orbitals)  ")
    	  end 
    end
    
   
    # `Base.string(basis::Basis)`  ... provides a String notation for variable basis::Basis.
    function Base.string(basis::Basis) 
    	if  basis.isDefined
    	    sa = "Atomic basis:  $(length(basis.csfs)) CSF defined for $(basis.NoElectrons) elecrons."
    	else
    	    sa = "Atomic basis:  Undefined."
    	end
    	return(sa)
    end
    

    """
    `struct  ManyElectron.Level`  
        ... defines a type for an atomic level in terms of its quantum number, energy and with regard to a relativistic basis.
    
        + J            ::AngularJ64       ... Total angular momentum J.
        + M            ::AngularM64       ... Total projection M, only defined if a particular sublevel is referred to.
        + parity       ::Parity           ... Parity of the level.
        + index        ::Int64            ... Index of this level in its original multiplet, or 0.
        + energy       ::Float64          ... energy
        + relativeOcc  ::Float64          ... Relative occupation of this level, if involved in the evolution of a cascade.
        + hasStateRep  ::Bool             ... Determines whether this level 'point' to a physical representation of an state
                                              (i.e. an basis and corresponding mixing coefficient vector), or just to empty instances, otherwise.
        + basis        ::Basis            ... relativistic atomic basis
        + mc           ::Vector{Float64}  ... Vector of mixing coefficients w.r.t basis.
    """
    struct  Level
        J		     ::AngularJ64
        M		     ::AngularM64
    	  parity  	     ::Parity
    	  index  	     ::Int64
    	  energy  	     ::Float64
    	  relativeOcc    ::Float64
    	  hasStateRep    ::Bool
    	  basis	     ::Basis
    	  mc		     ::Vector{Float64}
    end 
    
    
    """
    `ManyElectron.Level()`  ... constructor an empty level without useful data.
    """
    function Level()
        Level( AngularJ64(0),  AngularM64(0), Parity("+"), 0, 0., 0., false, Basis(), Vector{Float64}[] )
    end
    
    
    """
    `ManyElectron.Level(J::AngularJ64, parity::Parity, energy::Float64, relativeOcc::Float64)`  
        ... constructor an atomic level without a state representation: hasStateRep == false.
    """
    function Level(J::AngularJ64, parity::Parity, energy::Float64, relativeOcc::Float64)
        Level(J, J, parity, 0, energy, relativeOcc, false, Basis(), Vector{Float64}[] )
    end
    
    
    # `Base.show(io::IO, level::Level)`  ... prepares a proper printout of the variable level::Level.
    function Base.show(io::IO, level::Level) 
        println(io, "Level: J = $(level.J), M = $(level.M), parity = $(level.parity), index = $(level.index) ")
        println(io, "energy:         $(level.energy)  ")
        println(io, "relativeOcc:    $(level.relativeOcc)  ")
        println(io, "hasStateRep:    $(level.hasStateRep)  ")
        println(io, "basis:           (level.basis)  ")
        println(io, "mc:             $(level.mc)  ")
    end
    
    
    """
    `struct  ManyElectron.Multiplet`  
        ... defines a type for an ordered list of atomic levels; has only the default constructor.
    
    	  + name  	 ::String	         ... A name associated to the multiplet.
    	  + levels	 ::Array{Level,1}   ... A list of levels (pointers).
    """
    struct  Multiplet
        name       ::String
    	  levels  	 ::Array{Level,1}
    end 


    """
     `ManyElectron.Multiplet()`  ... constructor for providing an 'empty' instance of this struct.
    """
    function Multiplet()
        Multiplet("", Level[] )
    end
    
    
    """
    `ManyElectron.Multiplet("from Ratip2012", cslFilename::String, rwfFilename::String, mixFilename::String)`  
        ... to construct an instance::Multiplet from the Grasp92/Ratip2012 .csl, .rwf and .mix files.  Some consistency 
            checks are made and the method terminates with an error message if the files don't fit together. 
        
    + `("from Grasp2018", cslFilename::String, rwfFilename::String, mixFilename::String)` 
        ... to do the same but for files from Grasp2018.
    """
    function Multiplet(sa::String, cslFilename::String, rwfFilename::String, mixFilename::String)
        if  sa == "from Ratip2012"
           name      = "from "*cslFilename[1:end-4] 
           basis     = ManyElectron.Basis("from Grasp2013", cslFilename, rwfFilename)
           multiplet = Basics.readMixFileRelci(mixFilename, basis)
        else  error("Unrecognized keystrings:  $sa")
        end

        return( multiplet ) 
    end
    
    
    # `Base.show(io::IO, multiplet::Multiplet)`  ... prepares a proper printout of the variable multiplet::Multiplet.
    function Base.show(io::IO, multiplet::Multiplet) 
        println(io, "name:        $(multiplet.name)  ")
        println(io, "levels:      $(multiplet.levels)  ")
    end
    
    
    """
    `struct  ManyElectron.MultipletSettings`  
        ... a struct for defining the atomic interactions to be incorporated into the representation of a multiplet.
    
    	+ Coulomb		           ::Bool 		      ... logical flag to include Coulomb interactions.
    	+ eeInteraction             ::AbstractEeInteraction 		      ... logical flag to include Breit interactions.
    	+ QED			           ::Bool 		      ... logical flag to include QED interactions.
    	+ diagonalizationMethod    ::String		      ... method for diagonalizing the matrix.
    	+ selectLevels  	       ::Bool 		      ... true, if specific level (number)s have been selected for computation.
    	+ selectedLevels	       ::Array{Int64,1}	  ... Level number that have been selected.
    	+ selectSymmetries	       ::Bool 		      ... true, if specific level symmetries have been selected for computation.
    	+ selectededSymmetries     ::Array{LevelSymmetry,1}   ... Level symmetries that have been selected.
    
    	### `Methods for diagonalization`  
    
    	* `eigval`	       The internal Julia method for diagonalizing a quadratic matrix
    	*
    """
    struct  MultipletSettings
        Coulomb 		         ::Bool
        eeInteraction             ::AbstractEeInteraction
        QED			             ::Bool
    	diagonalizationMethod	 ::String
        selectLevels		     ::Bool
        selectedLevels  	     ::Array{Int64,1}
        selectSymmetries	     ::Bool
        selectedSymmetries	     ::Array{LevelSymmetry,1}
    end
    
    
    # Functions/methods that are later added to the module ManyElectron
    function provideSubshellStates                                  end
    
end # module
    
    
