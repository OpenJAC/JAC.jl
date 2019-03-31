
"""
`module  JAC.ManyElectron`  ... a submodel of JAC that contains all structs and methods for many-electron CSF, basis and levels; 
                                it is using JAC, JAC.Radial.
"""
module ManyElectron

    using  JAC, JAC.Radial
    export Configuration, ConfigurationR, AsfSettings, CsfR, Basis, Level, Multiplet


    """
    `struct  Configuration`  ... defines a struct for a non-relativistic electron configuration that is fully speficied by its shell notations 
                                 (such as "1s", "2s", "2p", ....) and the corresponding occupation numbers (>= 0). An electron configuration 
             is independent of any order of the shells, and a zero occupation is assumed for all shells that do not appear as a key in the shell 
             (dictionary). Therefore, the number of keys do not allow any conclusion about the underlying orbital space of any considered 
             computation that include more than just a single configuration.

        + shells         ::Dict{Shell,Int64}   ... Dictionary that maps shells to their occupation.
        + NoElectrons    ::Int64               ... No. of electrons.
    """
    struct  Configuration
        shells           ::Dict{Shell,Int64}  
        NoElectrons      ::Int64               
    end


    """
    `JAC.Configuration(sa::String)`  ... constructor for a given configuration string, such as "[He]", "[Ne]", "[Ne] 3s 3p^6"  or 
                                         "1s 2p^6 3s^2 3p".
    """
    function Configuration(sa::String)
        sax = strip(sa)
        shellOccList = split(sax, " ")
        println("shellOccList = $shellOccList")
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


    """
    `Base.show(io::IO, conf::Configuration)`  ... prepares a proper printout of the non-relativistic configuration conf::Configuration.
    """
    function Base.show(io::IO, conf::Configuration)
        sa = "Configuration: " * string(conf)
        print(io, sa)
    end


    """
    `Base.string(conf::Configuration)`  ... provides a String notation for the variable conf::Configuration.
    """
    function Base.string(conf::Configuration)
        wa = keys(conf.shells);   va = values(conf.shells)
        wb = JAC.give("ordered shell list: non-relativistic", 7)

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
    `Base.:(==)(confa::Configuration, confb::Configuration)`  ... compares (recursively) two non-relativistic configurations and return true
                                                                  if all subfields are equal, and false otherwise
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
    `struct  ConfigurationR`  ... defines a type for a relativistic electron configuration that is fully speficied by its subshell notations 
                                  (such as "1s_1/2", "2s_1/2", "2p_3/2", ....) and the corresponding occupation numbers (>= 0). An 
             electron configuration is independent of any order of the subshells, and a zero occupation is assumed for all shells that do not 
             appears as a key. Therefore, the number of keys do not allow any conclusion about the underlying orbital space of any considered 
             computation that include more than a single configuration.

        + subshells      ::Dict{Subshell,Int64}  ... Dictionary that maps subshells to their occupation.
        + NoElectrons    ::Int64                 ... No. of electrons.
    """
    struct  ConfigurationR
        subshells        ::Dict{Subshell,Int64}  
        NoElectrons      ::Int64               
    end


    """
    `Base.show(io::IO, conf::ConfigurationR)`  ... prepares a proper printout of the relativistic configuration conf::ConfigurationR.
    """
    function Base.show(io::IO, conf::ConfigurationR)
        sa = "Configuration: " * string(conf)
        print(io, sa)
    end


    """
    `Base.string(conf::ConfigurationR)`  ... provides a String notation for the variable conf::ConfigurationR.
    """
    function Base.string(conf::ConfigurationR)
        wa = keys(conf.subshells);   va = values(conf.subshells)
        wb = JAC.give("ordered subshell list: relativistic", 7)

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
    `Base.:(==)(confa::ConfigurationR, confb::ConfigurationR)`  ... compares (recursively) two relativistic configurations and return true
                                                                    if all subfields are equal, and false otherwise
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
    `struct  AsfSettings`  ... a struct for defining the settings for the atomic state functions, i.e. the self-consistent-field (SCF) 
                               and CI computations

        + generateScf          ::Bool               ... True, if a SCF need to be generated, and false otherwise (frozen orbitals).
        + breitScf             ::Bool               ... True, if Breit interaction is to be included into the SCF computations.
        + methodScf            ::String             ... Specify the SCF method: ["AL", "OL", "EOL", "meanDFS", "meanHS"].
        + startScf             ::String             ... Specify how the start orbitals are obtained ["fromNRorbitals", "fromGrasp", "hydrogenic"].
        + orbitalFileScf       ::String             ... Filename of orbitals, if taken from Grasp.
        + levelsScf            ::Array{Int64,1}     ... Levels on which the optimization need to be carried out.
        + maxIterationsScf     ::Int64              ... maximum number of SCF iterations
        + áccuracyScf          ::Float64            ... convergence criterion for the SCF field.
        + shellSequenceScf     ::Array{Subshell,1}  ... Sequence of subshells to be optimized.
        
    	+ coulombCI            ::Bool               ... logical flag to include Coulomb interactions.
    	+ breitCI              ::Bool               ... logical flag to include Breit interactions.
    	+ qedCI                ::Bool               ... logical flag to include QED interactions.
    	+ methodCI             ::String             ... method for diagonalizing the matrix.
    	+ selectLevelsCI       ::Bool               ... true, if specific level (number)s have been selected for computation.
    	+ selectedLevelsCI     ::Array{Int64,1}     ... Level number that have been selected.
    	+ selectSymmetriesCI   ::Bool               ... true, if specific level symmetries have been selected for computation.
    	+ selectedSymmetriesCI ::Array{LevelSymmetry,1}       ... Level symmetries that have been selected.
    """
    struct  AsfSettings
        generateScf            ::Bool 
        breitScf               ::Bool  
        methodScf              ::String  
        startScf               ::String 
        orbitalFileScf         ::String  
        levelsScf              ::Array{Int64,1}
        maxIterationsScf       ::Int64 
        accuracyScf            ::Float64   
        shellSequenceScf       ::Array{Subshell,1}
        #
    	coulombCI		       ::Bool 
    	breitCI		           ::Bool 
    	qedCI			       ::Bool 	
    	methodCI               ::String	
    	selectLevelsCI         ::Bool 
    	selectedLevelsCI       ::Array{Int64,1}
    	selectSymmetriesCI     ::Bool 	
    	selectedSymmetriesCI   ::Array{LevelSymmetry,1} 
     end


    """
    `JAC.AsfSettings()`  ... constructor for setting the default values.
    """
    function AsfSettings()
    	AsfSettings(false, false, "meanDFS", "hydrogenic", "", Int64[1], 24, 1.0e-6, Subshell[],   
    	            true, false, false, "eigen", false, Int64[], false, LevelSymmetry[] )
    end
    
    
    """
    `Base.show(io::IO, settings::AsfSettings)`  ... prepares a proper printout of the settings::AsfSettings.
    """
    function Base.show(io::IO, settings::AsfSettings)
    	  println(io, "generateScf:          $(settings.generateScf)  ")
    	  println(io, "breitScf:             $(settings.breitScf)  ")
    	  println(io, "methodScf:            $(settings.methodScf)  ")
    	  println(io, "startScf:             $(settings.startScf)  ")
    	  println(io, "orbitalFileScf:       $(settings.orbitalFileScf)  ")
    	  println(io, "levelsScf:            $(settings.levelsScf)  ")
    	  println(io, "maxIterationsScf:     $(settings.maxIterationsScf)  ")
    	  println(io, "accuracyScf:          $(settings.accuracyScf)  ")
    	  println(io, "shellSequenceScf:     $(settings.shellSequenceScf)  ")
    	  #
    	  println(io, "coulombCI:            $(settings.coulombCI)  ")
    	  println(io, "breitCI:              $(settings.breitCI)  ")
    	  println(io, "qedCI:                $(settings.qedCI)  ")
    	  println(io, "methodCI:             $(settings.methodCI)  ")
    	  println(io, "selectLevelsCI:       $(settings.selectLevelsCI)  ")
    	  println(io, "selectedLevelsCI:     $(settings.selectedLevelsCI)  ")
    	  println(io, "selectSymmetriesCI:   $(settings.selectSymmetriesCI)  ")
    	  println(io, "selectedSymmetriesCI: $(settings.selectedSymmetriesCI)  ")

    end
    
    
    """
    `Base.string(settings::AsfSettings)`  ... provides a String notation for the variable settings::AsfSettings.
    """
    function Base.string(settings::AsfSettings)
    	  error("Not yet implemented.")
    	  sa = "Asf settings: maximum No. of iterations = $(settings.maxIterationsScf), accuracy = (settings.accuracyScf)"
    	  return( sa )
    end
    
    
#=
    """
    `struct  ScfSettings`  ... a struct for defining the settings for self-consistent-field (SCF) computations

        + generateScf       ::Bool               ... True, if a SCF need to be generated, and false if just the start orbitals should be applied.
        + startOrbitals     ::String             ... Specify how the start orbitals are obtained ["fromNRorbitals", "fromGrasp", "hydrogenic"].
        + rwfFilename       ::String             ... Filename of orbitals, if taken from Grasp.
        + levels            ::Array{Int64,1}     ... Levels on which the optimization need to be carried out.
        + includeBreit      ::Bool               ... True, if the Breit interaction is included into the SCF generation.
        + maxIterations     ::Int64              ... maximum number of SCF iterations
        + scfAccuracy       ::Float64            ... convergence criterion for the SCF field.
        + iterationSequence ::Array{Subshell,1}  ... Sequence of subshells to be optimized.
    """
    struct  ScfSettings
        generateScf         ::Bool
        startOrbitals       ::String
        rwfFilename         ::String
        levels              ::Array{Int64,1}
        includeBreit        ::Bool
        maxIterations       ::Int64
        scfAccuracy         ::Float64
        iterationSequence   ::Array{Subshell,1}
     end


    """
    `JAC.ScfSettings()`  ... constructor for setting the default values.
    """
    function ScfSettings()
    	ScfSettings(false, "hydrogenic", "", Int64[1], false, 24, 1.0e-8, Subshell[] )
    end
    
    
    """
    `Base.show(io::IO, settings::ScfSettings)`  ... prepares a proper printout of the settings::ScfSettings.
    """
    function Base.show(io::IO, settings::ScfSettings)
    	  println(io, "generateScf:	 $(settings.generateScf)  ")
    	  println(io, "startOrbitals:	 $(settings.startOrbitals)  ")
    	  println(io, "rwfFilename:	 $(settings.rwfFilename)  ")
    	  println(io, "levels:		 $(settings.levels)  ")
    	  println(io, "includeBreit:	 $(settings.includeBreit)  ")
    	  println(io, "maxIterations:	 $(settings.maxIterations)  ")
    	  println(io, "scfAccuracy:	 $(settings.scfAccuracy)  ")
    	  println(io, "iterationSequence:  $(settings.iterationSequence)  ")
    end
    
    
    """
    `Base.string(settings::ScfSettings)`  ... provides a String notation for the variable settings::ScfSettings.
    """
    function Base.string(settings::ScfSettings)
    	  error("Not yet implemented.")
    	  sa = "SCF settings: maximum No. of iterations = $(settings.maxIterations), accuracy = (settings.scfAccuracy)"
    	  return( sa )
    end                                                                      =#
    
    
    """
    `struct  CsfR`  ... defines a type for a relativistic configuration state function (CSF) in terms of its orbitals (sequence of orbitals), 
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
    `JAC.CsfR(J::AngularJ64, parity::Parity)`  ... simple constructor for given J and parity, and where standardOrbitals is set to false.
    """
    function CsfR(J::AngularJ64, parity::String)
    	  CsfR(false, J, parity, Int64[], Int64[], AngularJ64[], AngularJ64[], Subshell[] )
    end
    
    
    """
    `JAC.CsfR(sa::String)`  ... simple constructor for given closed-shell occupation (such as "[Ne]"), and where standardOrbitals is set to true.
    """
    function CsfR(sa::String)
    	  wa = subshellsFromClosedShellConfiguration(sa)
    
    	  occupation = Int64[];	  subshellJ = AngularJ64[];	subshellX  = AngularJ64[]
    	  subshells  = Subshell[];  seniority = Int64[]
    	
        for  sb  in  wa
    	      #x wb = JAC.SubshellQuantumNumbers(string(sb));   kappa = wb[2];   occ = wb[4] + 1
    	      occ = subshell_2j(sb) + 1
    	      occupation = push!(occupation, occ)
    	      seniority  = push!(seniority,    0) 
    	      subshellJ  = push!(subshellJ, AngularJ64(0) ) 
    	      subshellX  = push!(subshellX, AngularJ64(0) ) 
    	      subshells  = push!(subshells, sb ) 
    	  end
    	  J = AngularJ64(0);    parity = plus
    
    	  CsfR(false, J, parity, occupation, seniority, subshellJ, subshellX, subshells)
    end
    
    
    """
    `Base.show(io::IO, csf::CsfR)`  ... prepares a proper printout of the variable orbital::CsfR.
    """
    function Base.show(io::IO, csf::CsfR) 
    	  print(io, string(csf) ) 
    end
    
    
    """
    `Base.string(csf::CsfR)`  ... provides a String notation for csf::CsfR.
    """
    function Base.string(csf::CsfR) 
    	  sa = "\n   CSF: "
        if    csf.useStandardSubshells

    	      for  i in 1:length(csf.occupation)
    		    sa = sa * JAC.subshellStateString(string(JAC.JAC_STANDARD_SUBSHELL_LIST[i]), csf.occupation[i], csf.seniority[i], 
                                                  csf.subshellJ[i], csf.subshellX[i])
    		    sa = sa * ", "
    	      end
    	      sa = sa * ": J=" * string(csf.J) * string(csf.parity)

            else
    	      for  i in 1:length(csf.subshells)
    		    sa = sa * JAC.subshellStateString(string(csf.subshells[i]), csf.occupation[i], csf.seniority[i], 
                                                  csf.subshellJ[i], csf.subshellX[i])
    		    sa = sa * ", "
    	      end
    	      sa = sa * ": J=" * string(csf.J) * string(csf.parity)
    	  end
    	  return( sa )
    end
    

    """
    `Base.:(==)(csfa::CsfR, csfb::CsfR)`  ... compares (recursively) two relativistic CSFs and return true if all subfields are equal, 
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
    `JAC.CsfRGrasp92(subshells::Array{Subshell,1}, coreSubshells::Array{Subshell,1}, sa::String, sb::String, sc::String)`  ... to construct a CsfR 
    	   from 3 Grasp-Strings and the given list of subshells and core subshells.
    """
    function CsfRGrasp92(subshells::Array{Subshell,1}, coreSubshells::Array{Subshell,1}, sa::String, sb::String, sc::String)
    	occupation = Int64[];	 seniority = Int64[];	 subshellJ = AngularJ64[];    subshellX = AngularJ64[];    subshellsx = Subshell[]

    	# Define the occupation and quantum numbers of the closed shells
    	for i in 1:length(coreSubshells)
    	    push!(occupation, JAC.subshell_2j( coreSubshells[i] ) +1);    
    	    push!(seniority, 0);     push!(subshellJ, AngularJ64(0));	 push!(subshellX, AngularJ64(0)) 
    	end
    
    	scc = sc * "	  ";   scx = " "
    	i  = -1;   ishell = length(coreSubshells)
    	while true
    	    i = i + 1;    sax = sa[9i+1:9i+9];     sbx = sb[9i+1:9i+9];    scx = strip( scc[9i+6:9i+14] )
    
    	    sh = strip( sax[1:5] );    sh = JAC.subshellGrasp(sh);    occ = parse( sax[7:8] )
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
    		    wa = JAC.provide("subshell states: antisymmetric, seniority", sh, occ);    nu = -1
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
    `JAC.CsfRTransformToStandardSubshells(csf::CsfR, subshells::Array{Subshell,1})`  ... constructor to re-define a CSF with regard to the 
    	   given list of subshells. It incorporates 'empty' subshell state and their coupling (if necessary) but does not support a re-arrangement 
         of the subshell order. This constructor is used to 'combine' individual CsfR into a common list (such as in a Basis) and, thus, use
         StandardSubshells = true is set. This method terminates with an error message if the given subshell list requires a re-coupling 
         of the CSF.
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
    `struct  Basis`  ... defines a type for a relativistic atomic basis, including the (full) specification of the 
    			  configuration space and radial orbitals.

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
    `JAC.ManyElectron.Basis()`  ... constructor for an 'empty' instance::Basis with isDefined == false.
    """
    function Basis()
    	  Basis(false, 0, Subshell[], CsfR[], Subshell[], Dict{Subshell, Orbital}() )
    end
    
    
    """
    `JAC.ManyElectron.Basis("from Grasp2013", cslFilename::String, rwfFilename::String)`  ... to construct an instance::Basis from 
                             the Grasp92/Grasp2013 .csl and .rwf files.
    """
    function Basis(sa::String, cslFilename::String, rwfFilename::String)
        if  sa == "from Grasp2013"
           isDefined = true 
           #
           println("JAC.ManyElectron.Basis-aa: warning ... The standard grid is set; make an interactive request if not appropriate.") 
           grid     = JAC.Radial.Grid("grid: exponential")
           basis    = JAC.read("CSF list: Grasp92", cslFilename)
           orbitals = JAC.readOrbitalFileGrasp92(rwfFilename, grid)
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
    
   
    """
    `Base.string(basis::Basis)`  ... provides a String notation for variable basis::Basis.
    """
    function Base.string(basis::Basis) 
    	if  basis.isDefined
    	    sa = "Atomic basis:  $(length(basis.csfs)) CSF defined for $(basis.NoElectrons) elecrons."
    	else
    	    sa = "Atomic basis:  Undefined."
    	end
    	return(sa)
    end
    

    """
    `struct  Level`  ... defines a type for an atomic level in terms of its quantum number, energy and with regard to a relativistic basis.
    
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
    `JAC.Level()`  ... constructor an empty level without useful data.
    """
    function Level()
        Level( AngularJ64(0),  AngularM64(0), Parity("+"), 0, 0., 0., false, Basis(), Vector{Float64}[] )
    end
    
    
    """
    `JAC.Level(J::AngularJ64, parity::Parity, energy::Float64, relativeOcc::Float64)`  
        ... constructor an atomic level without a state representation: hasStateRep == false.
    """
    function Level(J::AngularJ64, parity::Parity, energy::Float64, relativeOcc::Float64)
        Level(J, J, parity, 0, energy, relativeOcc, false, Basis(), Vector{Float64}[] )
    end
    
    
    """
    `Base.show(io::IO, level::Level)`  ... prepares a proper printout of the variable level::Level.
    """
    function Base.show(io::IO, level::Level) 
        println(io, "Level: J = $(level.J), M = $(level.M), parity = $(level.parity), index = $(level.index) ")
        println(io, "energy:         $(level.energy)  ")
        println(io, "relativeOcc:    $(level.relativeOcc)  ")
        println(io, "hasStateRep:    $(level.hasStateRep)  ")
        println(io, "basis:           (level.basis)  ")
        println(io, "mc:             $(level.mc)  ")
    end
    
    
    """
    `struct  Multiplet`  ... defines a type for an ordered list of atomic levels; has only the default constructor.
    
    	  + name  	 ::String	         ... A name associated to the multiplet.
    	  + levels	 ::Array{Level,1}   ... A list of levels (pointers).
    """
    struct  Multiplet
        name       ::String
    	  levels  	 ::Array{Level,1}
    end 


    """
     `JAC.ManyElectron.Multiplet()`  ... constructor for providing an 'empty' instance of this struct.
    """
    function Multiplet()
        Multiplet("", Level[] )
    end
    
    
    """
    `JAC.ManyElectron.Multiplet("from Ratip2012", cslFilename::String, rwfFilename::String, mixFilename::String)`  ... to construct an 
         instance::Multiplet from the Grasp92/Ratip2012 .csl, .rwf and .mix files.  Some consistency checks are made and the method terminates 
         with an error message if the files don't fit together. 
        
        + `("from Grasp2018", cslFilename::String, rwfFilename::String, mixFilename::String)` ... to do the same but for files from Grasp2018.
    """
    function Multiplet(sa::String, cslFilename::String, rwfFilename::String, mixFilename::String)
        if  sa == "from Ratip2012"
           name      = "from "*cslFilename[1:end-4] 
           basis     = JAC.ManyElectron.Basis("from Grasp2013", cslFilename, rwfFilename)
           multiplet = JAC.readMixFileRelci(mixFilename, basis)
        else  error("Unrecognized keystrings:  $sa")
        end

        return( multiplet ) 
    end
    
    
    """
    `Base.show(io::IO, multiplet::Multiplet)`  ... prepares a proper printout of the variable multiplet::Multiplet.
    """
    function Base.show(io::IO, multiplet::Multiplet) 
        println(io, "name:        $(multiplet.name)  ")
        println(io, "levels:      $(multiplet.levels)  ")
    end
    
    
    """
    `struct  MultipletSettings`  ... a struct for defining the atomic interactions to be incorporated into the representation of a multiplet.
    
    	+ Coulomb		       ::Bool 		          ... logical flag to include Coulomb interactions.
    	+ Breit 		       ::Bool 		          ... logical flag to include Breit interactions.
    	+ QED			       ::Bool 		          ... logical flag to include QED interactions.
    	+ diagonalizationMethod  ::String		          ... method for diagonalizing the matrix.
    	+ selectLevels  	       ::Bool 		          ... true, if specific level (number)s have been selected for computation.
    	+ selectedLevels	       ::Array{Int64,1}	          ... Level number that have been selected.
    	+ selectSymmetries	 ::Bool 		          ... true, if specific level symmetries have been selected for computation.
    	+ selectededSymmetries   ::Array{LevelSymmetry,1}   ... Level symmetries that have been selected.
    
    	### `Methods for diagonalization`  
    
    	* `eigval`	       The internal Julia method for diagonalizing a quadratic matrix
    	*
    """
    struct  MultipletSettings
    	  Coulomb 		       ::Bool
        Breit			 ::Bool
    	  QED			       ::Bool
    	  diagonalizationMethod	 ::String
        selectLevels		 ::Bool
    	  selectedLevels  	 ::Array{Int64,1}
        selectSymmetries	 ::Bool
    	  selectedSymmetries	 ::Array{LevelSymmetry,1}
    end

#=    
    """
    `JAC.MultipletSettings()`  ... constructor for a standard instance, ie. by just including Coulomb interactions.
    """
    function MultipletSettings()
    	  MultipletSettings(true, false, false, "eigval", false, Int64[], false, LevelSymmetry[] )
    end

    
    """
    `JAC.MultipletSettings(coulomb::Bool, breit::Bool, qed::Bool)`  ... constructor for a standard instance, ie. by just defining the atomic 
                                                                        interactions.
    """
    function MultipletSettings(coulomb::Bool, breit::Bool, qed::Bool)
    	  MultipletSettings(coulomb, breit, qed, "eigval", false, Int64[], false, LevelSymmetry[] )
    end


    """
    `Base.show(io::IO, settings::MultipletSettings)`  ... prepares a proper printout of the settings::MultipletSettings.
    """
    function Base.show(io::IO, settings::MultipletSettings)
    	println(io, "Coulomb:		       $(settings.Coulomb)  ")
    	println(io, "Breit:		       $(settings.Breit)  ")
    	println(io, "QED:		             $(settings.QED)  ")
    	println(io, "diagonalizationMethod:  $(settings.diagonalizationMethod)  ")
    	println(io, "selectLevels:	       $(settings.selectLevels)  ")
    	println(io, "selectedLevels:	       $(settings.selectedLevels)  ")
    	println(io, "selectSymmetries:       $(settings.selectSymmetries)  ")
    	println(io, "selectedSymmetries:     $(settings.selectedSymmetries)  ")
    end
    
    
    """
    `Base.string(settings::MultipletSettings)`  ... provides a String notation for the variable settings::MultipletSettings.
    """
    function Base.string(settings::MultipletSettings)
    	error("Not yet implemented.")
    	sa = " "
    	if   settings.Coulomb	 sa = sa * "Coulomb  "   end
    	if   settings.Breit	 sa = sa * "Breit  "	 end
    	if   settings.QED	 sa = sa * "QED  "	 end
    	sa = "Atomic interactions (to be) included: " * sa
    	return( sa )
    end                                        =#
    
end # module
    
    
