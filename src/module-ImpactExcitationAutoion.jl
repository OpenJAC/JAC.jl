
"""
`module  JAC.ImpactExcitationAutoion`  ... a submodel of JAC that contains all methods for computing electron impact excitation cross 
                                           sections and collision strengths; it is using JAC, JAC.ManyElectron, JAC.Radial, JAC.ImpactExcitation, 
                                           JAC.AutoIonization.
"""
module ImpactExcitationAutoion 

    using Printf, JAC, JAC.ManyElectron, JAC.Radial, JAC.ImpactExcitation, JAC.AutoIonization
    global JAC_counter = 0


    """
    `struct  ImpactExcitationAutoion.Settings`  ... defines a type for the details and parameters of computing electron-impact 
                                                    excitation-autoionization pathways |i(N)>  --> |m(N)>  --> |f(N-1)>.

        + electronEnergies        ::Array{Float64,1}                   ... List of impact-energies of the incoming elecgtrons.
        + printBeforeComputation  ::Bool                               ... True, if all energies and lines are printed before their evaluation.
        + selectPathways          ::Bool                               ... True if particular pathways are selected for the computations.
        + selectedPathways        ::Array{Tuple{Int64,Int64,Int64},1}  ... List of list of pathways, given by tupels (inital, inmediate, final).
        + maxKappa                ::Int64                              ... Maximum kappa value of partial waves to be included.
    """
    struct Settings
        electronEnergies          ::Array{Float64,1}
        printBeforeComputation    ::Bool
        selectPathways            ::Bool
        selectedPathways          ::Array{Tuple{Int64,Int64,Int64},1}
        maxKappa                  ::Int64 
    end 


    """
    `JAC.ImpactExcitationAutoion.Settings()`  ... constructor for the default values of electron-impact excitation-autoionizaton settings.
    """
    function Settings()
       Settings( Float64[], false, false, Tuple{Int64,Int64,Int64}[], 0)
    end


    # `Base.show(io::IO, settings::ImpactExcitationAutoion.Settings)`  
    #   ... prepares a proper printout of the variable settings::ImpactExcitationAutoion.Settings.  
    function Base.show(io::IO, settings::ImpactExcitationAutoion.Settings) 
        println(io, "electronEnergies:         $(settings.electronEnergies)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectPathways:           $(settings.selectPathways)  ")
        println(io, "selectedPathways:         $(settings.selectedPathways)  ")
        println(io, "maxKappa:                 $(settings.maxKappa)  ")
    end


    """
    `struct  ImpactExcitationAutoion.Channel`  ... defines a type for a electron-impact excitaton & autoionization channel that specifies 
                                                   all quantum numbers, phases and amplitudes.

        + excitationChannel  ::JAC.ImpactExcitation.Channel      ... Channel that describes the electron-impact excitation process.
        + augerChannel       ::JAC.AutoIonization.Channel                 ... Channel that describes the subsequent Auger/autoionization process.
    """
    struct  Channel
        excitationChannel    ::JAC.ImpactExcitation.Channel
        augerChannel         ::JAC.AutoIonization.Channel
    end 


    """
    `struct  ImpactExcitationAutoion.Pathway`  ... defines a type for a electron-impact excitation pathway that may include the definition 
                                                   of different excitation and autoionization channels and their corresponding amplitudes.

        + initialLevel        ::Level       ... initial-(state) level
        + intermediateLevel   ::Level       ... intermediate-(state) level
        + finalLevel          ::Level       ... final-(state) level
        + electronInEnergy    ::Float64     ... energy of the (incoming) electron
        + electronOutEnergy   ::Float64     ... energy of the (outgoing, scattered) electron
        + electronAugerEnergy ::Float64     ... energy of the (emitted Auger) electron
        + crossSection        ::Float64     ... total cross section of this pathway
        + hasChannels         ::Bool        ... Determines whether the individual excitation and autoionization channels are defined in terms of 
                                                their free-electron kappa's, phases and the total angular momentum/parity as well as the 
                                                amplitude, or not.
        + channels            ::Array{ImpactExcitationAutoion.Channel,1}  ... List of channels of this pathway.
    """
    struct  Pathway
        initialLevel          ::Level
        intermediateLevel     ::Level
        finalLevel            ::Level
        electronInEnergy      ::Float64
        electronOutEnergy     ::Float64
        electronAugerEnergy   ::Float64
        crossSection          ::Float64 
        hasChannels           ::Bool
        channels              ::Array{ImpactExcitationAutoion.Channel,1}
    end 


    """
    `JAC.ImpactExcitationAutoion.Pathway()`  ... constructor for an electron-impact excitation-autoionization pathway between a specified 
         initial, intermediate and final level.
    """
    function Line()
        Pathway(Level(), Level(), Level(), 0., 0., 0., 0., false, ImpactExcitationAutoion.Channel[] )
    end


    # `Base.show(io::IO, pathway::ImpactExcitationAutoion.Pathway)`  
    #   ... prepares a proper printout of the variable pathway::ImpactExcitationAutoion.Pathway.
    function Base.show(io::IO, pathway::ImpactExcitationAutoion.Pathway) 
        println(io, "initialLevel:               $(pathway.initialLevel)  ")
        println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
        println(io, "finalLevel:                 $(pathway.finalLevel)  ")
        println(io, "electronInEnergy:           $(pathway.electronInEnergy)  ")
        println(io, "electronOutEnergy:          $(pathway.electronOutEnergy)  ")
        println(io, "electronAugerEnergy:        $(pathway.electronAugerEnergy)  ")
        println(io, "crossSection:               $(pathway.crossSection)  ")
        println(io, "hasChannels:                $(pathway.hasChannels)  ")
        println(io, "channels:                   $(pathway.channels)  ")
    end


    """
    `JAC.ImpactExcitationAutoion.computeAmplitudesProperties(pathway::ImpactExcitationAutoion.Pathway, grid::Radial.Grid, 
                                                             settings::ImpactExcitationAutoion.Settings)`  ... to compute all amplitudes and 
         properties of the given line; a line::ImpactExcitationAutoion.Line is returned for which the amplitudes and properties have now been 
         evaluated.
    """
    function  computeAmplitudesProperties(pathway::ImpactExcitationAutoion.Pathway, grid::Radial.Grid, settings::ImpactExcitationAutoion.Settings)
        global JAC_counter
        newChannels = ImpactExcitationAutoion.Channel[]
        for channel in pathway.channels
            # Generate a continuum orbital
            JAC_counter = JAC_counter + 1
            if   JAC_counter < 20   println("ImpactExcitationAutoion.computeAmplitudesProperties-aa: warning ... no cont. orb. is generated.")  end
            initialPhase  = 0.;    finalPhase = 0;    augerPhase = 0
            # Define a proper continuum basis from the finalLevel.basis and the continuum orbital
            JAC_counter = JAC_counter + 1
            if   JAC_counter < 20   println("ImpactExcitationAutoion.computeAmplitudesProperties-ab: warning ... no cont. basis generated.") end
            # Compute the transition matrix for the continuum and the initial-state basis
            JAC_counter = JAC_counter + 1
            if   JAC_counter < 20   println("ImpactExcitationAutoion.computeAmplitudesProperties-ac: warning ... no trans-matrix computed.") end
            # matrix    = JAC.PhotoExcitationAutoion.computeMatrix(channel.multipole, channel.gauge, line.omega, line.finalLevel.basis, 
            #                                                      line.initialLevel.basis, grid, settings)
            # amplitude = line.finalLevel.mc * matrix * line.initialLevel.mc 
            amplitude = 1.0 
            exChannel = channel.excitationChannel
            eChannel  = ImpactExcitation.Channel( exChannel.initialKappa, exChannel.finalKappa, exChannel.symmetry, 
                                                  initialPhase, finalPhase, Complex(0.))
            aChannel  = AutoIonization.Channel( channel.augerChannel.kappa, channel.augerChannel.symmetry, augerPhase, Complex(0.))
            push!( newChannels, ImpactExcitationAutoion.Channel(eChannel, aChannel) )
        end
        # Calculate the totalRate 
        JAC_counter = JAC_counter + 1
        if   JAC_counter < 20   println("ImpactExcitationAutoion.computeAmplitudesProperties-ba: warning ... crossSection set to -1.") end
        crossSection = 1.
        pathway = ImpactExcitationAutoion.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, 
                                                   pathway.electronInEnergy, pathway.electronOutEnergy, pathway.electronAugerEnergy, 
                                                   crossSection, true, newChannels)
        return( pathway )
    end


    """
    `JAC.ImpactExcitationAutoion.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                 grid::Radial.Grid, settings::ImpactExcitationAutoion.Settings; output=true)`  ... to compute the 
         electron-impact-excitation-autoionization amplitudes and all properties as requested by the given settings. A list of 
         lines::Array{ImpactExcitationAutoion.Lines} is returned.
    """
    function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                              settings::ImpactExcitationAutoion.Settings; output=true)
        println("")
        printstyled("JAC.ImpactExcitationAutoion.computePathways(): The computation of e-impact excitation-autoionization cross sections starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        pathways = JAC.ImpactExcitationAutoion.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.ImpactExcitationAutoion.displayPathways(pathways)    end
        # Calculate all amplitudes and requested properties
        newPathways = ImpactExcitationAutoion.Pathway[]
        for  pathway in pathways
            newPathway = JAC.ImpactExcitationAutoion.computeAmplitudesProperties(pathway, grid, settings) 
            push!( newPathways, newPathway)
        end
        # Print all results to screen
        JAC.ImpactExcitationAutoion.displayResults(pathways)
        #
        if    output    return( pathways )
        else            return( nothing )
        end
    end


    """
    `JAC.ImpactExcitationAutoion.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                   settings::ImpactExcitationAutoion.Settings)`  ... to determine a list of electron-impact
         excitation-autoionization pathways between the levels from the given initial-, intermediate- and final-state multiplets and by taking 
         into account the particular selections and settings for this computation; an Array{ImpactExcitationAutoion.Line,1} is returned. Apart 
         from the level specification, all physical properties are set to zero during the initialization process.  
    """
    function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                settings::ImpactExcitationAutoion.Settings)
        if    settings.selectPathways    selectPathways = true;   selectedPathways = JAC.determine("selected pathways", settings.selectedPathways)
        else                             selectPathways = false
        end
    
        pathways = ImpactExcitationAutoion.Pathway[]
        for  i = 1:length(initialMultiplet.levels)
            for  n = 1:length(intermediateMultiplet.levels)
                for  f = 1:length(finalMultiplet.levels)
                    if  selectPathways  &&  !(haskey(selectedPathways, (i,n,f)) )    continue   end
                    for  en in settings.electronEnergies
                        eInEnergy  = en
                        eOutEnergy = en - (intermediateMultiplet.levels[n].energy - initialMultiplet.levels[i].energy)
                        eAEnergy   = intermediateMultiplet.levels[n].energy - finalMultiplet.levels[f].energy
                        if  eOutEnergy < 0.   ||   eAEnergy < 0    continue    end

                        channels = JAC.ImpactExcitationAutoion.determineChannels(finalMultiplet.levels[f], intermediateMultiplet.levels[n], 
                                                                                 initialMultiplet.levels[i], settings) 
                        push!( pathways, ImpactExcitationAutoion.Pathway(initialMultiplet.levels[i], intermediateMultiplet.levels[n], 
                                                                         finalMultiplet.levels[f], eInEnergy, eOutEnergy, eAEnergy, 
                                                                         0., true, channels) )
                    end
                end
            end
        end
        return( pathways )
    end


    """
    `JAC.ImpactExcitationAutoion.determineChannels(finalLevel::Level, intermediateLevel::Level, initialLevel::Level, 
                                                   settings::ImpactExcitationAutoion.Settings)`  ... to determine a list of electron-impact
                                                   excitation-autoionization Channels for a pathway from the initial to an intermediate and 
                                                   to a final level, and by taking into account the particular settings of for this computation; 
                                                   an Array{ImpactExcitationAutoion.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, intermediateLevel::Level, initialLevel::Level, settings::ImpactExcitationAutoion.Settings)
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        symn      = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
        # Determine first the electron-impact excitation channels
        eChannels = ImpactExcitation.Channel[];   
        for  inKappa = -(settings.maxKappa+1):settings.maxKappa
            if  inKappa == 0    continue    end
            symtList = JAC.AngularMomentum.allowedTotalSymmetries(symi, inKappa)
            for  symt in symtList
                outKappaList = JAC.AngularMomentum.allowedKappaSymmetries(symt, symn)
                for  outKappa in outKappaList
                    push!(eChannels, ImpactExcitation.Channel( inKappa, outKappa, symt, 0., 0.,Complex(0.)) )
                end
            end
        end

        # Determine next the AutoIonization channels
        aChannels = AutoIonization.Channel[];   
        kappaList = JAC.AngularMomentum.allowedKappaSymmetries(symn, symf)
        for  kappa in kappaList
            push!(aChannels, AutoIonization.Channel(kappa, symn, 0., Complex(0.)) )
        end

        # Now combine all these channels
        channels  = ImpactExcitationAutoion.Channel[]; 
        for    ec in eChannels  
            for    a in aChannels
                push!(channels,  ImpactExcitationAutoion.Channel(ec, a) )    
            end
        end
 
        return( channels )  
    end


    """
    `JAC.ImpactExcitationAutoion.displayPathways(pathways::Array{ImpactExcitationAutoion.Line,1})`  ... to display a list of pathways and 
         channels that have been selected due to the prior settings. A neat table of all selected transitions and energies is printed but 
         nothing is returned otherwise.
    """
    function  displayPathways(pathways::Array{ImpactExcitationAutoion.Pathway,1})
        println(" ")
        println("  Selected electron-impact excitation-autoionization pathways:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(179))
        sa = "     ";   sb = "     "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=2);            sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=3);
        sa = sa * JAC.TableStrings.center(58, "Energies  " * JAC.TableStrings.inUnits("energy");           na=4);              
        sb = sb * JAC.TableStrings.center(58, "  m--i        m--f     e^-(in)     e^-(out)    e^-(Auger)"; na=4)
        sa = sa * JAC.TableStrings.flushleft(57, "List of partial waves and total symmetries"; na=4)  
        sb = sb * JAC.TableStrings.flushleft(57, "(partial_in -> partial_out) J_total^P [partial_Auger] "; na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(179)) 
        #   
        for  pathway in pathways
            sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
            en_mf = pathway.intermediateLevel.energy - pathway.finalLevel.energy
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", en_mi))   * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", en_mf))   * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", pathway.electronInEnergy))    * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", pathway.electronOutEnergy))   * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", pathway.electronAugerEnergy)) * "    "
            kappaSymmetryList = Tuple{Int64,Int64,Int64,LevelSymmetry}[]
            for  i in 1:length(pathway.channels)
                eChannel = pathway.channels[i].excitationChannel;    aChannel = pathway.channels[i].augerChannel;  
                push!( kappaSymmetryList, (eChannel.initialKappa, eChannel.finalKappa, aChannel.kappa, eChannel.symmetry) )
            end
            wa = JAC.TableStrings.kappaKappaKappaSymmetryTupels(57, kappaSymmetryList)
            if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 2:length(wa)
                sb = JAC.TableStrings.hBlank( length(sa) )  
                if       i == 7   sb = sb * "... " * string( length(kappaSymmetryList) ) * " in total"
                elseif   i >  7   continue
                else     sb = sb * wa[i]
                end
                println( sb )
            end
        end
        println("  ", JAC.TableStrings.hLine(179))
        #
        return( nothing )
    end


    """
    `JAC.ImpactExcitationAutoion.displayResults(pathways::Array{ImpactExcitationAutoion.Line,1})`  ... to list all results, energies, 
         cross sections, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(pathways::Array{ImpactExcitationAutoion.Pathway,1})
        println(" ")
        println("  Electron-impact excitation-autoionization cross sections:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(133))
        sa = "     ";   sb = "     "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=2);            sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=3);
        sa = sa * JAC.TableStrings.center(58, "Energies  " * JAC.TableStrings.inUnits("energy");           na=3);              
        sb = sb * JAC.TableStrings.center(58, "  m--i        m--f     e^-(in)     e^-(out)    e^-(Auger)"; na=3)
        sa = sa * JAC.TableStrings.center(16, "Cross sections"; na=2);       
        sb = sb * JAC.TableStrings.center(16, JAC.TableStrings.inUnits("cross section"); na=2)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(133)) 
        #   
        for  pathway in pathways
            sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            en_mi = pathway.intermediateLevel.energy - pathway.initialLevel.energy
            en_mf = pathway.intermediateLevel.energy - pathway.finalLevel.energy
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", en_mi))   * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", en_mf))   * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", pathway.electronInEnergy))      * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", pathway.electronOutEnergy))     * "  "
            sa = sa * @sprintf("%.4e", Basics.convert("energy: from atomic", pathway.electronAugerEnergy))   * "     "
            sa = sa * @sprintf("%.6e", Basics.convert("cross section: from atomic", pathway.crossSection))   * "  "
            println(sa)
        end
        println("  ", JAC.TableStrings.hLine(133))
        #
        return( nothing )
    end

end # module


