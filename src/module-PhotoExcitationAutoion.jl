
"""
`module  JAC.PhotoExcitationAutoion`  ... a submodel of JAC that contains all methods for computing photo-excitation-autoionization cross 
                                          sections and rates; it is using JAC, JAC.ManyElectron, JAC.Radial, JAC.Radiative, JAC.Auger.
"""
module PhotoExcitationAutoion 

    using Printf, JAC, JAC.ManyElectron, JAC.Radial, JAC.Radiative, JAC.Auger
    global JAC_counter = 0


    """
    `struct  PhotoExcitationAutoion.Settings`  ... defines a type for the details and parameters of computing photon-impact 
                                                   excitation-autoionization pathways |i(N)>  --> |m(N)>  --> |f(N-1)>.

        + multipoles              ::Array{JAC.EmMultipole,1}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{JAC.UseGauge,1}              ... Specifies the gauges to be included into the computations.
        + printBeforeComputation  ::Bool                               ... True, if all energies and lines are printed before their evaluation.
        + selectPathways          ::Bool                               ... True if particular pathways are selected for the computations.
        + selectedPathways        ::Array{Tuple{Int64,Int64,Int64},1}  ... List of list of pathways, given by tupels (inital, inmediate, final).
        + maxKappa                ::Int64                              ... Maximum kappa value of partial waves to be included.
    """
    struct Settings
        multipoles                ::Array{JAC.EmMultipole,1}
        gauges                    ::Array{JAC.UseGauge,1} 
        printBeforeComputation    ::Bool
        selectPathways            ::Bool
        selectedPathways          ::Array{Tuple{Int64,Int64,Int64},1}
        maxKappa                  ::Int64 
    end 


    """
    `JAC.PhotoExcitationAutoion.Settings()`  ... constructor for the default values of photon-impact excitation-autoionizaton settings.
    """
    function Settings()
        Settings( JAC.EmMultipole[], UseGauge[], false,  false, Tuple{Int64,Int64,Int64}[], 0)
    end


    """
    `Base.show(io::IO, settings::PhotoExcitationAutoion.Settings)`  ... prepares a proper printout of the variable 
                                                                        settings::PhotoExcitationAutoion.Settings.  
    """
    function Base.show(io::IO, settings::PhotoExcitationAutoion.Settings) 
        println(io, "multipoles:              $(settings.multipoles)  ")
        println(io, "gauges:                  $(settings.gauges)  ")
        println(io, "printBeforeComputation:  $(settings.printBeforeComputation)  ")
        println(io, "selectPathways:          $(settings.selectPathways)  ")
        println(io, "selectedPathways:        $(settings.selectedPathways)  ")
        println(io, "maxKappa:                $(settings.maxKappa)  ")
    end



    """
    `struct  JAC.PhotoExcitationAutoion.Channel`  ... defines a type for a photon-impact excitaton & autoionization channel that specifies 
                                                      all quantum numbers, phases and amplitudes.

        + excitationChannel  ::JAC.Radiative.Channel       ... Channel that describes the photon-impact excitation process.
        + augerChannel       ::JAC.Auger.Channel           ... Channel that describes the subsequent Auger/autoionization process.
    """
    struct  Channel
        excitationChannel    ::JAC.Radiative.Channel
        augerChannel         ::JAC.Auger.Channel
    end 


    """
    `struct  JAC.PhotoExcitationAutoion.Pathway`  ... defines a type for a photon-impact excitation pathway that may include the definition 
                                                      of different excitation and autoionization channels and their corresponding amplitudes.

        + initialLevel        ::Level                  ... initial-(state) level
        + intermediateLevel   ::Level                  ... intermediate-(state) level
        + finalLevel          ::Level                  ... final-(state) level
        + photonEnergy        ::Float64                 ... energy of the (incoming) electron
        + electronEnergy      ::Float64                 ... energy of the (finally outgoing, scattered) electron
        + crossSection        ::EmProperty              ... total cross section of this pathway
        + hasChannels         ::Bool                    ... Determines whether the individual excitation and autoionization channels are defined 
                                                            in terms of their multipole, gauge, free-electron kappa, phases and the total 
                                                            angular momentum/parity as well as the amplitude, or not.
        + channels            ::Array{PhotoExcitationAutoion.Channel,1}  ... List of channels of this pathway.
    """
    struct  Pathway
        initialLevel          ::Level
        intermediateLevel     ::Level
        finalLevel            ::Level
        photonEnergy          ::Float64
        electronEnergy        ::Float64
        crossSection          ::EmProperty
        hasChannels           ::Bool
        channels              ::Array{PhotoExcitationAutoion.Channel,1}
    end 


    """
    `JAC.PhotoExcitationAutoion.Pathway()`  ... 'empty' constructor for an photon-impact excitation-autoionization pathway between a specified 
                                                initial, intermediate and final level.
    """
    function Pathway()
        Pathway(Level(), Level(), Level(), 0., 0., 0., false, PhotoExcitationAutoion.Channel[] )
    end


    """
    `Base.show(io::IO, pathway::PhotoExcitationAutoion.Pathway)`  ... prepares a proper printout of the variable 
                                                                      pathway::PhotoExcitationAutoion.Pathway.
    """
    function Base.show(io::IO, pathway::PhotoExcitationAutoion.Pathway) 
        println(io, "initialLevel:               $(pathway.initialLevel)  ")
        println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
        println(io, "finalLevel:                 $(pathway.finalLevel)  ")
        println(io, "photonEnergy                $(pathway.photonEnergy)  ") 
        println(io, "electronEnergy              $(pathway.electronEnergy)  ")
        println(io, "crossSection:               $(pathway.crossSection)  ")
        println(io, "hasChannels:                $(pathway.hasChannels)  ")
        println(io, "channels:                   $(pathway.channels)  ")
    end



    """
    `JAC.PhotoExcitationAutoion.computeAmplitudesProperties(pathway::PhotoExcitationAutoion.Pathway, grid::Radial.Grid, 
                                                            settings::PhotoExcitationAutoion.Settings)` ... to compute all amplitudes and 
         properties of the given line; a line::PhotoExcitationAutoion.Line is returned for which the amplitudes and properties have now been 
         evaluated.
    """
    function  computeAmplitudesProperties(pathway::PhotoExcitationAutoion.Pathway, grid::Radial.Grid, settings::PhotoExcitationAutoion.Settings)
        global JAC_counter
        newChannels = PhotoExcitationAutoion.Channel[]
        for channel in pathway.channels
            # Generate a continuum orbital
            JAC_counter = JAC_counter + 1
            if   JAC_counter < 20   println("PhotoExcitationAutoion.computeAmplitudesProperties-aa: warning ... no cont. orb. is generated.")  end
            phase  = 0.
            # Define a proper continuum basis from the finalLevel.basis and the continuum orbital
            JAC_counter = JAC_counter + 1
            if   JAC_counter < 20   println("PhotoExcitationAutoion.computeAmplitudesProperties-ab: warning ... no cont. basis generated.") end
            # Compute the transition matrix for the continuum and the initial-state basis
            JAC_counter = JAC_counter + 1
            if   JAC_counter < 20   println("PhotoExcitationAutoion.computeAmplitudesProperties-ac: warning ... no trans-matrix computed.") end
            # matrix    = JAC.PhotoExcitationAutoion.computeMatrix(channel.multipole, channel.gauge, line.omega, line.finalLevel.basis, 
            #                                                      line.initialLevel.basis, grid, settings)
            # amplitude = line.finalLevel.mc * matrix * line.initialLevel.mc 
            amplitude = 1.0 
            eChannel  = Radiative.Channel( channel.excitationChannel.multipole, channel.excitationChannel.gauge, 0.)
            aChannel  = Auger.Channel( channel.augerChannel.kappa, channel.augerChannel.symmetry, 0., Complex(0.))
            push!( newChannels, PhotoExcitationAutoion.Channel(eChannel, aChannel) )
        end
        # Calculate the totalRate 
        JAC_counter = JAC_counter + 1
        if   JAC_counter < 20   println("PhotoExcitationAutoion.computeAmplitudesProperties-ba: warning ... crossSection set to -1.") end
        crossSection = EmProperty(-1., -1.)
        pathway = PhotoExcitationAutoion.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, 
                                                  pathway.photonEnergy, pathway.electronEnergy, crossSection, true, newChannels)
        return( pathway )
    end



    """
    `JAC.PhotoExcitationAutoion.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                grid::Radial.Grid, settings::PhotoExcitation.Settings; output=true)`  ... to compute the 
         photo-excitation-autoionization amplitudes and all properties as requested by the given settings. A list of 
         lines::Array{PhotoExcitationAutoion.Lines} is returned.
    """
    function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                              settings::PhotoExcitationAutoion.Settings; output=true)
        println("")
        printstyled("JAC.PhotoExcitationAutoion.computePathways(): The computation of photo-excitation-autoionization amplitudes starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        pathways = JAC.PhotoExcitationAutoion.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.PhotoExcitationAutoion.displayPathways(pathways)    end
        # Calculate all amplitudes and requested properties
        newPathways = PhotoExcitationAutoion.Pathway[]
        for  pathway in pathways
            newPathway = JAC.PhotoExcitationAutoion.computeAmplitudesProperties(pathway, grid, settings) 
            push!( newPathways, newPathway)
        end
        # Print all results to screen
        JAC.PhotoExcitationAutoion.displayResults(stdout, pathways)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary    JAC.PhotoExcitationAutoion.displayResults(iostream, pathways)   end
        #
        if    output    return( pathways )
        else            return( nothing )
        end
    end


    """
    `JAC.PhotoExcitationAutoion.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                  settings::PhotoExcitationAutoion.Settings)`  ... to determine a list of photoexcitation-
         autoionization pathways between the levels from the given initial-, intermediate- and final-state multiplets and by taking into account 
         the particular selections and settings for this computation; an Array{PhotoExcitationAutoion.Line,1} is returned. Apart from the 
         level specification, all physical properties are set to zero during the initialization process.  
    """
    function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                settings::PhotoExcitationAutoion.Settings)
        if    settings.selectPathways    selectPathways = true;   selectedPathways = JAC.determine("selected pathways", settings.selectedPathways)
        else                             selectPathways = false
        end
    
        pathways = PhotoExcitationAutoion.Pathway[]
        for  i = 1:length(initialMultiplet.levels)
            for  n = 1:length(intermediateMultiplet.levels)
                for  f = 1:length(finalMultiplet.levels)
                    if  selectPathways  &&  !(haskey(selectedPathways, (i,n,f)) )    continue   end
                    ##x println("PhotoExcitationAutoion.determineLines-aa: angular i = $i, f = $f")
                    pEnergy = intermediateMultiplet.levels[n].energy - initialMultiplet.levels[i].energy
                    eEnergy = intermediateMultiplet.levels[n].energy - finalMultiplet.levels[f].energy
                    if  pEnergy < 0.   ||   eEnergy < 0    continue    end

                    channels = JAC.PhotoExcitationAutoion.determineChannels(finalMultiplet.levels[f], intermediateMultiplet.levels[n], 
                                                                            initialMultiplet.levels[i], settings) 
                    push!( pathways, PhotoExcitationAutoion.Pathway(initialMultiplet.levels[i], intermediateMultiplet.levels[i], 
                                                            finalMultiplet.levels[f], pEnergy, eEnergy, EmProperty(0., 0.), true, channels) )
                end
            end
        end
        return( pathways )
    end


    """
    `JAC.PhotoExcitationAutoion.determineChannels(finalLevel::Level, intermediateLevel::Level, initialLevel::Level, 
                                                  settings::PhotoExcitationAutoion.Settings)`  ... to determine a list of photoexcitation-
         autoionization Channels for a pathway from the initial to and intermediate and to a final level, and by taking into account the 
         particular settings of for this computation; an Array{PhotoExcitationAutoion.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, intermediateLevel::Level, initialLevel::Level, settings::PhotoExcitationAutoion.Settings)
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        symn      = LevelSymmetry(intermediateLevel.J, intermediateLevel.parity)
        # Determine first the radiative channels
        rChannels = Radiative.Channel[];   
        for  mp in settings.multipoles
            if   JAC.AngularMomentum.isAllowedMultipole(symi, mp, symn)
                hasMagnetic = false
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp)[1] == 'E'  &&   gauge == JAC.UseCoulomb      push!(rChannels, Radiative.Channel(mp, JAC.Coulomb,   0.) )
                    elseif string(mp)[1] == 'E'  &&   gauge == JAC.UseBabushkin    push!(rChannels, Radiative.Channel(mp, JAC.Babushkin, 0.) )  
                    elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)               push!(rChannels, Radiative.Channel(mp, JAC.Magnetic,  0.) );
                                                        hasMagnetic = true; 
                    end 
                end
            end
        end

        # Determine next the Auger channels
        aChannels = Auger.Channel[];   
        kappaList = JAC.AngularMomentum.allowedKappaSymmetries(symn, symf)
        for  kappa in kappaList
            push!(aChannels, Auger.Channel(kappa, symi, 0., Complex(0.)) )
        end

        # Now combine all these channels
        channels  = PhotoExcitationAutoion.Channel[]; 
        for    r in rChannels  
            for    a in aChannels    
                push!(channels,  PhotoExcitationAutoion.Channel(r, a) )    
            end
        end
 
        return( channels )  
    end


    """
    `JAC.PhotoExcitationAutoion.displayPathways(pathways::Array{PhotoExcitationAutoion.Line,1})`  ... to display a list of pathways and 
         channels that have been selected due to the prior settings. A neat table of all selected transitions and energies is printed but 
         nothing is returned otherwise.
    """
    function  displayPathways(pathways::Array{PhotoExcitationAutoion.Pathway,1})
        println(" ")
        println("  Selected photo-excitation-autoionization pathways:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(170))
        sa = "     ";   sb = "     "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=2);            sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=3);
        sa = sa * JAC.TableStrings.center(14, "Energy m-i"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(14, "Energy m-f"; na=3);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=3)
        sa = sa * JAC.TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
        sb = sb * JAC.TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(170)) 
        #   
        for  pathway in pathways
            sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", pathway.photonEnergy))   * "    "
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", pathway.electronEnergy)) * "   "
            kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
            for  i in 1:length(pathway.channels)
                eChannel = pathway.channels[i].excitationChannel;    aChannel = pathway.channels[i].augerChannel;  
                push!( kappaMultipoleSymmetryList, (aChannel.kappa, eChannel.multipole, eChannel.gauge, aChannel.symmetry) )
            end
            ##x println("PhotoExcitationAutoion-diplayLines-ad: kappaMultipoleSymmetryList = ", kappaMultipoleSymmetryList)
            wa = JAC.TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
            if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end  
            for  i = 2:length(wa)
                sb = JAC.TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
        end
        println("  ", JAC.TableStrings.hLine(170))
        #
        return( nothing )
    end


    """
    `JAC.PhotoExcitationAutoion.displayResults(stream::IO, pathways::Array{PhotoExcitationAutoion.Line,1})`  ... to list all results, energies, 
         cross sections, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, pathways::Array{PhotoExcitationAutoion.Pathway,1})
        println(stream, " ")
        println(stream, "  Photo-excitation & autoionization cross sections:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(135))
        sa = "     ";   sb = "     "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=2);            sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=2);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=3);    sb = sb * JAC.TableStrings.center(23, "i  --  m  --  f"; na=3);
        sa = sa * JAC.TableStrings.center(14, "Energy m-i"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(14, "Energy m-f"; na=3);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=3)
        sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(30, "Cou -- Cross sections -- Bab"; na=2);       
        sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("cross section")*"          "*
                                              JAC.TableStrings.inUnits("cross section"); na=2)
        println(stream, sa);    println(stream, sb);    println("  ", JAC.TableStrings.hLine(135)) 
        #   
        for  pathway in pathways
            sa  = "  ";    isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", pathway.photonEnergy))   * "    "
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", pathway.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ch in pathway.channels
                multipoles = push!( multipoles, ch.excitationChannel.multipole)
            end
            multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles) * "          "
            sa = sa * JAC.TableStrings.flushleft(11, mpString[1:10];  na=3)
            sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", pathway.crossSection.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", pathway.crossSection.Babushkin))   * "    "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(135))
        #
        return( nothing )
    end

end # module
