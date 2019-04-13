
"""
`module  JAC.PhotoExcitationFluores`  ... a submodel of JAC that contains all methods for computing photo-excitation-autoionization cross 
                                          sections and rates; it is using JAC, JAC.ManyElectron, JAC.Radial, JAC.PhotoEmission, JAC.AutoIonization.
"""
module PhotoExcitationFluores 

    using Printf, JAC, JAC.ManyElectron, JAC.Radial, JAC.PhotoEmission, JAC.AutoIonization
    global JAC_counter = 0


    """
    `struct  PhotoExcitationFluores.Settings`  ... defines a type for the details and parameters of computing photon-impact 
                                                   excitation-autoionization pathways |i(N)>  --> |m(N)>  --> |f(N-1)>.

        + multipoles              ::Array{JAC.EmMultipole,1}           ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{JAC.UseGauge,1}              ... Specifies the gauges to be included into the computations.
        + printBeforeComputation  ::Bool                               ... True, if all energies and lines are printed before their evaluation.
        + selectPathways          ::Bool                               ... True if particular pathways are selected for the computations.
        + selectedPathways        ::Array{Tuple{Int64,Int64,Int64},1}  ... List of list of pathways, given by tupels (inital, inmediate, final).
    """
    struct Settings
        multipoles                ::Array{JAC.EmMultipole,1}
        gauges                    ::Array{JAC.UseGauge,1} 
        printBeforeComputation    ::Bool
        selectPathways            ::Bool
        selectedPathways          ::Array{Tuple{Int64,Int64,Int64},1}
     end 


    """
    `JAC.PhotoExcitationFluores.Settings()`  ... constructor for the default values of photon-impact excitation-autoionizaton settings.
    """
    function Settings()
        Settings( JAC.EmMultipole[], UseGauge[], false,  false, Tuple{Int64,Int64,Int64}[])
    end


    # `Base.show(io::IO, settings::PhotoExcitationFluores.Settings)` 
    # 	 ... prepares a proper printout of the variable settings::PhotoExcitationFluores.Settings.  
    function Base.show(io::IO, settings::PhotoExcitationFluores.Settings) 
        println(io, "multipoles:              $(settings.multipoles)  ")
        println(io, "gauges:                  $(settings.gauges)  ")
        println(io, "printBeforeComputation:  $(settings.printBeforeComputation)  ")
        println(io, "selectPathways:          $(settings.selectPathways)  ")
        println(io, "selectedPathways:        $(settings.selectedPathways)  ")
    end


    #==
    """
    `struct  JAC.PhotoExcitationFluores.Channel`  ... defines a type for a photon-impact excitaton & autoionization channel that specifies 
                                                      all quantum numbers, phases and amplitudes.

        + excitationChannel  ::JAC.PhotoEmission.Channel       ... Channel that describes the photon-impact excitation process.
        + augerChannel       ::JAC.AutoIonization.Channel           ... Channel that describes the subsequent Auger/autoionization process.
    """
    struct  Channel
        excitationChannel    ::JAC.PhotoEmission.Channel
        augerChannel         ::JAC.AutoIonization.Channel
    end ==#


    """
    `struct  JAC.PhotoExcitationFluores.Pathway`  ... defines a type for a photon-impact excitation pathway that may include the definition 
                                                      of different excitation and autoionization channels and their corresponding amplitudes.

        + initialLevel        ::Level                  ... initial-(state) level
        + intermediateLevel   ::Level                  ... intermediate-(state) level
        + finalLevel          ::Level                  ... final-(state) level
        + excitEnergy         ::Float64                ... photon excitation energy of this pathway
        + fluorEnergy         ::Float64                ... photon fluorescence energy of this pathway
        + crossSection        ::EmProperty             ... total cross section of this pathway
        + hasChannels         ::Bool                   ... Determines whether the individual excitation and fluorescence channels are defined 
                                                           in terms of their multipole, gauge as well as the amplitude, or not.
        + excitChannels       ::Array{JAC.PhotoEmission.Channel,1}  ... List of excitation channels of this pathway.
        + fluorChannels       ::Array{JAC.PhotoEmission.Channel,1}  ... List of fluorescence channels of this pathway.
    """
    struct  Pathway
        initialLevel          ::Level
        intermediateLevel     ::Level
        finalLevel            ::Level
        excitEnergy           ::Float64
        fluorEnergy           ::Float64
        crossSection          ::EmProperty
        hasChannels           ::Bool
        excitChannels         ::Array{JAC.PhotoEmission.Channel,1}
        fluorChannels         ::Array{JAC.PhotoEmission.Channel,1}
    end 


    """
    `JAC.PhotoExcitationFluores.Pathway()`  ... 'empty' constructor for an photon-impact excitation-autoionization pathway between a specified 
                                                initial, intermediate and final level.
    """
    function Pathway()
        Pathway(Level(), Level(), Level(), 0., 0., 0., false, PhotoExcitationFluores.Channel[] )
    end


    # `Base.show(io::IO, pathway::PhotoExcitationFluores.Pathway)` 
    #		 ... prepares a proper printout of the variable pathway::PhotoExcitationFluores.Pathway.
    function Base.show(io::IO, pathway::PhotoExcitationFluores.Pathway) 
        println(io, "initialLevel:               $(pathway.initialLevel)  ")
        println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
        println(io, "finalLevel:                 $(pathway.finalLevel)  ")
        println(io, "excitEnergy                 $(pathway.excitEnergy)  ") 
        println(io, "fluorEnergy                 $(pathway.fluorEnergy)  ")
        println(io, "crossSection:               $(pathway.crossSection)  ")
        println(io, "hasChannels:                $(pathway.hasChannels)  ")
        println(io, "excitChannels:              $(pathway.excitChannels)  ")
        println(io, "fluorChannels:              $(pathway.fluorChannels)  ")
    end


    """
    `JAC.PhotoExcitationFluores.computeAmplitudesProperties(pathway::PhotoExcitationFluores.Pathway, grid::Radial.Grid, 
                                                            settings::PhotoExcitationFluores.Settings)` 
        ... to compute all amplitudes and properties of the given line; a pathway::PhotoExcitationFluores.Pathway is returned for which the 
            amplitudes and properties have now been evaluated.
    """
    function  computeAmplitudesProperties(pathway::PhotoExcitationFluores.Pathway, grid::Radial.Grid, settings::PhotoExcitationFluores.Settings)
        # Compute all excitation channels
        neweChannels = PhotoEmission.Channel[]
        for eChannel in pathway.excitChannels
            amplitude   = JAC.PhotoEmission.amplitude("absorption", eChannel.multipole, eChannel.gauge, pathway.excitEnergy, 
                                                  pathway.intermediateLevel, pathway.initialLevel, grid)
             push!( neweChannels, PhotoEmission.Channel( eChannel.multipole, eChannel.gauge, amplitude))
        end
        # Compute all fluorescence channels
        newfChannels = PhotoEmission.Channel[]
        for fChannel in pathway.fluorChannels
            amplitude   = JAC.PhotoEmission.amplitude("emission", fChannel.multipole, fChannel.gauge, pathway.fluorEnergy, 
                                                  pathway.finalLevel, pathway.intermediateLevel, grid)
            push!( newfChannels, PhotoEmission.Channel( fChannel.multipole, fChannel.gauge, amplitude))
        end
        crossSection = EmProperty(-1., -1.)
        pathway = PhotoExcitationFluores.Pathway( pathway.initialLevel, pathway.intermediateLevel, pathway.finalLevel, pathway.excitEnergy, 
                                                  pathway.fluorEnergy, crossSection, true, neweChannels, newfChannels)
        return( pathway )
    end



    """
    `JAC.PhotoExcitationFluores.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                grid::Radial.Grid, settings::PhotoExcitationFluores.Settings; output=true)`  
        ... to compute the photo-excitation-fluorescence amplitudes and all properties as requested by the given settings. A list of
            lines::Array{PhotoExcitationFluores.Lines} is returned.
    """
    function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                              settings::PhotoExcitationFluores.Settings; output=true)
        println("")
        printstyled("JAC.PhotoExcitationFluores.computePathways(): The computation of photo-excitation-fluorescence amplitudes starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------------------------------------ \n", color=:light_green)
        println("")
        pathways = JAC.PhotoExcitationFluores.determinePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, settings)
        # Display all selected pathways before the computations start
        if  settings.printBeforeComputation    JAC.PhotoExcitationFluores.displayPathways(pathways)    end
        # Calculate all amplitudes and requested properties
        newPathways = PhotoExcitationFluores.Pathway[]
        ##x println("computePathways: NO-pathwayS = $(length(pathways)) ")
        for  pathway in pathways
            ##x println("computePathways: pathway = $pathway ")
            newPathway = JAC.PhotoExcitationFluores.computeAmplitudesProperties(pathway, grid, settings) 
            push!( newPathways, newPathway)
        end
        # Print all results to screen
        JAC.PhotoExcitationFluores.displayResults(stdout, pathways, settings)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary   JAC.PhotoExcitationFluores.displayResults(iostream, pathways, settings)     end
        #
        if    output    return( pathways )
        else            return( nothing )
        end
    end


    """
    `JAC.PhotoExcitationFluores.determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                  settings::PhotoExcitationFluores.Settings)`  
        ... to determine a list of dielectronic-recombination pathways between the levels from the given initial-, intermediate- and 
            final-state multiplets and by taking into account the particular selections and settings for this computation; 
            an Array{PhotoExcitationFluores.Pathway,1} is returned. Apart from the level specification, all physical properties are 
            set to zero during the initialization process.  
    """
    function  determinePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                settings::PhotoExcitationFluores.Settings)
        if    settings.selectPathways    selectPathways = true;    
              selectedPathways = JAC.determineSelectedPathways(settings.selectedPathways, initialMultiplet, intermediateMultiplet, finalMultiplet)
        else                             selectPathways = false
        end
    
        pathways = PhotoExcitationFluores.Pathway[]
        for  i = 1:length(initialMultiplet.levels)
            for  n = 1:length(intermediateMultiplet.levels)
                for  f = 1:length(finalMultiplet.levels)
                    if  selectPathways  &&  !((i,n,f) in selectedPathways )    continue   end
                    eEnergy = intermediateMultiplet.levels[n].energy - initialMultiplet.levels[i].energy
                    fEnergy = intermediateMultiplet.levels[n].energy - finalMultiplet.levels[f].energy
                    if  eEnergy < 0.   ||   fEnergy < 0.    continue    end

                    rSettings = JAC.PhotoEmission.Settings( settings.multipoles, settings.gauges, false, false, false, Tuple{Int64,Int64}[], 0., 0., 0.)
                    eChannels = JAC.PhotoEmission.determineChannels(intermediateMultiplet.levels[n], initialMultiplet.levels[i], rSettings) 
                    fChannels = JAC.PhotoEmission.determineChannels(finalMultiplet.levels[f], intermediateMultiplet.levels[n], rSettings) 
                    push!( pathways, PhotoExcitationFluores.Pathway(initialMultiplet.levels[i], intermediateMultiplet.levels[i], finalMultiplet.levels[f], 
                                                                    eEnergy, fEnergy, EmProperty(0., 0.), true, eChannels, fChannels) )
                end
            end
        end
        return( pathways )
    end


    """
    `JAC.PhotoExcitationFluores.displayResults(stream::IO, pathways::Array{PhotoExcitationFluores.Pathway,1}, 
                                               settings::PhotoExcitationFluores.Settings)`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing 
            is returned otherwise.
    """
    function  displayResults(stream::IO, pathways::Array{PhotoExcitationFluores.Pathway,1}, settings::PhotoExcitationFluores.Settings)
        println(stream, " ")
        println(stream, "  Partial excitation & fluorescence cross sections:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(150))
        sa = "    ";   sb = "    "
        sa = sa * JAC.TableStrings.center(23, "Levels"; na=2);            sb = sb * JAC.TableStrings.center(23, "i  --  e  --  f"; na=2);          
        sa = sa * JAC.TableStrings.center(23, "J^P symmetries"; na=0);    sb = sb * JAC.TableStrings.center(23, "i  --  e  --  f"; na=2);
        sa = sa * JAC.TableStrings.center(30, "Energies  " * JAC.TableStrings.inUnits("energy"); na=4);              
        sb = sb * JAC.TableStrings.center(30, "r--i        r--f  "; na=4)
        sa = sa * JAC.TableStrings.center(30, "Multipoles"; na=4);        
        sb = sb * JAC.TableStrings.center(38, "i--r        r--f  "; na=4)
        sa = sa * JAC.TableStrings.center(36, "Cou -- cross sections -- Bab";        na=2)
        sb = sb * JAC.TableStrings.center(36, JAC.TableStrings.inUnits("cross section") * "          " * 
                                              JAC.TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(150)) 
        #   
        for  pathway in pathways
            sa  = " ";     isym = LevelSymmetry( pathway.initialLevel.J,      pathway.initialLevel.parity)
                           msym = LevelSymmetry( pathway.intermediateLevel.J, pathway.intermediateLevel.parity)
                           fsym = LevelSymmetry( pathway.finalLevel.J,        pathway.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.levels_imf(pathway.initialLevel.index, pathway.intermediateLevel.index, 
                                                                              pathway.finalLevel.index); na=3)
            sa = sa * JAC.TableStrings.center(23, JAC.TableStrings.symmetries_imf(isym, msym, fsym);  na=4)
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", pathway.excitEnergy))      * "   "
            sa = sa * @sprintf("%.4e", JAC.convert("energy: from atomic", pathway.fluorEnergy))      * "   "
            #
            multipoles = EmMultipole[]
            for  ech in pathway.excitChannels
                multipoles = push!( multipoles, ech.multipole)
            end
            multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles)       * "                   "
            sa = sa * JAC.TableStrings.flushleft(16, mpString[1:16];  na=2)
            #
            multipoles = EmMultipole[]
            for  fch in pathway.fluorChannels
                multipoles = push!( multipoles, fch.multipole)
            end
            multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles)       * "                   "
            sa = sa * JAC.TableStrings.flushleft(16, mpString[1:16];  na=2)
            sa = sa * @sprintf("%.4e", JAC.convert("cross section: from atomic", pathway.crossSection.Coulomb))     * "  "
            sa = sa * @sprintf("%.4e", JAC.convert("cross section: from atomic", pathway.crossSection.Babushkin))   * "  "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(150))
        #
        return( nothing )
    end

end # module
