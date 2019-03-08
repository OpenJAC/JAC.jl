
"""
`module  JAC.PhotoExcitation`  ... a submodel of JAC that contains all methods for computing photo-excitation properties between some initial 
                                   and final-state multiplets; it is using JAC, JAC.ManyElectron, JAC.Radial and JAC.Radiative.
"""
module PhotoExcitation

    using Printf, JAC, JAC.ManyElectron, JAC.Radial, JAC.Radiative
    global JAC_counter = 0


    """
    `struct  PhotoExcitation.Settings`  ... defines a type for the details and parameters of computing photo-excitation  lines.

        + multipoles              ::Array{EmMultipole,1}         ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge,1}            ... Specifies the gauges to be included into the computations.
        + calcForStokes           ::Bool                         ... True, if the excitation cross sections are to be calculated (and false otherwise)
                                                                     for given Stokes parameter of the incident plane-wave photons.
        + calcPhotonDm            ::Bool                         ... True, if the photon density matrix of a subsequently emitted fluorescence photon is to be
                                                                     calculated and false otherwise. 
        + calcTensors             ::Bool                         ... True, if the statistical tensors of the excited atom are to be calculated and false otherwise. 
        + printBeforeComputation  ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + selectLines             ::Bool                         ... True, if lines are selected individually for the computations.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines, given by tupels (inital-level, final-level).
        + photonEnergyShift       ::Float64                      ... An overall energy shift for all photon energies.
        + mimimumPhotonEnergy     ::Float64                      ... minimum transition energy for which (photon) transitions are included into the
                                                                     computation.
        + maximumPhotonEnergy     ::Float64                      ... maximum transition energy for which (photon) transitions are included.
        + stokes                  ::ExpStokes                    ... maximum transition energy for which (photon) transitions are included.
    """
    struct Settings 
        multipoles                ::Array{EmMultipole,1}
        gauges                    ::Array{UseGauge,1}
        calcForStokes             ::Bool
        calcPhotonDm              ::Bool  
        calcTensors               ::Bool  
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1}
        photonEnergyShift         ::Float64
        mimimumPhotonEnergy       ::Float64   
        maximumPhotonEnergy       ::Float64     
        stokes                    ::ExpStokes 
    end 


    """
    `JAC.PhotoExcitation.Settings()`  ... 'empty' constructor for the default values of photo-excitation line computations
    """
    function Settings()
        Settings(EmMultipole[], UseGauge[], false, false, false, false, false, Array{Tuple{Int64,Int64},1}[], 0., 0., 0.)
    end


    """
    `Base.show(io::IO, settings::PhotoExcitation.Settings)`  ... prepares a proper printout of the variable settings::PhotoExcitation.Settings.
    """
    function Base.show(io::IO, settings::PhotoExcitation.Settings) 
        println(io, "multipoles:               $(settings.multipoles)  ")
        println(io, "use-gauges:               $(settings.gauges)  ")
        println(io, "calcForStokes:            $(settings.calcForStokes)  ")
        println(io, "calcPhotonDm:             $(settings.calcPhotonDm)  ")
        println(io, "calcTensor:               $(settings.calcTensor)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLines:              $(settings.selectLines)  ")
        println(io, "selectedLines:            $(settings.selectedLines)  ")
        println(io, "photonEnergyShift:        $(settings.photonEnergyShift)  ")
        println(io, "mimimumPhotonEnergy:      $(settings.mimimumPhotonEnergy)  ")
        println(io, "maximumPhotonEnergy:      $(settings.maximumPhotonEnergy)  ")
        println(io, "stokes:                   $(settings.stokes)  ")
    end


    """
    `struct  JAC.PhotoExcitation.Line`  ... defines a type for a photo-excitation line that may include the definition of sublines and their 
                                            corresponding amplitudes.

        + initialLevel   ::Level                       ... initial-(state) level
        + finalLevel     ::Level                       ... final-(state) level
        + omega          ::Float64                     ... Transition frequency of this line; can be shifted w.r.t. the level energies.
        + crossSection   ::EmProperty                  ... Total cross section of this line.
        + staTensor      ::Array{TensorComp,1}         ... Array of statistical tensor components rho_kq
        + hasSublines    ::Bool                        ... Determines whether the individual sublines are defined in terms of their 
                                                           multipolarity, amplitude, or not; cf. JAC.Radiative.Channel
        + channels       ::Array{JAC.Radiative.Channel,1}  ... List of radiative (photon) channels
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        omega            ::Float64
        crossSection     ::EmProperty
        staTensor        ::Array{TensorComp,1}
        hasSublines      ::Bool
        channels         ::Array{Radiative.Channel,1}
    end 


    """
    `JAC.PhotoExcitation.Line(initialLevel::Level, finalLevel::Level, photonRate::Float64)`  ... constructor an photo-excitation line between a
                                                                                                 specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, omega::Float64, crossSection::EmProperty)
       Line(initialLevel, finalLevel, omega, crossSection, EmProperty(0., 0.), false, Radiative.Channel[])
    end


    """
    `Base.show(io::IO, line::PhotoExcitation.Line)`  ... prepares a proper printout of the variable line::PhotoExcitation.Line.
    """
    function Base.show(io::IO, line::PhotoExcitation.Line) 
        println(io, "initialLevel:         $(line.initialLevel)  ")
        println(io, "finalLevel:           $(line.finalLevel)  ")
        println(io, "omega:                $(line.omega)  ")
        println(io, "crossSection:         $(line.crossSection)  ")
        println(io, "staTensor:            $(line.staTensor)  ")
        println(io, "hasSublines:          $(line.hasSublines)  ")
        println(io, "channels:             $(line.channels)  ")
    end


    """
    `JAC.PhotoExcitation.computeAmplitudesProperties(line::PhotoExcitation.Line, grid::Radial.Grid, settings::PhotoExcitation.Settings)`  
         ... to compute all amplitudes and properties of the given line; a line::PhotoExcitation.Line is returned for which the amplitudes and 
         properties have now been evaluated.
    """
    function  computeAmplitudesProperties(line::PhotoExcitation.Line, grid::Radial.Grid, settings::PhotoExcitation.Settings)
        newChannels = Radiative.Channel[]
        for  channel  in  line.channels
            ##x println("channel = $channel ")
            amplitude = JAC.Radiative.amplitude("absorption", channel.multipole, channel.gauge, line.omega, line.finalLevel, line.initialLevel, grid)
            push!( newChannels, Radiative.Channel( channel.multipole, channel.gauge, amplitude) )
        end
        # Calculate the crossSection and stastistical tensors if requested
        crossSection = EmProperty(-1., -1.)    # Cross sections for linear-polarized plane-wave photons
        staTensor    = JAC.TensorComp[]
        if  settings.calcTensors
            push!(staTensor, JAC.TensorComp(0, 0, 1.))
        end
        line = PhotoExcitation.Line( line.initialLevel, line.finalLevel, line.omega, crossSection, staTensor, true, newChannels)
        return( line )
    end


    """
    `JAC.PhotoExcitation.computeCrossSection(line::PhotoExcitation.Line, stokes::ExpStokes)`  
         ... to compute the excitation cross section for the excitation of unpolarized atoms by plane-wave photons, whose polarization is described
             by the given (experimental) Stokes parameters. A cs::EmProperty is returned.
    """
    function  computeCrossSection(line::PhotoExcitation.Line, stokes::ExpStokes)
        csCoulomb = -2.;    csBabushkin = -2.
        cs = EmProperty( csCoulomb, csBabushkin )
        return( cs )
    end


    """
    `JAC.PhotoExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                      settings::PhotoExcitation.Settings; output=true)`  ... to compute the photo-excitation amplitudes 
         and all properties as requested by the given settings. A list of lines::Array{PhotoExcitation.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                           settings::PhotoExcitation.Settings; output=true)
        println("")
        printstyled("JAC.PhotoExcitation.computeLines(): The computation of the excitation cross sections, etc. starts now ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        ##x JAC.define("standard grid", grid)
        lines = JAC.PhotoExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.PhotoExcitation.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = PhotoExcitation.Line[]
        for  line in lines
            newLine = JAC.PhotoExcitation.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        JAC.PhotoExcitation.displayCrossSections(stdout, lines, settings)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary    JAC.PhotoExcitation.displayCrossSections(iostream, lines, settings)   end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `JAC.PhotoExcitation.computeMatrix_obsolete(Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalBasis::Basis, initialBasis::Basis, grid::Radial.Grid,
                                       settings::PhotoExcitation.Settings)`  ... to compute the transition matrix 
         (O^Mp_rs) = (<finalCSF_r|| O^Mp (omega; gauge) ||initialCSF_s>) of the Mp multiplet field for the given transition energy and gauge, 
         and between the CSF_r from the finalBasis and the CSF_s from the initialBasis. A (non-quadratic) matrix::Array{Float64,2} with 
         dimensions [length(finalBasis.csfs) x length(initialBasis.csfs)] is returned. Note that this transition matrix is specific to transitions 
         of the given multipolarity and omega.
    """
    function computeMatrix_obsolete(Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalBasis::Basis, initialBasis::Basis, grid::Radial.Grid, 
                           settings::PhotoExcitation.Settings)
        nf = length(finalBasis.csfs);    ni = length(initialBasis.csfs)
  
        print("Compute radiative T^L matrix of dimension $nf x $ni for the given initial- and final-state bases ...")
        matrix = zeros(Float64, nf, ni)
        for  r = 1:nf
            for  s = 1:ni
                wa = compute("angular coefficients: 1-p, Ratip2013", Mp.L, finalBasis.csfs[r], initialBasis.csfs[s])
                me = 0.
                for  coeff in wa
                    jj = JAC.subshell_2j(finalBasis.orbitals[coeff.a].subshell)
                    me = me + coeff.T * sqrt( jj + 1) * JAC.InteractionStrength.multipole(Mp, gauge, omega, finalBasis.orbitals[coeff.a],  
                                                                                                            initialBasis.orbitals[coeff.b], grid)
                end
                matrix[r,s] = me
            end
        end 
        println("   ... done.")
        return( matrix )
    end


    """
    `JAC.PhotoExcitation.computeStatisticalTensor(k::Int64, q::Int64, line::PhotoExcitation.Line, stokes::ExpStokes)`  
         ... to compute the statistical tensor (component) rho_{k,q} of the final level for the excitation of unpolarized atoms by 
             plane-wave photons, whose polarization is described by the given (experimental) Stokes parameters. 
             A rho_kq::EmProperty is returned.
    """
    function  computeStatisticalTensor(k::Int64, q::Int64, line::PhotoExcitation.Line, stokes::ExpStokes)
        if  k == q  ||  k == -q   rhoCoulomb = -3.;    rhoBabushkin = -3.    else   rhoCoulomb = 0.;    rhoBabushkin = 0.    end
        rho = EmProperty( rhoCoulomb, rhoBabushkin )
        return( rho )
    end


    """
    `JAC.PhotoExcitation.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoExcitation.Settings)`  ... to determine a list 
         of PhotoExcitation.Channel for a transitions from the initial to final level and by taking into account the particular settings of 
         for this computation; an Array{PhotoExcitation.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoExcitation.Settings)
        channels = Radiative.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp in settings.multipoles
            if   JAC.AngularMomentum.isAllowedMultipole(symi, mp, symf)
                hasMagnetic = false
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp)[1] == 'E'  &&   gauge == JAC.UseCoulomb      push!(channels, Radiative.Channel(mp, JAC.Coulomb,   0.) )
                    elseif string(mp)[1] == 'E'  &&   gauge == JAC.UseBabushkin    push!(channels, Radiative.Channel(mp, JAC.Babushkin, 0.) )  
                    elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)               push!(channels, Radiative.Channel(mp, JAC.Magnetic,  0.) );
                                                        hasMagnetic = true; 
                    end 
                end
            end
        end
        ##x println("PhotoExcitation.determineChannels-aa: channels = $channels ")
        return( channels )  
    end


    """
    `JAC.PhotoExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoExcitation.Settings)`  
         ... to determine a list of photo-excitation Line's for transitions between the levels from the given initial- and final-state multiplets 
         and by taking into account the particular selections and settings for this computation; an Array{PhotoExcitation.Line,1} is returned. 
         Apart from the level specification, all physical properties are set to zero during the initialization process.  
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoExcitation.Settings)
        if    settings.selectLines    selectLines   = true;   selectedLines = JAC.determine("selected lines", settings.selectedLines)
        else                          selectLines   = false
        end
    
        lines = PhotoExcitation.Line[]
        for  i = 1:length(initialMultiplet.levels)
            for  f = 1:length(finalMultiplet.levels)
                if  selectLines  &&  !(haskey(selectedLines, (i,f)) )    continue   end
                ##x println("PhotoExcitation.determineLines-aa: angular i = $i, f = $f")
                omega    = abs( finalMultiplet.levels[f].energy - initialMultiplet.levels[i].energy) + settings.photonEnergyShift
                if   omega == 0.  ||  omega < settings.mimimumPhotonEnergy  ||  omega > settings.maximumPhotonEnergy    continue   end  

                channels = JAC.PhotoExcitation.determineChannels(finalMultiplet.levels[f], initialMultiplet.levels[i], settings)
                if   length(channels) == 0   continue   end
                push!( lines, PhotoExcitation.Line(initialMultiplet.levels[i], finalMultiplet.levels[f], omega, EmProperty(0., 0.), 
                                                   TensorComp[], true, channels) )
            end
        end
        return( lines )
    end


    """
    `JAC.PhotoExcitation.displayLines(lines::Array{PhotoExcitation.Line,1})`  ... to display a list of lines and channels that have been 
         selected due to the prior settings. A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{PhotoExcitation.Line,1})
        println(" ")
        println("  Selected photo-excitation lines:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(120))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.flushleft(30, "List of multipoles"; na=4);             sb = sb * JAC.TableStrings.hBlank(34)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(120)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", line.omega)) * "    "
            mpGaugeList = Tuple{JAC.EmMultipole,JAC.EmGauge}[]
            for  i in 1:length(line.channels)
                push!( mpGaugeList, (line.channels[i].multipole, line.channels[i].gauge) )
            end
            ##x println("PhotoExcitation.diplayLines-ad: mpGaugeList = ", mpGaugeList)
            sa = sa * JAC.TableStrings.multipoleGaugeTupels(50, mpGaugeList)
            println( sa )
        end
        println("  ", JAC.TableStrings.hLine(120))
        #
        return( nothing )
    end


    """
    `JAC.PhotoExcitation.displayCrossSections(stream::IO, lines::Array{PhotoExcitation.Line,1}, settings::PhotoExcitation.Settings)`  ... to list all results, 
         energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayCrossSections(stream::IO, lines::Array{PhotoExcitation.Line,1}, settings::PhotoExcitation.Settings)
        println(stream, " ")
        println(stream, "  Photo-excitation cross sections for (completely) linearly-polarized plane-wave photons:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(105))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(14, "Energy"   ; na=3);               
        sb = sb * JAC.TableStrings.center(14,JAC.TableStrings.inUnits("energy"); na=3)
        sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(30, "Cou -- Cross section -- Bab"; na=2);       
        sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("cross section")*"          "*
                                              JAC.TableStrings.inUnits("cross section"); na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(105)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", line.omega)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles) * "          "
            sa = sa * JAC.TableStrings.flushleft(11, mpString[1:10];  na=3)
            sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", line.crossSection.Coulomb))     * "    "
            sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", line.crossSection.Babushkin))   * "    "
            ##x sa = sa * @sprintf("%.6e", line.anisotropy.Coulomb)      * "    "
            ##x sa = sa * @sprintf("%.6e", line.anisotropy.Babushkin)    * "    "
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(105))
        #
        #
        if  settings.calcForStokes
            stokes = settings.stokes
            println(stream, " ")
            println(stream, "  Photo-excitation cross sections for incident plane-wave photons with given $stokes:")
            println(stream, " ")
            println(stream, "  ", JAC.TableStrings.hLine(105))
            sa = "  ";   sb = "  "
            sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
            sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
            sa = sa * JAC.TableStrings.center(14, "Energy"   ; na=3);               
            sb = sb * JAC.TableStrings.center(14,JAC.TableStrings.inUnits("energy"); na=3)
            sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=3);                        sb = sb * JAC.TableStrings.hBlank(14)
            sa = sa * JAC.TableStrings.center(30, "Cou -- Cross section -- Bab"; na=2);       
            sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("cross section")*"          "*
                                                  JAC.TableStrings.inUnits("cross section"); na=2)
            println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(105)) 
            for  line in lines
                sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                               fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
                sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", line.omega)) * "    "
                multipoles = EmMultipole[]
                for  ch in line.channels
                    multipoles = push!( multipoles, ch.multipole)
                end
                multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles) * "          "
                sa = sa * JAC.TableStrings.flushleft(11, mpString[1:10];  na=3)
                crossSection = JAC.PhotoExcitation.computeCrossSection(line, stokes)
                sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", crossSection.Coulomb))     * "    "
                sa = sa * @sprintf("%.6e", JAC.convert("cross section: from atomic", crossSection.Babushkin))   * "    "
                println(stream, sa)
            end
            println(stream, "  ", JAC.TableStrings.hLine(105))
        end
        #
        #
        if  settings.calcTensors
            stokes = settings.stokes
            println(stream, " ")
            println(stream, "  Statistical tensors rho_kq  and alignment parameters A_kq for the excitation by incident plane-wave photons")
            println(stream, "  with given $stokes:")
            println(stream, " ")
            println(stream, "  ", JAC.TableStrings.hLine(145))
            sa = "  ";   sb = "  "
            sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
            sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
            sa = sa * JAC.TableStrings.center(14, "Energy"   ; na=3);               
            sb = sb * JAC.TableStrings.center(14,JAC.TableStrings.inUnits("energy"); na=3)
            sa = sa * JAC.TableStrings.center(10, "Multipoles"; na=1);                        sb = sb * JAC.TableStrings.hBlank(12)
            sa = sa * JAC.TableStrings.center(12, "     k    q    "; na=0);                   sb = sb * JAC.TableStrings.hBlank(14)
            sa = sa * JAC.TableStrings.center(26, "Cou -- rho_kq -- Bab"; na=5);              sb = sb * JAC.TableStrings.hBlank(14)   
            sa = sa * JAC.TableStrings.center(26, "Cou --  A_kq  -- Bab"; na=2);              sb = sb * JAC.TableStrings.hBlank(14)   
            println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(145)) 
            for  line in lines
                sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                               fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
                sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", line.omega)) * "    "
                multipoles = EmMultipole[]
                for  ch in line.channels
                    multipoles = push!( multipoles, ch.multipole)
                end
                multipoles = unique(multipoles);   mpString = JAC.TableStrings.multipoleList(multipoles) * "          "
                sa = sa * JAC.TableStrings.flushleft(11, mpString[1:10];  na=3)
                println(stream, sa);   sc = JAC.TableStrings.hBlank( length(sa)-1 )
                tc00   = JAC.PhotoExcitation.computeStatisticalTensor(0, 0, line, stokes)
                for  k = 2:2
                    for  q = -k:k
                        tc = JAC.PhotoExcitation.computeStatisticalTensor(k, q, line, stokes)
                        if  abs(tc.Coulomb) == abs(tc.Babushkin) == 0.    continue    end
                        #
                        aCoulomb = tc.Coulomb / tc00.Coulomb;    aBabushkin = tc.Babushkin / tc00.Babushkin
                        sd = sc * JAC.TableStrings.level(k) * JAC.TableStrings.level(q) * "   "
                        sd = sd * @sprintf("%.6e", tc.Coulomb)     * "  "
                        sd = sd * @sprintf("%.6e", tc.Babushkin)   * "    "
                        sd = sd * @sprintf("%.6e", aCoulomb)       * "  "
                        sd = sd * @sprintf("%.6e", aBabushkin)     * "    "
                        println(stream, sd)
                    end
                end
            end
            println(stream, "  ", JAC.TableStrings.hLine(145))
        end
        #
        #
        if  settings.calcTensors   &&   settings.calcPhotonDm
            stokes = settings.stokes
            println(stream, " ")
            println(stream, "  Reduced density matrix of the fluorescence photons at selected solid angles for the prior excitation by incident plane-wave photons")
            println(stream, "  with given $stokes:")
            println(stream, " ")
            println(stream, "  Not yet implemented !!")
        end
        #
        return( nothing )
    end

end # module

