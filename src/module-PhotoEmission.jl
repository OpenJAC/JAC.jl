
"""
`module  JAC.PhotoEmission`  ... a submodel of JAC that contains all methods for computing Einstein coefficients, oscillator strength, etc. between 
                             some initial and final-state multiplets; it is using JAC, JAC.ManyElectron, JAC.Radial.
"""
module PhotoEmission

    using Printf, JAC, JAC.BasicTypes, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


    """
    `struct  PhotoEmission.Settings`  ... defines a type for the details and parameters of computing radiative lines.

        + multipoles              ::Array{EmMultipoles}     ... Specifies the (radiat. field) multipoles to be included.
        + gauges                  ::Array{UseGauge}         ... Gauges to be included into the computations.
        + calcAnisotropy          ::Bool                    ... True, if the anisotropy (structure) functions are to be 
                                                                calculated and false otherwise 
        + printBeforeComputation  ::Bool                    ... True, if all energies and lines are printed before comput.
        + selectLines             ::Bool                    ... True, if lines are selected individually for the computat.
        + selectedLines           ::Array{Tuple{Int64,Int64},1}  ... List of lines as tupels (inital-level, final-level).
        + photonEnergyShift       ::Float64                 ... An overall energy shift for all photon energies.
        + mimimumPhotonEnergy     ::Float64                 ... minimum transition energy for which (photon) transitions 
                                                                are included into the computation.
        + maximumPhotonEnergy     ::Float64                 ... maximum transition energy for which (photon) transitions 
                                                                are included.
    """
    struct Settings 
        multipoles                ::Array{EmMultipole,1}
        gauges                    ::Array{UseGauge}
        calcAnisotropy            ::Bool         
        printBeforeComputation    ::Bool
        selectLines               ::Bool
        selectedLines             ::Array{Tuple{Int64,Int64},1}
        photonEnergyShift         ::Float64
        mimimumPhotonEnergy       ::Float64   
        maximumPhotonEnergy       ::Float64     
    end 


    """
    `JAC.PhotoEmission.Settings()`  ... constructor for the default values of radiative line computations
    """
    function Settings()
        Settings(EmMultipole[E1], UseGauge[BasicTypes.UseCoulomb], false, false, false, Array{Tuple{Int64,Int64},1}[], 0., 0., 0.)
    end


    # `Base.show(io::IO, settings::PhotoEmission.Settings)`  ... prepares a proper printout of the variable settings::PhotoEmissionSettings.
    function Base.show(io::IO, settings::PhotoEmission.Settings) 
        println(io, "PhotoEmission.Settings(with multipoles = $(settings.multipoles), use-gauges = $(settings.gauges), "   *
                    "calcAnisotropy = $(settings.calcAnisotropy), " )
        println(io, "                 printBeforeComputation = $(settings.printBeforeComputation), selectLines = $(settings.selectLines), " *
                    "selectedLines = $(settings.selectedLines), photonEnergyShift = $(settings.photonEnergyShift) " )
        println(io, "                 mimimumPhotonEnergy = $(settings.mimimumPhotonEnergy), " *
                    "maximumPhotonEnergy = $(settings.maximumPhotonEnergy)) ") 
    end


    """
    `struct  Channel`  ... defines a type for a single radiative emission/absorption channel that specifies the multipole, gauge and amplitude.

        + multipole         ::EmMultipole        ... Multipole of the photon emission/absorption.
        + gauge             ::EmGauge            ... Gauge for dealing with the (coupled) radiation field.
        + amplitude         ::Complex{Float64}   ... Amplitude of this multiple channel.
    """
    struct Channel 
        multipole           ::EmMultipole
        gauge               ::EmGauge
        amplitude           ::Complex{Float64}
    end 


    # `Base.show(io::IO, channel::PhotoEmission.Channel)`  ... prepares a proper printout of the variable channel::PhotoEmission.Channel.
    function Base.show(io::IO, channel::PhotoEmission.Channel) 
        print(io, "PhotoEmission.Channel($(channel.multipole), $(channel.gauge), amp = $(channel.amplitude)) ") 
    end


    """
    `struct  Line`  ... defines a type for a radiative line that may include the definition of sublines and their corresponding amplitudes.

        + initialLevel   ::Level               ... initial-(state) level
        + finalLevel     ::Level               ... final-(state) level
        + omega          ::Float64             ... Transition frequency of this line; can be shifted w.r.t. the level energies.
        + photonRate     ::EmProperty          ... Total rate of this line.
        + angularBeta    ::EmProperty          ... Angular beta_2 coefficient.
        + hasSublines    ::Bool                ... Determines whether the sublines are defined in terms of their multipolarity, amplitude, or not.
        + channels       ::Array{PhotoEmission.Channel,1}  ... List of radiative (photon) channels
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        omega            ::Float64
        photonRate       ::EmProperty
        angularBeta      ::EmProperty
        hasSublines      ::Bool
        channels         ::Array{PhotoEmission.Channel,1}
    end 


    """
    `JAC.PhotoEmission.Line(initialLevel::Level, finalLevel::Level, photonRate::Float64)`  ... constructor an radiative line between a
                                                                                           specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, omega::Float64, photonRate::EmProperty)
       Line(initialLevel, finalLevel, omega, photonRate, EmProperty(0., 0.), false, 0, PhotoEmission.Channel[])
    end


    # `Base.show(io::IO, line::PhotoEmissionLine)`  ... prepares a proper printout of the variable line::PhotoEmission.Line.
    function Base.show(io::IO, line::PhotoEmission.Line) 
        println(io, "initialLevel:         $(line.initialLevel)  ")
        println(io, "finalLevel:           $(line.finalLevel)  ")
        println(io, "omega:                $(line.omega)  ")
        println(io, "photonRate:           $(line.photonRate)  ")
        println(io, "angularBeta:          $(line.angularBeta)  ")
        println(io, "hasSublines:          $(line.hasSublines)  ")
        println(io, "channels:             $(line.channels)  ")
    end


    """
    `JAC.PhotoEmission.amplitude(kind::String, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                             grid::Radial.Grid; display::Bool=false)`  ... to compute the kind = (absorption or emission) amplitude
                             <alpha_f J_f || O^(Mp, kind) || alpha_i J_i> for the interaction with  photon of multipolarity Mp and
                            for the given transition energy and gauge. A value::ComplexF64 is returned. The amplitude value is printed
                            to screen if display=true.
    """
    function amplitude(kind::String, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                       grid::Radial.Grid; display::Bool=false, printout::Bool=true)
        
        if      kind == "emission"
        #-------------------------
            nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
            if  printout   printstyled("Compute radiative $(Mp) matrix of dimension $nf x $ni in the initial- and final-state bases " *
                                       "for the transition [$(initialLevel.index)-$(finalLevel.index)] ... ", color=:light_green)    end
            matrix = zeros(ComplexF64, nf, ni)
            #
            for  r = 1:nf
                if  finalLevel.basis.csfs[r].J != finalLevel.J          ||  finalLevel.basis.csfs[r].parity   != finalLevel.parity    continue    end 
                for  s = 1:ni
                    if  initialLevel.basis.csfs[s].J != initialLevel.J  ||  initialLevel.basis.csfs[s].parity != initialLevel.parity  continue    end 
                    wa = JAC.compute("angular coefficients: 1-p, Grasp92", 0, Mp.L, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                    me = 0.
                    for  coeff in wa
                        ja = JAC.subshell_2j(finalLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(initialLevel.basis.orbitals[coeff.b].subshell)
                        MbaCheng  = JAC.InteractionStrength.MbaEmissionCheng(Mp, gauge, omega, finalLevel.basis.orbitals[coeff.a],  
                                                                                                 initialLevel.basis.orbitals[coeff.b], grid)
                        me = me + coeff.T * MbaCheng  ## sqrt( (ja + 1)/(jb + 1) ) * ##   0.707106781186548 *
                    end
                    ##x println("r = $r, s = $s, me = $me")
                    matrix[r,s] = me
                end
            end 
            if  printout   printstyled("done. \n", color=:light_green)    end
            amplitude = transpose(finalLevel.mc) * matrix * initialLevel.mc 
            #
            #
        elseif  kind == "absorption"
        #---------------------------
            iLevel = finalLevel;   fLevel = initialLevel
            amplitude = JAC.PhotoEmission.amplitude("emission", Mp, gauge, omega, fLevel, iLevel, grid) 
            amplitude = (amplitude)
        else    error("stop a")
        end
        
        if  display  
            println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                    " O^($Mp, $kind) ($omega a.u., $gauge) ||" *
                    " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
        end
        
        return( amplitude )
    end


    """
    `JAC.PhotoEmission.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::PhotoEmission.Settings; 
                                output=true)`  ... to compute the radiative transition amplitudes and all properties as requested by the 
         given settings. A list of lines::Array{PhotoEmission.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::PhotoEmission.Settings; output=true) 
        # Define a common subshell list for both multiplets
        subshellList = JAC.generate("subshells: ordered list for two bases", finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
        JAC.define("relativistic subshell list", subshellList; printout=true)
        ##x JAC.define("standard grid", grid)
        println("")
        printstyled("JAC.PhotoEmission.computeLines(): The computation of the transition amplitudes and properties starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        ##x println("basis = $(finalMultiplet.levels[1].basis) ")
        lines = JAC.PhotoEmission.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBeforeComputation    JAC.PhotoEmission.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = PhotoEmission.Line[]
        for  line in lines
            newLine = JAC.PhotoEmission.computeAmplitudesProperties(line, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        JAC.PhotoEmission.displayRates(stdout, newLines)
        if  settings.calcAnisotropy    JAC.PhotoEmission.displayAnisotropies(stdout, newLines)    end
        JAC.PhotoEmission.displayLifetimes(stdout, newLines)
        #
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary   JAC.PhotoEmission.displayRates(iostream, newLines)       
                           JAC.PhotoEmission.displayAnisotropies(stdout, newLines)
                           JAC.PhotoEmission.displayLifetimes(iostream, newLines)
        end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `JAC.PhotoEmission.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                       settings::PhotoEmission.Settings; output::Bool=true, printout::Bool=true)`  
        ... to compute the radiative transition amplitudes and all properties as requested by the given settings. The computations
            and printout is adapted for larger cascade computations by including only lines with at least one channel and by sending
            all printout to a summary file only. A list of lines::Array{PhotoEmission.Lines} is returned.
    """
    function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                  settings::PhotoEmission.Settings; output=true, printout::Bool=true) 
        # Define a common subshell list for both multiplets
        subshellList = JAC.generate("subshells: ordered list for two bases", finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
        JAC.define("relativistic subshell list", subshellList; printout=false)
        lines = JAC.PhotoEmission.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        # if  settings.printBeforeComputation    JAC.PhotoEmission.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = PhotoEmission.Line[]
        for  line in lines
            newLine = JAC.PhotoEmission.computeAmplitudesProperties(line, grid, settings, printout=printout) 
            #
            # Don't add this line if it does not contribute to the decay
            wa = 0.
            for  ch in newLine.channels   wa = wa + abs(ch.amplitude)^2    end
            if   wa == 0.    continue    end
            push!( newLines, newLine)
        end
        # Print all results to a summary file, if requested
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary   JAC.PhotoEmission.displayRates(iostream, newLines)    end
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end



    """
    `JAC.PhotoEmission.computeAmplitudesProperties(line::PhotoEmission.Line, grid::Radial.Grid, settings::Einstein.Settings; printout::Bool=true)`  
        ... to compute all amplitudes and properties of the given line; a line::Einstein.Line is returned for which the amplitudes and 
            properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::PhotoEmission.Line, grid::Radial.Grid, settings::PhotoEmission.Settings; printout::Bool=true)
        global JAC_counter
        newChannels = PhotoEmission.Channel[];    rateC = 0.;    rateB = 0.
        for channel in line.channels
            #
            amplitude = JAC.PhotoEmission.amplitude("emission", channel.multipole, channel.gauge, line.omega, 
                                                line.finalLevel, line.initialLevel, grid, printout=printout)
            #
            push!( newChannels, PhotoEmission.Channel( channel.multipole, channel.gauge, amplitude) )
            if       channel.gauge == JAC.Coulomb     rateC = rateC + abs(amplitude)^2
            elseif   channel.gauge == JAC.Babushkin   rateB = rateB + abs(amplitude)^2
            elseif   channel.gauge == JAC.Magnetic    rateB = rateB + abs(amplitude)^2;   rateC = rateC + abs(amplitude)^2
            end
        end
        #     
        # Calculate the photonrate and angular beta if requested 
        wa = 8.0pi * JAC.give("alpha") * line.omega / (JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) * 
                                                      (JAC.AngularMomentum.twoJ(line.finalLevel.J) + 1)
        photonrate  = EmProperty(wa * rateC, wa * rateB)    
        angularBeta = EmProperty(-9., -9.)
        line = PhotoEmission.Line( line.initialLevel, line.finalLevel, line.omega, photonrate, angularBeta, true, newChannels)
        return( line )
    end


    """
    `JAC.PhotoEmission.computeMatrix_obsolete(Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid,
                                 settings::PhotoEmission.Settings)`  ... to compute the transition matrix 
         (O^Mp_rs) = (<finalCSF_r|| O^Mp (omega; gauge) ||initialCSF_s>)   of the Mp multiplet field for the given transition energy and gauge, 
         and between the CSF_r from the basis of the finalLevel and the CSF_s from the basis of the initialLevel. A (non-quadratic) 
         matrix::Array{Float64,2} with dimensions [length(finalBasis.csfs) x length(initialBasis.csfs)] is returned. Note that this transition 
         matrix is specific to a particular transitions with the given multipolarity and omega.
    """
    function computeMatrix_obsolete(Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid, 
                           settings::PhotoEmission.Settings)
        nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
  
        print("Compute radiative $(Mp) matrix of dimension $nf x $ni in the initial- and final-state bases for the transition " *
              "[$(initialLevel.index)-$(finalLevel.index)] ... ")
        matrix = zeros(ComplexF64, nf, ni)
        for  r = 1:nf
            if  finalLevel.basis.csfs[r].J != finalLevel.J          ||  finalLevel.basis.csfs[r].parity   != finalLevel.parity    continue    end 
            for  s = 1:ni
                if  initialLevel.basis.csfs[s].J != initialLevel.J  ||  initialLevel.basis.csfs[s].parity != initialLevel.parity  continue    end 
                ##x wb = compute("angular coefficients: 1-p, Ratip2013",  Mp.L, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                wa = compute("angular coefficients: 1-p, Grasp92", 0, Mp.L, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                me = 0.
                for  coeff in wa
                    ja = JAC.subshell_2j(finalLevel.basis.orbitals[coeff.a].subshell)
                    jb = JAC.subshell_2j(initialLevel.basis.orbitals[coeff.b].subshell)
                    MbaCheng  = JAC.InteractionStrength.MbaEmissionCheng(Mp, gauge, omega, finalLevel.basis.orbitals[coeff.a],  
                                                                                   initialLevel.basis.orbitals[coeff.b], grid)
                    MbaAndrey = JAC.InteractionStrength.MbaEmissionAndrey(Mp, gauge, omega, finalLevel.basis.orbitals[coeff.a],  
                                                                                    initialLevel.basis.orbitals[coeff.b], grid)
                    ##x println("M-Cheng = $MbaCheng,  M-Andrey = $MbaAndrey,  a = $(coeff.a),  b = $(coeff.b),  T = $(coeff.T)")
                    me = me + coeff.T *   ## sqrt( (ja + 1)/(jb + 1) ) * ##   0.707106781186548 *
                              JAC.InteractionStrength.MbaEmissionCheng(Mp, gauge, omega, finalLevel.basis.orbitals[coeff.a],  
                                                                                 initialLevel.basis.orbitals[coeff.b], grid)
                end
                matrix[r,s] = me
            end
        end 
        println("done.")
        return( matrix )
    end


    """
    `JAC.PhotoEmission.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoEmission.Settings)`  ... to determine a list of 
         PhotoEmission.Channel for a transitions from the initial to final level and by taking into account the particular settings of for this 
         computation; an Array{PhotoEmission.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoEmission.Settings)
        channels = PhotoEmission.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp in settings.multipoles
            if   JAC.AngularMomentum.isAllowedMultipole(symi, mp, symf)
                hasMagnetic = false
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp)[1] == 'E'  &&   gauge == BasicTypes.UseCoulomb      push!(channels, PhotoEmission.Channel(mp, JAC.Coulomb,   0.) )
                    elseif string(mp)[1] == 'E'  &&   gauge == BasicTypes.UseBabushkin    push!(channels, PhotoEmission.Channel(mp, JAC.Babushkin, 0.) )  
                    elseif string(mp)[1] == 'M'  &&   !(hasMagnetic)                      push!(channels, PhotoEmission.Channel(mp, JAC.Magnetic,  0.) );
                                                        hasMagnetic = true; 
                    end 
                end
            end
        end
        return( channels )  
    end


    """
    `JAC.PhotoEmission.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoEmission.Settings)`  ... to determine a 
         list of PhotoEmission Line's for transitions between the levels from the given initial- and final-state multiplets and by taking into account 
         the particular selections and settings for this computation; an Array{PhotoEmission.Line,1} is returned. Apart from the level specification,
         all physical properties are set to zero during the initialization process.  
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoEmission.Settings)
        if    settings.selectLines    selectLines   = true
            selectedLines = JAC.determineSelectedLines(settings.selectedLines, initialMultiplet, finalMultiplet)
        else                          selectLines   = false
        end
    
        lines = PhotoEmission.Line[]
        for  i = 1:length(initialMultiplet.levels)
            for  f = 1:length(finalMultiplet.levels)
                if  selectLines  &&  !((i,f) in selectedLines )    continue   end
                omega = initialMultiplet.levels[i].energy - finalMultiplet.levels[f].energy   + settings.photonEnergyShift
                if  omega <= settings.mimimumPhotonEnergy  ||  omega > settings.maximumPhotonEnergy    continue   end  

                channels = JAC.PhotoEmission.determineChannels(finalMultiplet.levels[f], initialMultiplet.levels[i], settings) 
                if   length(channels) == 0   continue   end
                push!( lines, PhotoEmission.Line(initialMultiplet.levels[i], finalMultiplet.levels[f], omega, EmProperty(0., 0.), EmProperty(0., 0.),
                                             true, channels) )
            end
        end
        return( lines )
    end


    """
    `JAC.PhotoEmission.displayAnisotropies(stream::IO, lines::Array{PhotoEmission.Line,1})`  ... to list all energies and anisotropy parameters of 
         the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayAnisotropies(stream::IO, lines::Array{PhotoEmission.Line,1})
        println(stream, " ")
        println(stream, "  Anisotropy (structure) functions:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(153))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(14, "Energy"   ; na=4);               
        sb = sb * JAC.TableStrings.center(14,JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.flushleft(20, "Multipoles";   na=0);              
        sa = sa * JAC.TableStrings.center(14, "f_2 (Coulomb)";   na=3);       
        sa = sa * JAC.TableStrings.center(14, "f_2 (Babushkin)"; na=3);       
        sa = sa * JAC.TableStrings.center(14, "f_4 (Coulomb)";   na=3);       
        sa = sa * JAC.TableStrings.center(14, "f_4 (Babushkin)"; na=3);       
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(153)) 
        #   
        for  line in lines
            f2Coulomb   = 0.0im;    f2Babushkin   = 0.0im;    f4Coulomb = 0.0im;    f4Babushkin = 0.0im;   mpList = EmMultipole[]
            normCoulomb = 0.;       normBabushkin = 0.
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", line.omega)) * "    "
            # Now compute the anisotropy parameter from the amplitudes
            for  ch in line.channels
                if  !(ch.multipole  in  mpList)    push!(mpList, ch.multipole)    end
                for  chp  in line.channels
                    #
                    if  ch.gauge  in  [ EmGauge("Coulomb"), EmGauge("Magnetic")]
                        if  ch == chp    normCoulomb = normCoulomb + (abs(ch.amplitude)^2)    end
                        angL  = AngularJ64(ch.multipole.L);   p  = JAC.multipole_p(ch.multipole)
                        angLp = AngularJ64(ch.multipole.L);   pp = JAC.multipole_p(chp.multipole)
                        angJi = line.initialLevel.J;    angJf = line.finalLevel.J
                        #
                        f2Coulomb = f2Coulomb + 1.0im^(chp.multipole.L + pp - ch.multipole.L - p + 16)                         *
                                                JAC.AngularMomentum.phaseFactor([angJf, +1, angJi, +1, AngularJ64(3)])         *
                                                sqrt( (2ch.multipole.L+1) * (2chp.multipole.L+1) )                             * 
                                    JAC.AngularMomentum.ClebschGordan(angL, AngularM64(1), angLp, AngularM64(-1), AngularJ64(2), AngularM64(0) ) *
                                                (1. + (-1.)^(ch.multipole.L + p + chp.multipole.L + pp - 2) )                  *
                                                JAC.AngularMomentum.Wigner_6j(angL, angLp, AngularJ64(2), angJi, angJi, angJf) *
                                                conj(ch.amplitude) * chp.amplitude
                        #
                        f4Coulomb = f4Coulomb + 1.0im^(chp.multipole.L + pp - ch.multipole.L - p + 16)                         *
                                                JAC.AngularMomentum.phaseFactor([angJf, +1, angJi, +1, AngularJ64(5)])         *
                                                sqrt( (2ch.multipole.L+1) * (2chp.multipole.L+1) )                             * 
                                    JAC.AngularMomentum.ClebschGordan(angL, AngularM64(1), angLp, AngularM64(-1), AngularJ64(4), AngularM64(0) ) *
                                                (1. + (-1.)^(ch.multipole.L + p + chp.multipole.L + pp - 2) )                  *
                                                JAC.AngularMomentum.Wigner_6j(angL, angLp, AngularJ64(4), angJi, angJi, angJf) *
                                                conj(ch.amplitude) * chp.amplitude
                    end
                    #
                    if  ch.gauge  in  [ EmGauge("Babushkin"), EmGauge("Magnetic")]
                        if  ch == chp    normBabushkin = normBabushkin + (abs(ch.amplitude)^2)    end
                        angL  = AngularJ64(ch.multipole.L);   p  = JAC.multipole_p(ch.multipole)
                        angLp = AngularJ64(ch.multipole.L);   pp = JAC.multipole_p(chp.multipole)
                        angJi = line.initialLevel.J;    angJf = line.finalLevel.J
                        #
                        f2Babushkin = f2Babushkin + 1.0im^(chp.multipole.L + pp - ch.multipole.L - p + 16)                         *
                                                    JAC.AngularMomentum.phaseFactor([angJf, +1, angJi, +1, AngularJ64(3)])         *
                                                    sqrt( (2ch.multipole.L+1) * (2chp.multipole.L+1) )                             * 
                                      JAC.AngularMomentum.ClebschGordan(angL, AngularM64(1), angLp, AngularM64(-1), AngularJ64(2), AngularM64(0) ) *
                                                    (1. + (-1.)^(ch.multipole.L + p + chp.multipole.L + pp - 2) )                  *
                                                    JAC.AngularMomentum.Wigner_6j(angL, angLp, AngularJ64(2), angJi, angJi, angJf) *
                                                    conj(ch.amplitude) * chp.amplitude
                        #
                        f4Babushkin = f4Babushkin + 1.0im^(chp.multipole.L + pp - ch.multipole.L - p + 16)                         *
                                                    JAC.AngularMomentum.phaseFactor([angJf, +1, angJi, +1, AngularJ64(5)])         *
                                                    sqrt( (2ch.multipole.L+1) * (2chp.multipole.L+1) )                             * 
                                      JAC.AngularMomentum.ClebschGordan(angL, AngularM64(1), angLp, AngularM64(-1), AngularJ64(4), AngularM64(0) ) *
                                                    (1. + (-1.)^(ch.multipole.L + p + chp.multipole.L + pp - 2) )                  *
                                                    JAC.AngularMomentum.Wigner_6j(angL, angLp, AngularJ64(4), angJi, angJi, angJf) *
                                                    conj(ch.amplitude) * chp.amplitude
                    end
                end
            end
            f2Coulomb   = f2Coulomb   / normCoulomb   * sqrt(JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) / 2.
            f4Coulomb   = f4Coulomb   / normCoulomb   * sqrt(JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) / 2.
            f2Babushkin = f2Babushkin / normBabushkin * sqrt(JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) / 2.
            f4Babushkin = f4Babushkin / normBabushkin * sqrt(JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) / 2.
            #
            sa = sa * JAC.TableStrings.flushleft( 17, JAC.TableStrings.multipoleList(mpList);  na=3)
            sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", f2Coulomb.re);    na=3) 
            sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", f2Babushkin.re);  na=3)
            sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", f4Coulomb.re);    na=3)
            sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", f4Babushkin.re);  na=3)
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(153))
        #
        return( nothing )
    end


    """
    `JAC.PhotoEmission.displayLifetimes(stream::IO, lines::Array{PhotoEmission.Line,1})`  ... to list all lifetimes as derived from the selected
         lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayLifetimes(stream::IO, lines::Array{PhotoEmission.Line,1})
        # Determine all initial levels (and their level information) for printing the lifetimes
        ilevels = Int64[];   istr = String[]
        for  i = 1:length(lines)
            ii = lines[i].initialLevel.index;   
            if  !(ii in ilevels)   
               sa  = "  ";    sym = LevelSymmetry( lines[i].initialLevel.J, lines[i].initialLevel.parity )
               sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(lines[i].initialLevel.index); na=2)
               sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
               sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", lines[i].initialLevel.energy)) * "    "
               push!( ilevels, ii);    push!( istr, sa )   
            end
        end
        # Determine the lifetime (in a.u.) of the selected initial levels
        irates = JAC.EmProperty[]
        for  ii in  ilevels
            waCoulomb = 0.;    waBabushkin = 0.
            for  i = 1:length(lines)
                if   lines[i].initialLevel.index == ii    
                    waCoulomb   = waCoulomb   + lines[i].photonRate.Coulomb
                    waBabushkin = waBabushkin + lines[i].photonRate.Babushkin
                end
            end
            push!(irates, EmProperty( waCoulomb, waBabushkin) )
        end
        
        println(stream, " ")
        println(stream, "  PhotoEmission lifetimes (as derived from these computations):")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(111))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                              sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                              sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Level energy"   ; na=4);               
        sb = sb * JAC.TableStrings.center(14,JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(12, "Used Gauge"    ; na=7);                     sb = sb * JAC.TableStrings.hBlank(18)
        sa = sa * JAC.TableStrings.center(28, "Lifetime"; na=5);       
        sb = sb * JAC.TableStrings.center(28, "[a.u.]"*"          "*JAC.TableStrings.inUnits("time"); na=5)
        sa = sa * JAC.TableStrings.center(14, "Decay widths"; na=4);       
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(111)) 
        #    
        for  ii = 1:length(ilevels)
            sa = istr[ii]
            sa = sa * "Coulomb          " * @sprintf("%.8e",              1.0/irates[ii].Coulomb)     * "  "
            sa = sa * @sprintf("%.8e", JAC.convert("time: from atomic",   1.0/irates[ii].Coulomb) )   * "    "
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic",     irates[ii].Coulomb) )
            println(stream, sa)
            sa = repeat(" ", length(istr[ii]) )
            sa = sa * "Babushkin        " * @sprintf("%.8e",              1.0/irates[ii].Babushkin)   * "  "
            sa = sa * @sprintf("%.8e", JAC.convert("time: from atomic",   1.0/irates[ii].Babushkin) ) * "    "
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic",     irates[ii].Babushkin) )
            println(stream, sa)
        end
        println(stream, "  ", JAC.TableStrings.hLine(111))
        #
        return( nothing )
    end


    """
    `JAC.PhotoEmission.displayLines(lines::Array{PhotoEmission.Line,1})`  ... to display a list of lines and channels that have been selected due to
         the prior settings. A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{PhotoEmission.Line,1})
        println(" ")
        println("  Selected radiative lines:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(95))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.flushleft(30, "List of multipoles"; na=4);             sb = sb * JAC.TableStrings.hBlank(34)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(95)) 
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
            sa = sa * JAC.TableStrings.multipoleGaugeTupels(50, mpGaugeList)
            println( sa )
        end
        println("  ", JAC.TableStrings.hLine(95))
        println(" ")
        #
        return( nothing )
    end


    """
    `JAC.PhotoEmission.displayRates(stream::IO, lines::Array{PhotoEmission.Line,1})`  ... to list all results, energies, rates, etc. of the selected lines.
         A neat table is printed but nothing is returned otherwise.
    """
    function  displayRates(stream::IO, lines::Array{PhotoEmission.Line,1})
        println(stream, " ")
        println(stream, "  Einstein coefficients, transition rates and oscillator strengths:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(155))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * JAC.TableStrings.hBlank(20)
        sa = sa * JAC.TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * JAC.TableStrings.hBlank(22)
        sa = sa * JAC.TableStrings.center(14, "Energy"   ; na=4);               
        sb = sb * JAC.TableStrings.center(14,JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center( 9, "Multipole"; na=0);                         sb = sb * JAC.TableStrings.hBlank(10)
        sa = sa * JAC.TableStrings.center(11, "Gauge"    ; na=5);                         sb = sb * JAC.TableStrings.hBlank(16)
        sa = sa * JAC.TableStrings.center(30, "A--Einstein--B"; na=3);       
        sb = sb * JAC.TableStrings.center(30, JAC.TableStrings.inUnits("rate")*"          "*JAC.TableStrings.inUnits("rate"); na=2)
        sa = sa * JAC.TableStrings.center(13, "Osc. strength"    ; na=4);                 sb = sb * JAC.TableStrings.hBlank(19)
        sa = sa * JAC.TableStrings.center(14, "Decay widths"; na=4);       
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(155)) 
        #   
        for  line in lines
            for  ch in line.channels
                sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                               fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * JAC.TableStrings.center(18, JAC.TableStrings.symmetries_if(isym, fsym); na=4)
                sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", line.omega)) * "    "
                sa = sa * JAC.TableStrings.center(9,  string(ch.multipole); na=4)
                sa = sa * JAC.TableStrings.flushleft(11, string(ch.gauge);  na=2)
                chRate =  8pi * JAC.give("alpha") * line.omega / (JAC.AngularMomentum.twoJ(line.initialLevel.J) + 1) * (abs(ch.amplitude)^2) * 
                                                                 (JAC.AngularMomentum.twoJ(line.finalLevel.J) + 1)
                sa = sa * @sprintf("%.8e", JAC.recast("rate: radiative, to Einstein A",  line, chRate)) * "  "
                sa = sa * @sprintf("%.8e", JAC.recast("rate: radiative, to Einstein B",  line, chRate)) * "    "
                sa = sa * @sprintf("%.8e", JAC.recast("rate: radiative, to g_f",         line, chRate)) * "    "
                sa = sa * @sprintf("%.8e", JAC.recast("rate: radiative, to decay width", line, chRate)) * "      "
                ##x sa = sa * @sprintf("%.5e  %.5e", ch.amplitude.re, ch.amplitude.im) * "        "
                println(stream, sa)
            end
        end
        println(stream, "  ", JAC.TableStrings.hLine(155))
        #
        return( nothing )
    end

end # module

