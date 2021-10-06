
"""
`module  JAC.AutoIonization`  
    ... a submodel of JAC that contains all methods for computing Auger properties between some initial and final-state 
        multiplets.
"""
module AutoIonization

    using  Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, 
                   ..PlasmaShift, ..Radial, ..SpinAngular, ..TableStrings
    
    """
    `struct  Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing Auger lines.

        + calcAnisotropy      ::Bool               ... True, if the intrinsic alpha_2,4 angular parameters are to be 
                                                       calculated, and false otherwise.
        + printBefore         ::Bool               ... True, if all energies and lines are printed before their evaluation.
        + lineSelection       ::LineSelection      ... Specifies the selected levels, if any.
        + minAugerEnergy      ::Float64            ... Minimum energy of free (Auger) electrons to be included.
        + maxAugerEnergy      ::Float64            ... Maximum energy of free (Auger) electrons to be included.
        + maxKappa            ::Int64              ... Maximum kappa value of partial waves to be included.
        + operator            ::AbstractEeInteraction   
            ... Auger operator that is to be used for evaluating the Auger amplitudes; allowed values are: 
                CoulombInteraction(), BreitInteraction(), ...
    """
    struct Settings  <:  AbstractProcessSettings
        calcAnisotropy        ::Bool         
        printBefore           ::Bool 
        lineSelection         ::LineSelection 
        minAugerEnergy        ::Float64
        maxAugerEnergy        ::Float64
        maxKappa              ::Int64
        operator              ::AbstractEeInteraction 
    end 


    """
    `AutoIonization.Settings()`  ... constructor for the default values of AutoIonization line computations
    """
    function Settings()
        Settings(false, false, LineSelection(), 0., 10e5, 100, CoulombInteraction())
    end


    """
    `AutoIonization.Settings(set::AutoIonization.Settings;`
    
            calcAnisotropy=..,      printBefore=..,              
            minAugerEnergy=..,      maxAugerEnergy=..,          maxKappa=..,            operator=..)
                        
        ... constructor for modifying the given AutoIonization.Settings by 'overwriting' the previously selected parameters.
    """
    function Settings(set::AutoIonization.Settings;    
        calcAnisotropy::Union{Nothing,Bool}=nothing,            printBefore::Union{Nothing,Bool}=nothing, 
        lineSelection::Union{Nothing,LineSelection}=nothing, 
        minAugerEnergy::Union{Nothing,Float64}=nothing,         maxAugerEnergy::Union{Nothing,Float64}=nothing,
        maxKappa::Union{Nothing,Int64}=nothing,                 operator::Union{Nothing,String}=nothing)  
        
        if  calcAnisotropy  == nothing   calcAnisotropyx  = set.calcAnisotropy    else  calcAnisotropyx  = calcAnisotropy   end 
        if  printBefore     == nothing   printBeforex     = set.printBefore       else  printBeforex     = printBefore      end 
        if  lineSelection   == nothing   lineSelectionx   = set.lineSelection     else  lineSelectionx   = lineSelection    end 
        if  minAugerEnergy  == nothing   minAugerEnergyx  = set.minAugerEnergy    else  minAugerEnergyx  = minAugerEnergy   end 
        if  maxAugerEnergy  == nothing   maxAugerEnergyx  = set.maxAugerEnergy    else  maxAugerEnergyx  = maxAugerEnergy   end 
        if  maxKappa        == nothing   maxKappax        = set.maxKappa          else  maxKappax        = maxKappa         end 
        if  operator        == nothing   operatorx        = set.operator          else  operatorx        = operator         end 

        Settings( calcAnisotropyx, printBeforex, lineSelectionx, minAugerEnergyx, maxAugerEnergyx, maxKappax, operatorx)
    end


    # `Base.show(io::IO, settings::AutoIonization.Settings)`  ... prepares a proper printout of the variable settings::AutoIonization.Settings.
    function Base.show(io::IO, settings::AutoIonization.Settings) 
        println(io, "calcAnisotropy:                $(settings.calcAnisotropy)  ")
        println(io, "printBefore:                   $(settings.printBefore)  ")
        println(io, "lineSelection:                 $(settings.lineSelection)  ")
        println(io, "minAugerEnergy:                $(settings.minAugerEnergy)  ")
        println(io, "maxAugerEnergy:                $(settings.maxAugerEnergy)  ")
        println(io, "maxKappa:                      $(settings.maxKappa)  ")
        println(io, "operator:                      $(settings.operator)  ")
    end


    """
    `struct  Channel`   
        ... defines a type for a AutoIonization channel to help characterize a scattering (continuum) state of many 
            electron-states with a single free electron.

        + kappa          ::Int64                ... partial-wave of the free electron
        + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
        + phase          ::Float64              ... phase of the partial wave
        + amplitude      ::Complex{Float64}     ... Auger amplitude associated with the given channel.
    """
    struct  Channel
        kappa            ::Int64
        symmetry         ::LevelSymmetry
        phase            ::Float64
        amplitude        ::Complex{Float64}
    end


    """
    `struct  Line`  
        ... defines a type for a AutoIonization line that may include the definition of sublines and their 
            corresponding amplitudes.

        + initialLevel   ::Level           ... initial-(state) level
        + finalLevel     ::Level           ... final-(state) level
        + electronEnergy ::Float64         ... Energy of the (incoming free) electron.
        + totalRate      ::Float64         ... Total rate of this line.
        + angularAlpha   ::Float64         ... Angular alpha_2 coefficient.
        + hasChannels    ::Bool            
            ... Determines whether the individual scattering (sub-) channels are defined in terms of their free-electron energy, kappa 
                and the total angular momentum/parity as well as the amplitude, or not.
        + channels       ::Array{AutoIonization.Channel,1}  ... List of AutoIonization channels of this line.
    """
    struct  Line
        initialLevel     ::Level
        finalLevel       ::Level
        electronEnergy   ::Float64
        totalRate        ::Float64
        angularAlpha     ::Float64
        hasChannels      ::Bool
        channels         ::Array{AutoIonization.Channel,1}
    end 


    """
    `AutoIonization.Line(initialLevel::Level, finalLevel::Level, totalRate::Float64)`  
        ... constructor for an AutoIonization line between a specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, totalRate::Float64)
        Line(initialLevel, finalLevel, 0., totalRate, 0., false, AutoIonization.Channel[])
    end


    # `Base.show(io::IO, line::AutoIonization.Line)`  ... prepares a proper printout of the variable line::AutoIonization.Line.
    function Base.show(io::IO, line::AutoIonization.Line) 
        println(io, "initialLevel:           $(line.initialLevel)  ")
        println(io, "finalLevel:             $(line.finalLevel)  ")
        println(io, "electronEnergy:         $(line.electronEnergy)  ")
        println(io, "totalRate:              $(line.totalRate)  ")
        println(io, "angularAlpha:           $(line.angularAlpha)  ")
        println(io, "hasChannels:            $(line.hasChannels)  ")
        println(io, "channels:               $(line.channels)  ")
    end


    """
    `AutoIonization.amplitude(kind::String, channel::AutoIonization.Channel, continuumLevel::Level, initialLevel::Level, 
                              grid::Radial.Grid; printout::Bool=true)`  
        ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit()   Auger amplitude 
            <(alpha_f J_f, kappa) J_i || O^(Auger, kind) || alpha_i J_i>  due to the interelectronic interaction for the given 
            final and initial level. A value::ComplexF64 is returned.
    """
    function amplitude(kind::AbstractEeInteraction, channel::AutoIonization.Channel, continuumLevel::Level, initialLevel::Level, grid::Radial.Grid; 
                       printout::Bool=true)
        nt = length(continuumLevel.basis.csfs);    ni = length(initialLevel.basis.csfs);    partial = Subshell(9,channel.kappa)
        if  printout  printstyled("Compute ($kind) Auger matrix of dimension $nt x $ni in the continuum- and initial-state bases " *
                                  "for the transition [$(initialLevel.index)- ...] and for partial wave $(string(partial)[2:end]) ... ", 
                                  color=:light_green)    end
        matrix = zeros(ComplexF64, nt, ni)
        #
        if      kind in [ CoulombInteraction(), BreitInteraction(), CoulombBreit()]        ## pure V^Coulomb interaction
        #--------------------------------------------------------------------------
            for  r = 1:nt
                for  s = 1:ni
                    if  initialLevel.basis.csfs[s].J != initialLevel.J  ||  initialLevel.basis.csfs[s].parity != initialLevel.parity      continue    end 
                    ##x wa = compute("angular coefficients: e-e, Ratip2013", continuumLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                    # Calculate the spin-angular coefficients
                    if  Defaults.saRatip()
                        waR = compute("angular coefficients: e-e, Ratip2013", continuumLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                        wa  = waR       
                    end
                    if  Defaults.saGG()
                        subshellList = continuumLevel.basis.subshells
                        opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                        waG2 = SpinAngular.computeCoefficients(opa, continuumLevel.basis.csfs[r], initialLevel.basis.csfs[s], subshellList)
                        wa   = [1.0, waG2]
                    end
                    if  Defaults.saRatip() && Defaults.saGG() && true
                        if  length(waR[2]) != 0     println("\n>> Angular coeffients from Ratip2013   = $(waR[2]) ")    end
                        if  length(waG2)   != 0     println(  ">> Angular coeffients from SpinAngular = $waG2 ")        end
                    end
                    #
                    me = 0.
                    for  coeff in wa[2]
                        if   kind in [ CoulombInteraction(), CoulombBreit()]    
                            me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, 
                                                    continuumLevel.basis.orbitals[coeff.a], continuumLevel.basis.orbitals[coeff.b],
                                                    initialLevel.basis.orbitals[coeff.c],   initialLevel.basis.orbitals[coeff.d], grid)   end
                        if   kind in [ BreitInteraction(), CoulombBreit()]    
                            me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, 
                                                    continuumLevel.basis.orbitals[coeff.a], continuumLevel.basis.orbitals[coeff.b],
                                                    initialLevel.basis.orbitals[coeff.c],   initialLevel.basis.orbitals[coeff.d], grid)   end
                    end
                    matrix[r,s] = me
                end
            end 
            if  printout  printstyled("done. \n", color=:light_green)    end
            amplitude = transpose(continuumLevel.mc) * matrix * initialLevel.mc 
            amplitude = im^Basics.subshell_l(Subshell(101, channel.kappa)) * exp( -im*channel.phase ) * amplitude
            #
            #
         elseif  kind == "H-E"
        #--------------------
            iLevel = finalLevel;   fLevel = initialLevel
            amplitude = 0.;    error("stop a")
        else    error("stop b")
        end
        
        return( amplitude )
    end


    """
    `AutoIonization.channelAmplitude(kind::String, channel::AutoIonization.Channel, energy::Float64, finalLevel::Level, 
                                     initialLevel::Level, grid::Radial.Grid)`  
        ... to compute the kind = (CoulombInteraction(), BreitInteraction(), CoulombBreit())   Auger amplitude  
            <(alpha_f J_f, kappa) J_i || O^(Auger, kind) || alpha_i J_i>  due to the interelectronic interaction for the given final and 
            initial level. A newChannel::AutoIonization.Channel is returned.
    """
    function channelAmplitude(kind::String, channel::AutoIonization.Channel, energy::Float64, finalLevel::Level, 
                              initialLevel::Level, grid::Radial.Grid)
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(initialLevel, initialLevel.basis.subshells)
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(finalLevel, finalLevel.basis.subshells)
        cOrbital, phase  = Continuum.generateOrbital(energy, Subshell(101, channel.kappa), newfLevel, grid, contSettings)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
        newChannel = AutoIonization.Channel(channel.kappa, channel.symmetry, phase, 0.)
        amplitude = AutoIonization.amplitude(kind, channel, newcLevel, newiLevel, grid)

        newChannel = AutoIonization.Channel(newChannel.kappa, newChannel.symmetry, newChannel.phase, amplitude)    
        return( newChannel )
    end


    """
    `AutoIonization.computeAmplitudesProperties(line::AutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                settings::AutoIonization.Settings; printout::Bool=true)` 
        ... to compute all amplitudes and properties of the given line; a line::AutoIonization.Line is returned for which the amplitudes 
            and properties are now evaluated.
    """
    function computeAmplitudesProperties(line::AutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, settings::AutoIonization.Settings; 
                                         printout::Bool=true) 
        newChannels = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum);   rate = 0.
        for channel in line.channels
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, newiLevel.basis.subshells)
            newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
            cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
            newChannel = AutoIonization.Channel(channel.kappa, channel.symmetry, phase, 0.)
            amplitude  = AutoIonization.amplitude(settings.operator, newChannel, newcLevel, newiLevel, grid, printout=printout)
            rate       = rate + conj(amplitude) * amplitude
            push!( newChannels, AutoIonization.Channel(newChannel.kappa, newChannel.symmetry, newChannel.phase, amplitude) )
        end
        totalRate = 2pi* rate;   angularAlpha = 0.
        newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, angularAlpha, true, newChannels)
        #
        if  settings.calcAnisotropy    angularAlpha = AutoIonization.computeIntrinsicAlpha(2, newLine)
            newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, real(angularAlpha), true, newChannels)
        end
        
        return( newLine )
    end


    """
    `AutoIonization.computeAmplitudesPropertiesPlasma(line::AutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                      settings::PlasmaShift.AugerSettings; printout::Bool=true)`  
        ... to compute all amplitudes and properties of the given line but for the given plasma model; a line::AutoIonization.Line is returned 
            for which the amplitudes and properties are now evaluated.
    """
    function computeAmplitudesPropertiesPlasma(line::AutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                               settings::PlasmaShift.AugerSettings; printout::Bool=true) 
        newChannels = AutoIonization.Channel[];   contSettings = Continuum.Settings(false, nrContinuum);   rate = 0.
        for channel in line.channels
            newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
            newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
            newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, line.finalLevel.basis.subshells)
            @warn "Adapt a proper continuum orbital for the plasma potential"
            cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
            newChannel = AutoIonization.Channel(channel.kappa, channel.symmetry, phase, 0.)
            @warn "Adapt a proper Auger amplitude for the plasma e-e interaction"
            amplitude = 1.0
            # amplitude  = AutoIonization.amplitude(settings.operator, newChannel, newcLevel, newiLevel, grid)
            rate       = rate + conj(amplitude) * amplitude
            push!( newChannels, AutoIonization.Channel(newChannel.kappa, newChannel.symmetry, newChannel.phase, amplitude) )
        end
        totalRate = 2pi* rate;   angularAlpha = 0.
        newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, angularAlpha, true, newChannels)
        
        return( newLine )
    end



    """
    `AutoIonization.computeIntrinsicAlpha(k::Int64, line::AutoIonization.Line)`  
        ... to compute the intrinsic alpha_k anisotropy parameter for the given line. A value::Float64 is returned.
    """
    function  computeIntrinsicAlpha(k::Int64, line::AutoIonization.Line)
        if  !line.hasChannels   error("No channels are defined for the given AutoIonization.line.")                   end
        wn = 0.;    for  channel in line.channels    wn = wn + conj(channel.amplitude) * channel.amplitude   end
        wa = 0.;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J;
        for  cha  in line.channels  
            j = AngularMomentum.kappa_j(cha.kappa);    l = AngularMomentum.kappa_l(cha.kappa)
            for  chp  in line.channels  
                jp = AngularMomentum.kappa_j(chp.kappa);    lp = AngularMomentum.kappa_l(chp.kappa)
                wa = wa + sqrt( AngularMomentum.bracket([l, lp, j, jp]) ) *  
                          AngularMomentum.ClebschGordan(l, AngularM64(0), lp, AngularM64(0), AngularJ64(k), AngularM64(0)) *
                          AngularMomentum.Wigner_6j(Ji, j, Jf, jp, Ji, AngularJ64(k)) * 
                          AngularMomentum.Wigner_6j(l,  j, AngularJ64(1//2), jp, lp, AngularJ64(k)) * 
                          cha.amplitude * conj(chp.amplitude)
            end    
        end
        value = AngularMomentum.phaseFactor([Ji, +1, Jf, +1, AngularJ64(k), -1, AngularJ64(1//2)]) * 
                sqrt(Basics.twice(Ji) + 1) * wa / wn

        if  false
            # Calculate the value given by M.H. Chen, PRA 47 (1993) 3733; this does not agree exactly so far.
            wa = 0.
            for  cha  in line.channels  
                j = AngularMomentum.kappa_j(cha.kappa);    l = AngularMomentum.kappa_l(cha.kappa);    phase = cha.phase
                for  chp  in line.channels  
                    jp = AngularMomentum.kappa_j(chp.kappa);    lp = AngularMomentum.kappa_l(chp.kappa);    phasep = chp.phase
                    wa = wa + AngularMomentum.phaseFactor([Ji, +1, Jf, -1, AngularJ64(1//2)]) * im^(Float64(l) - Float64(lp)) * cos(phase - phasep) *
                              sqrt( AngularMomentum.bracket([l, lp, j, jp, AngularJ64(k), Ji]) ) *  
                              AngularMomentum.Wigner_3j(lp, l, AngularJ64(k), AngularJ64(0), AngularJ64(0), AngularJ64(0)) * 
                              AngularMomentum.Wigner_6j(j, jp, AngularJ64(k), lp, l, AngularJ64(1//2)) * 
                              AngularMomentum.Wigner_6j(Ji, Ji, AngularJ64(k), jp, j, Jf) * 
                              cha.amplitude * conj(chp.amplitude)
                end
            end
            println("*** Comparison value = $value    Chen-value = $(wa/wn)")
        end

        return( value )
    end



    """
    `AutoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                 settings::AutoIonization.Settings; output=true, printout::Bool=true)`  
        ... to compute the Auger transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{AutoIonization.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::AutoIonization.Settings; output=true, printout::Bool=true)
        println("")
        printstyled("AutoIonization.computeLines(): The computation of Auger rates and properties starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        lines = AutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    AutoIonization.displayLines(lines)    end  
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = AutoIonization.Line[]
        for  line in lines
            newLine = AutoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        AutoIonization.displayRates(stdout, newLines, settings)
        AutoIonization.displayLifetimes(stdout, newLines)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   AutoIonization.displayRates(iostream, newLines, settings);   AutoIonization.displayLifetimes(iostream, newLines)     end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end



    """
    `AutoIonization.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                        settings::AutoIonization.Settings; output::Bool=true, printout::Bool=true)`  
        ... to compute the Auger transition amplitudes and all properties as requested by the given settings. The computations
            and printout is adapted for large cascade computations by including only lines with at least one channel and by sending
            all printout to a summary file only. A list of lines::Array{AutoIonization.Lines} is returned.
    """
    function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                  settings::AutoIonization.Settings; output::Bool=true, printout::Bool=true)
        
        lines = AutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        # if  settings.printBefore    AutoIonization.displayLines(lines)    end  
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = AutoIonization.Line[]
        for  (i,line)  in  enumerate(lines)
            if  rem(i,10) == 0    println("> Auger line $i:  ... calculated ")    end
            newLine = AutoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings, printout=printout) 
            push!( newLines, newLine)
        end
        # Print all results to a summary file, if requested
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   AutoIonization.displayRates(iostream, newLines, settings)     end
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end



    """
    `AutoIonization.computeLinesFromOrbitals(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                             settings::AutoIonization.Settings, contOrbitals::Dict{Subshell, Orbital}; output::Bool=true, printout::Bool=true)`  
        ... to compute the Auger transition amplitudes and all properties as requested by the given settings but by using the given set of 
            continuum orbitals. The computations and printout is adapted for large cascade computations by including only lines with at least 
            one channel and by sending all printout to a summary file only. A list of lines::Array{AutoIonization.Lines} is returned.
    """
    function  computeLinesFromOrbitals(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                       settings::AutoIonization.Settings, contOrbitals::Dict{Subshell, Orbital}; output::Bool=true, printout::Bool=true)

        lines = AutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
        # Calculate all amplitudes and requested properties
        newLines = AutoIonization.Line[]
        for  (i,line)  in  enumerate(lines)
            if  rem(i,500) == 0    println("> Auger line $i:  ... calculated ")    end
            # Calculate the individual channels with the given orbitals
            newChannels = AutoIonization.Channel[];   rate = 0.
            for channel in line.channels
                newiLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
                newfLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, newiLevel.basis.subshells)
                sh         = Subshell(101, channel.kappa)
                if haskey(contOrbitals, sh)   cOrbital = contOrbitals[sh]      else    println(">>> skip Auger channel for $sh");   continue    end
                newiLevel  = Basics.generateLevelWithExtraSubshell(sh, newiLevel)
                newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
                newChannel = AutoIonization.Channel(channel.kappa, channel.symmetry, 0., 0.)
                amplitude  = AutoIonization.amplitude(settings.operator, newChannel, newcLevel, newiLevel, grid, printout=printout)
                rate       = rate + conj(amplitude) * amplitude
                push!( newChannels, AutoIonization.Channel(newChannel.kappa, newChannel.symmetry, newChannel.phase, amplitude) )
            end
            totalRate = 2pi* rate;   angularAlpha = 0.
            newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, angularAlpha, true, newChannels)
            ##
            ## if  settings.calcAnisotropy    angularAlpha = AutoIonization.computeIntrinsicAlpha(2, newLine)
            ##  newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, real(angularAlpha), true, newChannels)
            ## end
            push!( newLines, newLine)
        end
        # Print all results to a summary file, if requested
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   AutoIonization.displayRates(iostream, newLines, settings)     end
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end



    """
    `AutoIonization.computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                       settings::PlasmaShift.AugerSettings; output=true)`  
        ... to compute the Auger transition amplitudes and all properties as requested by the given settings. A list of 
            lines::Array{AutoIonization.Lines} is returned.
    """
    function  computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                 settings::PlasmaShift.AugerSettings; output=true)
        println("")
        printstyled("AutoIonization.computeLinesPlasma(): The computation of Auger rates starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        augerSettings = AutoIonization.Settings(false, settings.printBefore, settings.lineSelection, 0., 1.0e6, 100, CoulombInteraction())
        lines         = AutoIonization.determineLines(finalMultiplet, initialMultiplet, augerSettings)
        # Display all selected lines before the computations start
        if  settings.printBefore    AutoIonization.displayLines(lines)    end
        # Determine maximum energy and check for consistency of the grid
        maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
        nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
        # Calculate all amplitudes and requested properties
        newLines = AutoIonization.Line[]
        for  line in lines
            newLine = AutoIonization.computeAmplitudesPropertiesPlasma(line, nm, grid, nrContinuum, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        AutoIonization.displayRates(stdout, newLines, augerSettings)
        AutoIonization.displayLifetimes(stdout, newLines)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   AutoIonization.displayRates(iostream, newLines, augerSettings);   AutoIonization.displayLifetimes(iostream, newLines)     end
        #
        if    output    return( lines )
        else            return( nothing )
        end
    end


    """
    `AutoIonization.determineChannels(finalLevel::Level, initialLevel::Level, settings::AutoIonization.Settings)`  
        ... to determine a list of Auger Channel for a transitions from the initial to final level and by taking into account the particular 
            settings of for this computation; an Array{AutoIonization.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::AutoIonization.Settings)
        channels  = AutoIonization.Channel[];   
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        kappaList = AngularMomentum.allowedKappaSymmetries(symi, symf)
        for  kappa in kappaList
            if  abs(kappa) > settings.maxKappa      continue    end
            push!(channels, AutoIonization.Channel(kappa, symi, 0., Complex(0.)) )
        end
        return( channels )  
    end


    """
    `AutoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::AutoIonization.Settings)`  
        ... to determine a list of AutoIonization.Line's for transitions between levels from the initial- and final-state multiplets, and  
            by taking into account the particular selections and settings for this computation; an Array{AutoIonization.Line,1} is returned. 
            Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::AutoIonization.Settings)
        lines = AutoIonization.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    energy = iLevel.energy - fLevel.energy
                    if   energy < 0.01                                                             continue   end
                    if   energy < settings.minAugerEnergy  ||  energy > settings.maxAugerEnergy    continue   end  
                    channels = AutoIonization.determineChannels(fLevel, iLevel, settings) 
                    push!( lines, AutoIonization.Line(iLevel, fLevel, energy, 0., 0., true, channels) )
                end
            end
        end
        return( lines )
    end


    """
    `AutoIonization.displayLines(lines::Array{AutoIonization.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
            transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{AutoIonization.Line,1})
        nx = 150
        println(" ")
        println("  Selected Auger lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(14, "Energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(14, "Energy e_A"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.flushleft(37, "List of kappas and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(37, "partial (total J^P)                "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.electronEnergy))       * "   "
            kappaSymmetryList = Tuple{Int64,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaSymmetryList, (line.channels[i].kappa, line.channels[i].symmetry) )
            end
            sa = sa * TableStrings.kappaSymmetryTupels(80, kappaSymmetryList)
            println( sa )
        end
        println("  ", TableStrings.hLine(nx), "\n")
        #
        return( nothing )
    end


    """
    `AutoIonization.displayLifetimes(stream::IO, lines::Array{AutoIonization.Line,1})`  
        ... to list all lifetimes as associated with the selected lines. A neat table is printed but nothing is returned otherwise.
    """
    function  displayLifetimes(stream::IO, lines::Array{AutoIonization.Line,1})
        nx = 104
        println(stream, " ")
        println(stream, "  Auger lifetimes, total rates and widths:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level";    na=2);                           sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center( 8, "J^P";      na=4);                           sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(12, "Lifetime"; na=4);               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("time"); na=4)
        sa = sa * TableStrings.center(14, "Total rate"; na=6);               
        sb = sb * TableStrings.center(14,TableStrings.inUnits("rate"); na=4)
        sa = sa * TableStrings.center(42, "Widths"; na=2);       
        sb = sb * TableStrings.center(42, "  Hartrees         Kaysers           eV"; na=2)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        # 
        notYetDone = trues(1000)
        for  line in lines
            if  notYetDone[line.initialLevel.index]
                notYetDone[line.initialLevel.index] = false
                totalRate = 0.
                for  ln in lines
                    if  ln.initialLevel.index == line.initialLevel.index   totalRate = totalRate + ln.totalRate    end
                end
                sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                sa = sa * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
                sa = sa * TableStrings.center( 8, string(isym); na=4)
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",  1/totalRate))            * "     "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic",    totalRate))            * "      "
                sa = sa * @sprintf("%.6e", totalRate)                                                           * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to Kayser",  totalRate))  * "    "
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to eV",      totalRate))  * "    "
                println(stream, sa)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `AutoIonization.displayRates(stream::IO, lines::Array{AutoIonization.Line,1}, settings::AutoIonization.Settings)`  
        ... to list all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is returned 
            otherwise.
    """
    function  displayRates(stream::IO, lines::Array{AutoIonization.Line,1}, settings::AutoIonization.Settings)
        nx = 106
        println(stream, " ")
        if  settings.calcAnisotropy    println(stream, "  Auger rates and intrinsic angular parameters: \n")
        else                           println(stream, "  Auger rates (without angular parameters): \n")        end
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "Energy"   ; na=2);               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(14, "Electron energy"   ; na=2);               
        sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
        sa = sa * TableStrings.center(14, "Auger rate"; na=2);       
        sb = sb * TableStrings.center(14, TableStrings.inUnits("rate"); na=2)
        sa = sa * TableStrings.center(15, "alpha_2"; na=2);                           sb = sb * TableStrings.hBlank(18)     
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy))       * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.totalRate))              * "    "
            sa = sa * TableStrings.flushright(13, @sprintf("%.4e", line.angularAlpha))            * "    "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

end # module
