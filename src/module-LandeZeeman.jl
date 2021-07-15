
"""
`module  JAC.LandeZeeman`  
    ... a submodel of JAC that contains all methods for computing Lande factors and Zeeman properties for some level(s).
"""
module LandeZeeman

    using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial,
                  ..SpinAngular, ..TableStrings


    """
    `struct  LandeZeeman.SublevelJ`  ... defines a type to specify a magnetic sublevel with well-defined J.

        + M                      ::AngularM64        ... M_J-value
        + energy                 ::Float64           ... energy of this sublevel
    """
    struct SublevelJ 
        M                        ::AngularM64
        energy                   ::Float64
    end 


    # `Base.show(io::IO, Jsublevel::LandeZeeman.SublevelJ)`  ... prepares a proper printout of the variable Jsublevel::LandeZeeman.SublevelJ.
    function Base.show(io::IO, Jsublevel::LandeZeeman.SublevelJ) 
        println(io, "Sublevel [M=$(Jsublevel.M); energy = $(Jsublevel.energy)]")
    end



    # SublevelF defines a magnetic hyperfine sublevel
    """
    `struct  LandeZeeman.SublevelF`  ... defines a type to specify a magnetic hyperfine sublevel with well-defined F.

        + F                      ::AngularJ64        ... F-value
        + M                      ::AngularM64        ... M_F-value
        + energy                 ::Float64           ... energy of this sublevel
    """
    struct SublevelF
        F                        ::AngularJ64
        M                        ::AngularM64
        energy                   ::Float64
    end 


    # `Base.show(io::IO, Fsublevel::LandeZeeman.SublevelF)`  ... prepares a proper printout of the variable Fsublevel::LandeZeeman.SublevelF.
    function Base.show(io::IO, Fsublevel::LandeZeeman.SublevelF) 
        println(io, "Sublevel [F=$(Fsublevel.F); M=$(Fsublevel.M); energy = $(Fsublevel.energy)]")
    end


    """
    `struct  LandeZeeman.Outcome`  
        ... defines a type to keep the Lande Factors and Zeeman spittling parameters of a fine-structure level.

        + Jlevel                 ::Level             ... Fine-structure levels to which the results refer to.
        + LandeJ                 ::Float64           ... Lande-J factor of this fine-structure level 
        + amplitudeN1            ::Complex{Float64}  ... N1 amplitude
        + amplitudeDeltaN1       ::Complex{Float64}  ... Delta-N1 amplitude
        + nuclearI               ::AngularJ64        ... nuclear spin
        + hasJsublevels          ::Bool              ... True, if information about the magnetic sublevels are provided and false otherwise.
        + Jsublevels             ::Array{LandeZeeman.SublevelJ,1}   ... List of the magnetic fine-structure sublevels and data with 
                                                                        well-defined J-value 
        + hasFsublevels          ::Bool              ... True, if information about the magnetic hyperfine sublevels are provided and 
                                                         false otherwise.
        + Fsublevels             ::Array{LandeZeeman.SublevelF,1}   ... List of the magnetic hyperfine sublevels and data with 
                                                                        well-defined F-value 
    """
    struct Outcome 
        Jlevel                   ::Level
        LandeJ                   ::Float64
        amplitudeN1              ::Complex{Float64}
        amplitudeDeltaN1         ::Complex{Float64}
        nuclearI                 ::AngularJ64
        hasJsublevels            ::Bool
        Jsublevels               ::Array{LandeZeeman.SublevelJ,1}
        hasFsublevels            ::Bool
        Fsublevels               ::Array{LandeZeeman.SublevelF,1}
    end 


    """
    `LandeZeeman.Outcome()`  ... constructor for an `empty` instance of LandeZeeman.Outcome.
    """
    function Outcome()
        Outcome(Level(), 0., 0., 0., AngularJ64(0), false, SublevelJ[],  false, SublevelF[])
    end


    # `Base.show(io::IO, outcome::LandeZeeman.Outcome)`  ... prepares a proper printout of the variable LandeZeeman.Outcome.
    function Base.show(io::IO, outcome::LandeZeeman.Outcome) 
        println(io, "Jlevel:           $(outcome.Jlevel)  ")
        println(io, "LandeJ:           $(outcome.LandeJ)  ")
        println(io, "amplitudeN1:      $(outcome.amplitudeN1)  ")
        println(io, "amplitudeDeltaN1: $(outcome.amplitudeDeltaN1)  ")
        println(io, "nuclearI:         $(outcome.nuclearI)  ")
        println(io, "hasJsublevels:    $(outcome.hasJsublevels)  ")
        println(io, "Fsublevels:       $(outcome.Fsublevels)  ")
        println(io, "hasFsublevels:    $(outcome.hasFsublevels)  ")
        println(io, "Fsublevels:       $(outcome.Fsublevels)  ")
    end


    """
    `struct  LandeZeeman.Settings  <:  AbstractPropertySettings`  
        ... defines a type for the details and parameters of computing the Lande Factors and Zeeman spittling 
            of fine-structure levels.

        + calcLandeJ             ::Bool              ... True if Lande-J factors are to be calculated.
        + calcLandeF             ::Bool              ... True if Lande-F factors are to be calculated.
        + calcZeeman             ::Bool              ... True if Zeeman splittings are to be calculated.
        + includeSchwinger       ::Bool              ... True if Schwinger's QED correction ``\\Delta N^(1)`` is to be included, 
                                                         and false otherwise.
        + BField                 ::Float64           ... Strength of the magnetic field in [Tesla]
        + printBefore            ::Bool              ... True if a list of selected levels is printed before the actual computations start. 
        + levelSelection         ::LevelSelection    ... Specifies the selected levels, if any.
    """
    struct Settings  <:  AbstractPropertySettings 
        calcLandeJ               ::Bool
        calcLandeF               ::Bool
        calcZeeman               ::Bool
        includeSchwinger         ::Bool
        BField                   ::Float64
        printBefore              ::Bool 
        levelSelection           ::LevelSelection
     end 


    """
    `LandeZeeman.Settings()`  
        ... constructor for an `empty` instance of ZeemanSettings for the computation of isotope M and F parameters.
     """
     function Settings()
         Settings(false, false, false, false, 0., false, LevelSelection() )
    end


    # `Base.show(io::IO, settings::LandeZeeman.Settings)`  ... prepares a proper printout of the variable settings::LandeZeeman.Settings.
    function Base.show(io::IO, settings::LandeZeeman.Settings) 
        println(io, "calcLandeJ:               $(settings.calcLandeJ)  ")
        println(io, "calcLandeF:               $(settings.calcLandeF)  ")
        println(io, "calcZeeman:               $(settings.calcZeeman)  ")
        println(io, "includeSchwinger:         $(settings.includeSchwinger)  ")
        println(io, "BField:                   $(settings.BField)  ")
        println(io, "printBefore:              $(settings.printBefore)  ")
        println(io, "levelSelection:           $(settings.levelSelection)  ")
    end


    """
    `LandeZeeman.amplitude(kind::String, rLevel::Level, sLevel::Level, grid::Radial.Grid)` 
        ... to compute either the  N^(1) or Delta N^(1) Zeeman amplitude <alpha_r J_r || N^(1)) || alpha_s J_s>  
            for a given pair of levels. A value::ComplexF64 is returned.
    """
    function  amplitude(kind::String, rLevel::Level, sLevel::Level, grid::Radial.Grid)
        nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
        printstyled("Compute Zeeman $(kind[1:5]) matrix of dimension $nr x $ns ... ", color=:light_green)
        #
        if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
        #
        for  r = 1:nr
            for  s = 1:ns
                me = 0.
                if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                    rLevel.parity != sLevel.parity    continue    
                end 
                #
                if      kind == "N^(1) amplitude"
                #--------------------------------
                    ##x wa = compute("angular coefficients: 1-p, Grasp92", 0, 1, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    # Calculate the spin-angular coefficients
                    if  Defaults.saRatip()
                        waR = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 1, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                        wa  = waR       
                    end
                    if  Defaults.saGG()
                        subshellList = sLevel.basis.subshells
                        opa = SpinAngular.OneParticleOperator(1, plus, true)
                        waG = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList) 
                        wa  = waG
                    end
                    if  Defaults.saRatip() && Defaults.saGG() && true
                        if  length(waR) != 0     println("\n>> Angular coeffients from GRASP/MCT   = $waR ")    end
                        if  length(waG) != 0     println(  ">> Angular coeffients from SpinAngular = $waG ")    end
                    end
                    #
                    for  coeff in wa
                        ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = InteractionStrength.zeeman_n1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                        me = me + coeff.T * tamp  
                    end
                #
                elseif  kind == "Delta N^(1) amplitude"
                #--------------------------------------
                    ##x wa = compute("angular coefficients: 1-p, Grasp92", 0, 1, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    # Calculate the spin-angular coefficients
                    if  Defaults.saRatip()
                        waR = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 1, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                        wa  = waR       
                    end
                    if  Defaults.saGG()
                        subshellList = sLevel.basis.subshells
                        opa = SpinAngular.OneParticleOperator(1, plus, true)
                        waG = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList) 
                        wa  = waG
                    end
                    if  Defaults.saRatip() && Defaults.saGG() && true
                        if  length(waR) != 0     println("\n>> Angular coeffients from GRASP/MCT   = $waR ")    end
                        if  length(waG) != 0     println(  ">> Angular coeffients from SpinAngular = $waG ")    end
                    end
                    #
                    for  coeff in wa
                        ja = Basics.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = Basics.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = InteractionStrength.zeeman_Delta_n1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                        me = me + coeff.T * tamp  
                    end
                #
                else    error("stop a")
                end
                #
                matrix[r,s] = me
            end
        end
        printstyled("done.\n", color=:light_green)
        amplitude = transpose(rLevel.mc) * matrix * sLevel.mc 
        #
        return( amplitude )
    end



    """
    `LandeZeeman.computeAmplitudesProperties(outcome::LandeZeeman.Outcome, grid::Radial.Grid, settings::LandeZeeman.Settings)`  
        ... to compute all amplitudes and properties of for a given level; an outcome::LandeZeeman.Outcome is returned for 
            which the amplitudes and all requested properties are now evaluated explicitly.
    """
    function  computeAmplitudesProperties(outcome::LandeZeeman.Outcome, grid::Radial.Grid, settings::LandeZeeman.Settings)
        LandeJ = 0.0;   amplitudeN1 = 0.0;   amplitudeDeltaN1 = 0.0;    J = AngularMomentum.oneJ(outcome.Jlevel.J) 
        amplitudeN1 = LandeZeeman.amplitude("N^(1) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
        #
        if  settings.includeSchwinger
            amplitudeDeltaN1 = LandeZeeman.amplitude("Delta N^(1) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
        end
        #
        if  settings.calcLandeJ    LandeJ = 2*(amplitudeN1 + amplitudeDeltaN1) / sqrt(J*(J+1))    end

        
        newOutcome = LandeZeeman.Outcome( outcome.Jlevel, LandeJ, amplitudeN1, amplitudeDeltaN1, outcome.nuclearI, 
                                          outcome.hasJsublevels, outcome.Jsublevels, outcome.hasFsublevels, outcome.Fsublevels)
        return( newOutcome )
    end


    """
    `LandeZeeman.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::LandeZeeman.Settings; output=true)`  
        ... to compute (as selected) the Lande J or F factors for the levels of the given multiplet and as specified by the given settings. The 
            results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::LandeZeeman.Settings; output=true)
        println("")
        printstyled("LandeZeeman.computeOutcomes(): The computation of the Zeeman amplitudes and Lande factors starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        outcomes = LandeZeeman.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBefore    LandeZeeman.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = LandeZeeman.Outcome[]
        for  outcome in outcomes
            newOutcome = LandeZeeman.computeAmplitudesProperties(outcome, grid, settings) 
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        LandeZeeman.displayResults(stdout, newOutcomes, nm, settings)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary    LandeZeeman.displayResults(iostream, newOutcomes, nm, settings)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `LandeZeeman.determineOutcomes(multiplet::Multiplet, settings::LandeZeeman.Settings)`  
        ... to determine a list of Outcomes's for the computation of Lande factors and/or Zeeman splittings for 
            levels from the given multiplet. It takes into account the particular selections and settings. 
            An Array{LandeZeeman.Outcome,1} is returned. Apart from the level specification, all physical properties are 
            still set to zero during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::LandeZeeman.Settings) 
        # Define values that depend on the requested computations
        nuclearI = AngularJ64(0);    hasJsublevels = false;    Jsublevels = LandeZeeman.SublevelJ[]
                                     hasFsublevels = false;    Fsublevels = LandeZeeman.SublevelF[]

        outcomes = LandeZeeman.Outcome[]
        for  level  in  multiplet.levels
            if  Basics.selectLevel(level, settings.levelSelection)
                push!( outcomes, LandeZeeman.Outcome(level, 0., 0., 0., nuclearI, hasJsublevels, Jsublevels, hasFsublevels, Fsublevels) )
            end
        end
        return( outcomes )
    end


    """
    `LandeZeeman.displayOutcomes(outcomes::Array{LandeZeeman.Outcome,1})`  
        ... to display a list of levels that have been selected for the computations. A small neat table of all selected 
            levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{LandeZeeman.Outcome,1})
        nx = 43
        println(" ")
        println("  Selected Lande-Zeeman levels:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
        sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
        sa = sa * TableStrings.center(14, "Energy"; na=4);              
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.Jlevel.energy)) * "    "
            println( sa )
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `LandeZeeman.displayResults(stream::IO, outcomes::Array{LandeZeeman.Outcome,1}, nm::Nuclear.Model, settings::LandeZeeman.Settings)`  
        ... to display the energies, Lande factors, Zeeman amplitudes etc. for the selected levels. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{LandeZeeman.Outcome,1}, nm::Nuclear.Model, settings::LandeZeeman.Settings)
        #
        if  settings.calcLandeJ
            nx = 135
            println(stream, " ")
            println(stream, "  Lande g_J factors and Zeeman amplitudes:")
            println(stream, " ")
            println(stream, "  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
            sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
            sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
            sa = sa * TableStrings.center(14, "Energy"; na=5)              
            sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
            sa = sa * TableStrings.center(12, "Lande-J"; na=8);                           sb = sb * TableStrings.hBlank(20)          
            sa = sa * TableStrings.center(62, "N1   -- Zeeman Amplitudes --   Delta N1"    ; na=4) 
            sb = sb * TableStrings.center(62, "   Re              Im                 Re              Im"; na=5)
            println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
            #  
            for  outcome in outcomes
                sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
                sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                energy = 1.0
                sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))                    * "    "
                sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.LandeJ) )              * "    "
                sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeN1.re) )      * " "
                sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeN1.im) )      * "    "
                sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeDeltaN1.re) ) * " "
                sa = sa * TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeDeltaN1.im) ) * "    "
                println(stream, sa )
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        if  settings.calcLandeF
            nx = 80
            println(stream, " ")
            println(stream, "  Hyperfine levels and Lande g_F factors:")
            println(stream, " ")
            println(stream, "  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
            sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
            sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
            sa = sa * TableStrings.center(14, "Energy"; na=5)              
            sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=5)
            sa = sa * TableStrings.center(10, "F^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
            sa = sa * TableStrings.center(12, "Lande-F"; na=8);                           sb = sb * TableStrings.hBlank(20)          
            println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
            #
            for  outcome in outcomes
                sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
                sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                energy = 1.0
                sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))                    * "    "
                println(stream, sa )
                #
                sc = TableStrings.hBlank( length(sa) + 1 )
                Flist = oplus(nm.spinI, outcome.Jlevel.J)
                for  F in Flist
                    symf   = LevelSymmetry( F, outcome.Jlevel.parity);      Fx = AngularMomentum.oneJ(F)
                    Jx     = AngularMomentum.oneJ(outcome.Jlevel.J);    Ix = AngularMomentum.oneJ(nm.spinI)
                    LandeF = (Fx*(Fx+1) + Jx*(Jx+1) - Ix*(Ix+1)) / (2*Fx*(Fx+1)) * outcome.LandeJ
                    sa = sc * TableStrings.center(10, string(symf); na=4)
                    sa = sa * TableStrings.flushright(15, @sprintf("%.8e", LandeF) )   
                    println(stream, sa )
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        if  settings.calcZeeman
            nx = 135
            println(stream, " ")
            println(stream, "  Zeeman splittings of (hyper-) fine-structure levels:")
            println(stream, " ")
            println(stream, "  ", TableStrings.hLine(nx))
            println(stream, "  ", TableStrings.hLine(nx))
        end
        #
        return( nothing )
    end

end # module
