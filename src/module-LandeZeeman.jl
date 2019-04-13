
"""
`module  JAC.LandeZeeman`  ... a submodel of JAC that contains all methods for computing Lande factors and Zeeman properties for some level; 
                               it is using JAC, JAC.ManyElectron, JAC.Radial.
"""
module LandeZeeman

    using Printf, JAC, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0


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
    `struct  LandeZeeman.Outcome`  ... defines a type to keep the Lande Factors and Zeeman spittling parameters of a fine-structure level.

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
    `JAC.LandeZeeman.Outcome()`  ... constructor for an `empty` instance of LandeZeeman.Outcome.
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
    `struct  LandeZeeman.Settings`  ... defines a type for the details and parameters of computing the Lande Factors and Zeeman spittling 
                                        of fine-structure levels.

        + calcLandeJ             ::Bool              ... True if Lande-J factors are to be calculated.
        + calcLandeF             ::Bool              ... True if Lande-F factors are to be calculated.
        + calcZeeman             ::Bool              ... True if Zeeman splittings are to be calculated.
        + includeSchwinger       ::Bool              ... True if Schwinger's QED correction ``\\Delta N^(1)`` is to be included, and false otherwise.
        + BField                 ::Float64           ... Strength of the magnetic field in [Tesla]
        + printBeforeComputation ::Bool              ... True if a list of selected levels is printed before the actual computations start. 
        + selectLevels           ::Bool              ... true, if specific level (number)s have been selected for computation.
        + selectedLevels         ::Array{Int64,1}    ... Level number that have been selected.
    """
    struct Settings 
        calcLandeJ               ::Bool
        calcLandeF               ::Bool
        calcZeeman               ::Bool
        includeSchwinger         ::Bool
        BField                   ::Float64
        printBeforeComputation   ::Bool 
        selectLevels             ::Bool
        selectedLevels           ::Array{Int64,1}
    end 


    """
    `JAC.LandeZeeman.Settings()`  ... constructor for an `empty` instance of ZeemanSettings for the computation of isotope M and F parameters.
     """
     function Settings()
         Settings(false, false, false, false, 0., false, false, Level[])
    end


    # `Base.show(io::IO, settings::LandeZeeman.Settings)`  ... prepares a proper printout of the variable settings::LandeZeeman.Settings.
    function Base.show(io::IO, settings::LandeZeeman.Settings) 
        println(io, "calcLandeJ:               $(settings.calcLandeJ)  ")
        println(io, "calcLandeF:               $(settings.calcLandeF)  ")
        println(io, "calcZeeman:               $(settings.calcZeeman)  ")
        println(io, "includeSchwinger:         $(settings.includeSchwinger)  ")
        println(io, "BField:                   $(settings.BField)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    """
    `JAC.LandeZeeman.amplitude(kind::String, rLevel::Level, sLevel::Level, grid::Radial.Grid)` ... to compute either 
         the  N^(1) or Delta N^(1) Zeeman amplitude <alpha_r J_r || N^(1)) || alpha_s J_s>  for a given pair of levels. 
         A value::ComplexF64 is returned.
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
                    wa = compute("angular coefficients: 1-p, Grasp92", 0, 1, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.zeeman_n1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                        me = me + coeff.T * tamp  
                    end
                #
                elseif  kind == "Delta N^(1) amplitude"
                #--------------------------------------
                    wa = compute("angular coefficients: 1-p, Grasp92", 0, 1, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.zeeman_Delta_n1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
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
    `JAC.LandeZeeman.computeAmplitudesProperties(outcome::LandeZeeman.Outcome, grid::Radial.Grid, settings::LandeZeeman.Settings)`  ... to compute 
         all amplitudes and properties of for a given level; an outcome::LandeZeeman.Outcome is returned for which the amplitudes and all 
         requested properties are now evaluated explicitly.
    """
    function  computeAmplitudesProperties(outcome::LandeZeeman.Outcome, grid::Radial.Grid, settings::LandeZeeman.Settings)
        LandeJ = 0.0;   amplitudeN1 = 0.0;   amplitudeDeltaN1 = 0.0;    J = JAC.AngularMomentum.oneJ(outcome.Jlevel.J) 
        amplitudeN1 = JAC.LandeZeeman.amplitude("N^(1) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
        #
        if  settings.includeSchwinger
            amplitudeDeltaN1 = JAC.LandeZeeman.amplitude("Delta N^(1) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
        end
        #
        if  settings.calcLandeJ    LandeJ = 2*(amplitudeN1 + amplitudeDeltaN1) / sqrt(J*(J+1))    end

        
        newOutcome = LandeZeeman.Outcome( outcome.Jlevel, LandeJ, amplitudeN1, amplitudeDeltaN1, outcome.nuclearI, 
                                          outcome.hasJsublevels, outcome.Jsublevels, outcome.hasFsublevels, outcome.Fsublevels)
        return( newOutcome )
    end


    """
    `JAC.LandeZeeman.computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::LandeZeeman.Settings; output=true)`  
         ... to compute (as selected) the Lande J or F factors for the levels of the given multiplet and as specified by the given settings. The 
             results are printed in neat tables to screen but nothing is returned otherwise.
    """
    function computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::LandeZeeman.Settings; output=true)
        println("")
        printstyled("JAC.LandeZeeman.computeOutcomes(): The computation of the Zeeman amplitudes and Lande factors starts now ... \n", color=:light_green)
        printstyled("------------------------------------------------------------------------------------------------------------ \n", color=:light_green)
        println("")
        outcomes = JAC.LandeZeeman.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    JAC.LandeZeeman.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        newOutcomes = LandeZeeman.Outcome[]
        for  outcome in outcomes
            newOutcome = JAC.LandeZeeman.computeAmplitudesProperties(outcome, grid, settings) 
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        JAC.LandeZeeman.displayResults(stdout, newOutcomes, nm, settings)
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary    JAC.LandeZeeman.displayResults(iostream, newOutcomes, nm, settings)   end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end


    """
    `JAC.LandeZeeman.determineOutcomes(multiplet::Multiplet, settings::LandeZeeman.Settings)`  ... to determine a list of Outcomes's for the 
         computation of Lande factors and/or Zeeman splittings for levels from the given multiplet. It takes into account the particular 
         selections and settings. An Array{LandeZeeman.Outcome,1} is returned. Apart from the level specification, all physical properties are 
         still set to zero during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::LandeZeeman.Settings) 
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end

        # Define values that depend on the requested computations
        nuclearI = AngularJ64(0);    hasJsublevels = false;    Jsublevels = LandeZeeman.SublevelJ[]
                                     hasFsublevels = false;    Fsublevels = LandeZeeman.SublevelF[]

        outcomes = LandeZeeman.Outcome[]
        for  i = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
            push!( outcomes, LandeZeeman.Outcome(multiplet.levels[i], 0., 0., 0., nuclearI, 
                                                 hasJsublevels, Jsublevels, hasFsublevels, Fsublevels) )
        end
        return( outcomes )
    end


    """
    `JAC.LandeZeeman.displayOutcomes(outcomes::Array{LandeZeeman.Outcome,1})`  ... to display a list of levels that have been selected for 
         the computations. A small neat table of all selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{LandeZeeman.Outcome,1})
        println(" ")
        println("  Selected Lande-Zeeman levels:")
        println(" ")
        println("  ", JAC.TableStrings.hLine(43))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        println(sa);    println(sb);    println("  ", JAC.TableStrings.hLine(43)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", outcome.Jlevel.energy)) * "    "
            println( sa )
        end
        println("  ", JAC.TableStrings.hLine(43))
        #
        return( nothing )
    end


    """
    `JAC.LandeZeeman.displayResults(stream::IO, outcomes::Array{LandeZeeman.Outcome,1}, nm::JAC.Nuclear.Model, settings::LandeZeeman.Settings)`  
         ... to display the energies, Lande factors, Zeeman amplitudes etc. for the selected levels. A neat table is printed but nothing is 
             returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{LandeZeeman.Outcome,1}, nm::JAC.Nuclear.Model, settings::LandeZeeman.Settings)
        #
        if  settings.calcLandeJ
            println(stream, " ")
            println(stream, "  Lande g_J factors and Zeeman amplitudes:")
            println(stream, " ")
            println(stream, "  ", JAC.TableStrings.hLine(135))
            sa = "  ";   sb = "  "
            sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
            sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
            sa = sa * JAC.TableStrings.center(14, "Energy"; na=5)              
            sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=5)
            sa = sa * JAC.TableStrings.center(12, "Lande-J"; na=8);                           sb = sb * JAC.TableStrings.hBlank(20)          
            sa = sa * JAC.TableStrings.center(62, "N1   -- Zeeman Amplitudes --   Delta N1"    ; na=4) 
            sb = sb * JAC.TableStrings.center(62, "   Re              Im                 Re              Im"; na=5)
            println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(135)) 
            #  
            for  outcome in outcomes
                sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
                sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.Jlevel.index); na=2)
                sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
                energy = 1.0
                sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", energy))                    * "    "
                sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", outcome.LandeJ) )              * "    "
                sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeN1.re) )      * " "
                sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeN1.im) )      * "    "
                sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeDeltaN1.re) ) * " "
                sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", outcome.amplitudeDeltaN1.im) ) * "    "
                println(stream, sa )
            end
            println(stream, "  ", JAC.TableStrings.hLine(135))
        end
        #
        if  settings.calcLandeF
            println(stream, " ")
            println(stream, "  Hyperfine levels and Lande g_F factors:")
            println(stream, " ")
            println(stream, "  ", JAC.TableStrings.hLine(80))
            sa = "  ";   sb = "  "
            sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
            sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
            sa = sa * JAC.TableStrings.center(14, "Energy"; na=5)              
            sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=5)
            sa = sa * JAC.TableStrings.center(10, "F^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
            sa = sa * JAC.TableStrings.center(12, "Lande-F"; na=8);                           sb = sb * JAC.TableStrings.hBlank(20)          
            println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(80)) 
            #
            for  outcome in outcomes
                sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
                sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.Jlevel.index); na=2)
                sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
                energy = 1.0
                sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", energy))                    * "    "
                println(stream, sa )
                #
                sc = JAC.TableStrings.hBlank( length(sa) + 1 )
                Flist = JAC.oplus(nm.spinI, outcome.Jlevel.J)
                for  F in Flist
                    symf   = LevelSymmetry( F, outcome.Jlevel.parity);      Fx = JAC.AngularMomentum.oneJ(F)
                    Jx     = JAC.AngularMomentum.oneJ(outcome.Jlevel.J);    Ix = JAC.AngularMomentum.oneJ(nm.spinI)
                    LandeF = (Fx*(Fx+1) + Jx*(Jx+1) - Ix*(Ix+1)) / (2*Fx*(Fx+1)) * outcome.LandeJ
                    sa = sc * JAC.TableStrings.center(10, string(symf); na=4)
                    sa = sa * JAC.TableStrings.flushright(15, @sprintf("%.8e", LandeF) )   
                    println(stream, sa )
                end
            end
            println(stream, "  ", JAC.TableStrings.hLine(80))
        end
        #
        if  settings.calcZeeman
            println(stream, " ")
            println(stream, "  Zeeman splittings of (hyper-) fine-structure levels:")
            println(stream, " ")
            println(stream, "  ", JAC.TableStrings.hLine(135))
            println(stream, "  ", JAC.TableStrings.hLine(135))
        end
        #
        return( nothing )
    end

end # module
