
"""
`Basics.tabulate()`  ... tabulates the data from various objects due to different criteria.

  + `("multiplet: energies", multiplet::Multiplet; stream::IO=stdout)`  
                             ... to tabulate the energies of all levels of the given multiplet into a neat format; nothing is returned.
  + `("multiplet: energy relative to immediately lower level", multiplet::Multiplet; stream::IO=stdout)`  
                             ... to tabulate the energy splitting between neighboured levels of all levels of the given multiplet into a neat
                                 format; nothing is returned.
  + `("multiplet: energy of each level relative to lowest level", multiplet::Multiplet; stream::IO=stdout)`  
                             ... to tabulate the energy splitting of all levels with regard to the lowest level of the given multiplet into 
                                 a neat format; nothing is returned.
"""
function Basics.tabulate(sa::String, multiplet::Multiplet; stream::IO=stdout)
    if        sa == "multiplet: energies"
        println(stream, "\n  Eigenenergies:")
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " *  JAC.TableStrings.inUnits("energy")     
        println(stream, "\n", sb, "\n")
        for  i = 1:length(multiplet.levels)
            lev = multiplet.levels[i]
            en  = lev.energy;    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = " " * JAC.TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
            @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
        end

    elseif    sa == "multiplet: energy relative to immediately lower level"
        println(stream, "\n  Energy of each level relative to immediately lower level:")
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " * JAC.TableStrings.inUnits("energy")  
        println(stream, "\n", sb, "\n")
        for  i = 2:length(multiplet.levels)
            lev = multiplet.levels[i]
            en    = lev.energy - multiplet.levels[i-1].energy;    
            en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = " " * JAC.TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
            @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
        end

    elseif    sa == "multiplet: energy of each level relative to lowest level"
        println(stream, "\n  Energy of each level relative to lowest level:")
        sb = "  Level  J Parity          Hartrees       " * "             eV                   " * JAC.TableStrings.inUnits("energy")      
        println(stream, "\n", sb, "\n")
        for  i = 2:length(multiplet.levels)
            lev = multiplet.levels[i]
            en    = lev.energy - multiplet.levels[1].energy;    
            en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = " " * JAC.TableStrings.level(i) * "    " * string(LevelSymmetry(lev.J, lev.parity)) * "    "
            @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
        end
    else
        error("Unsupported keystring.")
    end

    return( nothing )  
end



"""
  + `("multiplet: energies", multiplet::JAC.Hfs.IJF_Multiplet; stream::IO=stdout)`  
                             ... to tabulate the energies of all hyperfine levels of the given multiplet into a neat format; nothing is returned.
  + `("multiplet: energy of each level relative to lowest level", multiplet::JAC.Hfs.IJF_Multiplet; stream::IO=stdout)`  
                             ... to tabulate the energy splitting of all levels with regard to the lowest level of the given multiplet into 
                                 a neat format; nothing is returned.
"""
function Basics.tabulate(sa::String, multiplet::JAC.Hfs.IJF_Multiplet; stream::IO=stdout)
    if        sa == "multiplet: energies"
        println(stream, "\n  Eigenenergies for nuclear spin I = $(multiplet.levelFs[1].I):")
        sb = "  Level  F Parity          Hartrees       " * "             eV                   " *  JAC.TableStrings.inUnits("energy")     
        println(stream, "\n", sb, "\n")
        for  i = 1:length(multiplet.levelFs)
            lev = multiplet.levelFs[i]
            en  = lev.energy;    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = " " * JAC.TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity)) * "    "
            @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
        end

    elseif    sa == "multiplet: energy of each level relative to lowest level"
        println(stream, "\n  Energy of each level relative to lowest level for nuclear spin I = $(multiplet.levelFs[1].I):")
        sb = "  Level  F Parity          Hartrees       " * "             eV                   " * JAC.TableStrings.inUnits("energy")      
        println(stream, "\n", sb, "\n")
        for  i = 2:length(multiplet.levelFs)
            lev = multiplet.levelFs[i]
            en    = lev.energy - multiplet.levelFs[1].energy;    
            en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = " " * JAC.TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity))  * "    "
            @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
        end
    else
        error("Unsupported keystring.")
    end

    return( nothing )  
end


"""
`Basics.tabulateKappaSymmetryEnergiesDirac(kappa::Int64, evalues::Array{Float64,1}, ns::Int64, nuclearModel::Nuclear.Model)`  ... tabulates the 
             eigenenergies for a given symmetry block kappa together with the corresponding Dirac energies for a point-like nucleus. The index 
             ns tells the number of 'negative-continnum' energies in the given evalues. nothing is returned.
"""
function Basics.tabulateKappaSymmetryEnergiesDirac(kappa::Int64, evalues::Array{Float64,1}, ns::Int64, nuclearModel::Nuclear.Model)
    Z = nuclearModel.Z
    # Determine the allowed principal quantum numbers n
    l = Basics.subshell_l(Subshell(101,kappa))
    println("  ", JAC.TableStrings.hLine(77))
    sa = "  "
    sa = sa * JAC.TableStrings.center( 7, "Index";            na=2)
    sa = sa * JAC.TableStrings.center(10, "Subshell";         na=3)
    sa = sa * JAC.TableStrings.center(17, "Energies [a.u.]";  na=2)
    sa = sa * JAC.TableStrings.center(17, "Dirac-E  [a.u.]";  na=2)
    sa = sa * JAC.TableStrings.center(17, "Delta-E / |E|";    na=2)
    println(sa)
    println("  ", JAC.TableStrings.hLine(77))
    for  i = ns+1:ns+8
        sa = " " * JAC.TableStrings.center( 6, JAC.TableStrings.level(i-ns); na=2)
        sa = sa *  JAC.TableStrings.flushright(10, string(Subshell(i-ns+l, kappa)); na=6)
        en = Basics.computeDiracEnergy(Subshell(i-ns+l, kappa), Z)
        if  evalues[i] >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", evalues[i])                       * "    "
        if  en         >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", en)                               * "    "
        if  evalues[i]-en >= 0.   sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", (evalues[i]-en)/abs(evalues[i])) * "    "
        println(sa)
    end
    println("      :       :    ")
    for  i = length(evalues)-5:length(evalues)
        sa = " " * JAC.TableStrings.center( 6, JAC.TableStrings.level(i-ns); na=2)
        sa = sa *  JAC.TableStrings.flushright(10, string(Subshell(i-ns+l, kappa)); na=6)
        en = Basics.computeDiracEnergy(Subshell(i-ns+l, kappa), Z)
        if  evalues[i] >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", evalues[i])                       * "    "
        if  en         >= 0.      sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", en)                               * "    "
        if  evalues[i]-en >= 0.   sa = sa * "+"   end;    sa = sa * @sprintf("%.8e", (evalues[i]-en)/abs(evalues[i])) * "    "
        println(sa)
    end
    println("  ", JAC.TableStrings.hLine(77))

    return( nothing )
end

