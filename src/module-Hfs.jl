
"""
`module  JAC.Hfs`  
... a submodel of JAC that contains all methods for computing HFS A and B coefficients, hyperfine sublevel representations, etc.
"""
module Hfs


using Printf, ..AngularMomentum, ..Basics,  ..Defaults, ..InteractionStrength, ..ManyElectron, ..Radial, ..Nuclear, 
                ..SpinAngular, ..TableStrings, ..PhotoEmission


"""
`struct  Hfs.InteractionMatrix`  ... defines a type for storing the T^1 and T^2 interaction matrices for a given basis.

    + calcT1   ::Bool               ... true, if the matrixT1 has been calculated and false otherwise.
    + calcT2   ::Bool               ... true, if the matrixT2 has been calculated and false otherwise.
    + matrixT1 ::Array{Float64,2}   ... T1 interaction matrix
    + matrixT2 ::Array{Float64,2}   ... T2 interaction matrix

"""
struct InteractionMatrix
    calcT1     ::Bool
    calcT2     ::Bool
    matrixT1   ::Array{Float64,2}
    matrixT2   ::Array{Float64,2}
end 


"""
`Hfs.InteractionMatrix()`  ... constructor for an `empty` instance of InteractionMatrix.
"""
function InteractionMatrix()
    InteractionMatrix(false, false, zeros(2,2), zeros(2,2))
end


# `Base.show(io::IO, im::Hfs.InteractionMatrix)`  ... prepares a proper printout of the variable InteractionMatrix.
function Base.show(io::IO, im::Hfs.InteractionMatrix) 
    println(io, "calcT1:           $(im.calcT1)  ")
    println(io, "calcT1:           $(im.calcT1)  ")
    println(io, "matrixT1:         $(im.matrixT1)  ")
    println(io, "matrixT2:         $(im.matrixT2)  ")
end


"""
`struct  Hfs.IJF_Vector`  ... defines a type for a IJF-coupled basis vector, here based on an ASF.

    + I        ::AngularJ64   ... Nuclear spin I.
    + F        ::AngularJ64   ... Total angular momentum F
    + levelJ   ::Level        ... Atomic level with well-defined total (electronic) angular momentum J

"""
struct IJF_Vector
    I          ::AngularJ64
    F          ::AngularJ64
    levelJ     ::Level
end 


"""
`Hfs.IJF_Vector()`  ... constructor for an `empty` instance of IJF_Vector.
"""
function IJF_Vector()
    IJF_Vector(AngularJ64(0), AngularJ64(0), Level())
end


# `Base.show(io::IO, ijfVector::Hfs.IJF_Vector)`  ... prepares a proper printout of the variable ijfVector.
function Base.show(io::IO, ijfVector::Hfs.IJF_Vector) 
    println(io, "I:           $(ijfVector.I)  ")
    println(io, "F:           $(ijfVector.F)  ")
    println(io, "levelJ:      $(ijfVector.levelJ)  ")
end


"""
`struct  Hfs.IJF_Basis`  ... defines a type for a IJF-coupled basis.

    + vectors   ::Array{IJF_Vector,1}   ... List of IJF_Vectors that form the basis.
    + basisJ    ::Basis                 ... Electronic basis that allows access to the electronic orbitals and CSF basis.
"""
struct IJF_Basis
    vectors     ::Array{IJF_Vector,1}
    basisJ      ::Basis
end 


"""
`Hfs.IJF_Basis()`  ... constructor for an `empty` instance of IJF_Basis.
"""
function IJF_Basis()
    IJF_Basis(IJF_Vector[], Basis())
end


# `Base.show(io::IO, ijfBasis::Hfs.IJF_Basis)`  ... prepares a proper printout of the variable ijfBasis.
function Base.show(io::IO, ijfBasis::Hfs.IJF_Basis) 
    println(io, "levelFs:           $(ijfBasis.levelFs)  ")
    println(io, "basisJ:            $(ijfBasis.basisJ)  ")
end


"""
`struct  Hfs.IJF_Level`  ... defines a type for a IJF-coupled level.

    + I              ::AngularJ64       ... Nuclear spin I.
    + F              ::AngularJ64       ... Total angular momentum F.
    + M              ::AngularM64       ... Total projection M, only defined if a particular sublevel is referred to.
    + parity         ::Parity           ... Parity of the level which corresponds to the electronic system.
    + energy         ::Float64          ... energy
    + basis          ::IJF_Basis        ... basis for this level
    + mc             ::Vector{Float64}  ... Vector of mixing coefficients w.r.t basis.
"""
struct IJF_Level
    I                ::AngularJ64
    F                ::AngularJ64
    M                ::AngularM64
    parity           ::Parity
    energy           ::Float64
    basis            ::IJF_Basis
    mc               ::Vector{Float64}
end 


"""
`Hfs.IJF_Level()`  ... constructor for an `empty` instance of IJF_Level.
"""
function IJF_Level()
    IJF_Level(AngularJ64(0), AngularJ64(0), AngularJ64(0), AngularM64(0), Basics.plus, 0., false, IJF_Basis(), Float64[])
end


# `Base.show(io::IO, ijfLevel::Hfs.IJF_Level)`  ... prepares a proper printout of the variable ijfLevel::Hfs.IJF_Level.
function Base.show(io::IO, ijfLevel::Hfs.IJF_Level) 
    println(io, "I:              $(ijfLevel.I)  ")
    println(io, "J:              $(ijfLevel.J)  ")
    println(io, "F:              $(ijfLevel.F)  ")
    println(io, "M:              $(ijfLevel.M)  ")
    println(io, "parity:         $(ijfLevel.parity)  ")
    println(io, "energy:         $(ijfLevel.energy)  ")
    println(io, "basis:          $(ijfLevel.basis)  ")
    println(io, "mc:             $(ijfLevel.mc)  ")
end


"""
`struct  Hfs.IJF_Multiplet`  ... defines a type for a multiplet of IJF-coupled levels.

    + name     ::String                ... A name associated to the multiplet.
    + levelFs  ::Array{IJF_Level,1}    ... List of IJF-coupled levels.

"""
struct IJF_Multiplet
    name       ::String
    levelFs    ::Array{IJF_Level,1}
end 


"""
`Hfs.IJF_Multiplet()`  ... constructor for an `empty` instance of Hfs.IJF_Multiplet.
"""
function IJF_Multiplet()
    IJF_Multiplet("", IJF_Level[])
end


# `Base.show(io::IO, ijfMultiplet::Hfs.IJF_Multiplet)`  ... prepares a proper printout of the variable ijfMultiplet::Hfs.IJF_Multiplet.
function Base.show(io::IO, ijfMultiplet::Hfs.IJF_Multiplet) 
    println(io, "name:           $(ijfMultiplet.name)  ")
    println(io, "levelFs:        $(ijfMultiplet.levelFs)  ")
end


"""
`struct  Hfs.Outcome`  
    ... defines a type to keep the outcome of a HFS computation, such as the HFS A and B coefficients as well 
        other results.

    + Jlevel                    ::Level            ... Atomic level to which the outcome refers to.
    + AIoverMu                  ::Float64          ... HFS A * I / mu value.
    + BoverQ                    ::Float64          ... HFS B / Q value
    + amplitudeT1               ::Complex{Float64} ... T1 amplitude
    + amplitudeT2               ::Complex{Float64} ... T2 amplitude
    + nuclearI                  ::AngularJ64       ... nuclear spin
    + hfsMultiplet              ::IJF_Multiplet    ... Multiplet of IJF levels as associated with the JLevel.
"""
struct Outcome 
    Jlevel                      ::Level 
    AIoverMu                    ::Float64
    BoverQ                      ::Float64
    amplitudeT1                 ::Complex{Float64}
    amplitudeT2                 ::Complex{Float64}
    nuclearI                    ::AngularJ64
    hfsMultiplet                ::IJF_Multiplet
end 


"""
`Hfs.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of HFS properties.
"""
function Outcome()
    Outcome(Level(), 0., 0., 0., 0., AngularJ64(0), IJF_Multiplet() )
end


# `Base.show(io::IO, outcome::Hfs.Outcome)`  ... prepares a proper printout of the variable outcome::Hfs.Outcome.
function Base.show(io::IO, outcome::Hfs.Outcome) 
    println(io, "Jlevel:                    $(outcome.Jlevel)  ")
    println(io, "AIoverMu:                  $(outcome.AIoverMu)  ")
    println(io, "BoverQ:                    $(outcome.BoverQ)  ")
    println(io, "amplitudeT1:               $(outcome.amplitudeT1)  ")
    println(io, "amplitudeT2:               $(outcome.amplitudeT2)  ")
    println(io, "nuclearI:                  $(outcome.nuclearI)  ")
    println(io, "hfsMultiplet:               (outcome.hfsMultiplet)  ")
end


"""
`struct  Settings  <:  AbstractPropertySettings`  ... defines a type for the details and parameters of computing HFS A and B coefficients.

    + calcT1                    ::Bool             ... True if T1-amplitudes (HFS A values) need to be calculated, and false otherwise.
    + calcT2                    ::Bool             ... True if T2-amplitudes (HFS B values) need to be calculated, and false otherwise.
    + calcNondiagonal           ::Bool             
        ... True if also (non-)diagonal hyperfine amplitudes are to be calculated and printed, and false otherwise.
    + calcIJFexpansion          ::Bool             
        ... True if the selected atomic levels are to be represented in a IJF-coupled basis, and false otherwise.
    + printBefore               ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
    + levelSelection            ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    calcT1                      ::Bool
    calcT2                      ::Bool
    calcNondiagonal             ::Bool 
    calcIJFexpansion            ::Bool 
    printBefore                 ::Bool 
    levelSelection              ::LevelSelection
end 


"""
`Hfs.Settings(; calcT1::Bool=true,` calcT2::Bool=false, calcNondiagonal::Bool=false, calcIJFexpansion::Bool=false, 
                    printBefore::Bool=false, levelSelection::LevelSelection=LevelSelection()) 
    ... keyword constructor to overwrite selected value of Einstein line computations.
"""
function Settings(; calcT1::Bool=true, calcT2::Bool=false, calcNondiagonal::Bool=false, calcIJFexpansion::Bool=false, 
                    printBefore::Bool=false, levelSelection::LevelSelection=LevelSelection())
    Settings(calcT1, calcT2, calcNondiagonal, calcIJFexpansion, printBefore, levelSelection)
end


# `Base.show(io::IO, settings::Hfs.Settings)`  ... prepares a proper printout of the variable settings::Hfs.Settings.
function Base.show(io::IO, settings::Hfs.Settings) 
    println(io, "calcT1:                   $(settings.calcT1)  ")
    println(io, "calcT2:                   $(settings.calcT2)  ")
    println(io, "calcNondiagonal:          $(settings.calcNondiagonal)  ")
    println(io, "calcIJFexpansion:         $(settings.calcIJFexpansion)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


function Base.isless(x::Hfs.IJF_Level, y::Hfs.IJF_Level)
    return x.energy < y.energy
end


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

    
"""
`Basics.sortByEnergy(multiplet::Hfs.IJF_Multiplet)`  
    ... to sort all hyperfine levels in the multiplet into a sequence of increasing energy; a (new) multiplet::Hfs.IJF_Multiplet 
        is returned.
"""
function Basics.sortByEnergy(multiplet::Hfs.IJF_Multiplet)
    sortedLevels = Base.sort( multiplet.levelFs , lt=Base.isless)
    newLevels = Hfs.IJF_Level[];   index = 0
    for lev in sortedLevels
        index = index + 1
        push!(newLevels, Hfs.IJF_Level(lev.I, lev.F, lev.M, lev.parity, lev.energy, lev.basis, lev.mc) )
    end
    
    newMultiplet = Hfs.IJF_Multiplet(multiplet.name, newLevels)
    
    return( newMultiplet )  
end


"""
`Basics.tabulate(sa::String, multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)`  
    ... tabulates the energies from the multiplet due to different criteria.

+ `("multiplet: energies", multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)`  
    ... to tabulate the energies of all hyperfine levels of the given multiplet into a neat format; nothing is returned.
+ `("multiplet: energy of each level relative to lowest level", multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)`  
    ... to tabulate the energy splitting of all levels with regard to the lowest level of the given multiplet into 
        a neat format; nothing is returned.
"""
function Basics.tabulate(sa::String, multiplet::Hfs.IJF_Multiplet; stream::IO=stdout)
    if        sa == "multiplet: energies"
        println(stream, "\n  Eigenenergies for nuclear spin I = $(multiplet.levelFs[1].I):")
        sb = "  Level  F Parity          Hartrees       " * "             eV                   " *  TableStrings.inUnits("energy")     
        println(stream, "\n", sb, "\n")
        for  i = 1:length(multiplet.levelFs)
            lev = multiplet.levelFs[i]
            en  = lev.energy;    en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity)) * "    "
            @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
        end

    elseif    sa == "multiplet: energy of each level relative to lowest level"
        println(stream, "\n  Energy of each level relative to lowest level for nuclear spin I = $(multiplet.levelFs[1].I):")
        sb = "  Level  F Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")      
        println(stream, "\n", sb, "\n")
        for  i = 2:length(multiplet.levelFs)
            lev = multiplet.levelFs[i]
            en    = lev.energy - multiplet.levelFs[1].energy;    
            en_eV = Defaults.convertUnits("energy: from atomic to eV", en);    en_requested = Defaults.convertUnits("energy: from atomic", en)
            sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity))  * "    "
            @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
        end
    else
        error("Unsupported keystring.")
    end

    return( nothing )  
end


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`Hfs.amplitude(kind::String, rLevel::Level, sLevel::Level, grid::Radial.Grid; printout::Bool=true)` 
    ... to compute either the  T^(1) or T^(2) hyperfine amplitude <alpha_r J_r || T^(n)) || alpha_s J_s>  
        for a given pair of levels. A value::ComplexF64 is returned.
"""
function  amplitude(kind::String, rLevel::Level, sLevel::Level, grid::Radial.Grid; printout::Bool=true)
    #
    if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
    nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
    if  printout   printstyled("Compute hyperfine $(kind[1:5]) matrix of dimension $nr x $ns ... \n", color=:light_green)   end
    #
    for  r = 1:nr
        for  s = 1:ns
            me = 0.
            if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                rLevel.parity != sLevel.parity    continue    
            end 
            #
            if      kind == "T^(1) amplitude"
            #--------------------------------
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
                    tamp  = InteractionStrength.hfs_tM1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                    #
                    ## println("**  <$(coeff.a) || t1 || $(coeff.b)>  = $(coeff.T * tamp)   = $(coeff.T) * $tamp" )
                    me = me + coeff.T * tamp  
                end
            #
            elseif  kind == "T^(2) amplitude"
            #--------------------------------
                # Calculate the spin-angular coefficients
                if  Defaults.saRatip()
                    waR = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 2, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    wa  = waR       
                end
                if  Defaults.saGG()
                    subshellList = sLevel.basis.subshells
                    opa = SpinAngular.OneParticleOperator(2, plus, true)
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
                    tamp  = InteractionStrength.hfs_tE2(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                    #
                    ## println("**  <$(coeff.a) || t2 || $(coeff.b)>  = $(coeff.T * tamp)   = $(coeff.T) * $tamp" )
                    me = me + coeff.T * tamp  
                end
            #
            else    error("stop a")
            end
            #
            matrix[r,s] = me
            println(">>> HF interaction matrix:  <$r || $kind || $s> = $me" )
        end
    end
    if  printout   printstyled("done.\n", color=:light_green)   end
    amplitude = transpose(rLevel.mc) * matrix * sLevel.mc 
    #
    return( amplitude )
end
"""
`Hfs.computeInteractionAmplitudeM(mp::EmMultipole, leftIsomer::Nuclear.Isomer, rightIsomer::Nuclear.Isomer)` 
    ... to compute the hyperfine interaction amplitude (<leftIsomer || M^(mp)) || rightIsomer>) for the interaction of two
        nuclear levels; this ME is geometrically fixed if the left and right isomer are the same, and it depends
        on the nuclear ME otherwise. An amplitude::ComplexF64 is returned.
"""
function  computeInteractionAmplitudeM(mp::EmMultipole, leftIsomer::Nuclear.Isomer, rightIsomer::Nuclear.Isomer)
    amplitude = 1.
    # Calculate the geometrical factor if the left- and right-hand isomer is the same
    if  leftIsomer == rightIsomer
        floatI = Basics.twice(leftIsomer.spinI) / 2.
        if       mp == M1       amplitude = leftIsomer.mu * sqrt( (floatI + 1)*(2*floatI+1) / floatI)
        elseif   mp == E2       amplitude = leftIsomer.Q / 2 * sqrt( (floatI + 1)*(2*floatI+1) * (2*floatI + 3)/ (floatI * (2*floatI -1)) )
        else   error("stop a; mp = $mp")
        end
    else
        if mp in leftIsomer.multipoleM && mp in rightIsomer.multipoleM
            lidx = findall(==(mp), leftIsomer.multipoleM)
            ridx = findall(==(mp), rightIsomer.multipoleM)
            if lidx != ridx
                error("stop a; leftIsomer.multipoleM != rightIsomer.multipoleM ") 
            else
                if length(lidx) == 1 && length(ridx) == 1 
                    amplitude = (leftIsomer.elementM[lidx[1]] + rightIsomer.elementM[ridx[1]]) / 2
                    if rightIsomer.energy < leftIsomer.energy
                        amplitude =amplitude *(-1)^(Basics.twice(rightIsomer.spinI)/2-Basics.twice(leftIsomer.spinI)/2)
                    end
                else error("stop b; leftIsomer.multipoleM setting error")
                end
            end
        else
            amplitude = 0.
        end  
    end  
    return( amplitude )
end

"""
`Hfs.computeInteractionAmplitudeT(mp::EmMultipole, aLevel::Level, bLevel, grid::Radial.Grid)` 
    ... to compute the T^(mp) interaction matrices for the given basis, i.e. (<aLevel || T^(mp) || bLevel>).
        Both levels must refer to the same basis. A me::ComplexF64 is returned.
"""
function  computeInteractionAmplitudeT(mp::EmMultipole, aLevel::Level, bLevel, grid::Radial.Grid)
    #
    ncsf = length(aLevel.basis.csfs);  me = ComplexF64(0.)
    if  ncsf != length(bLevel.basis.csfs)  ||  aLevel.basis.subshells != bLevel.basis.subshells
        error("stop a: both levels must refer to the same electronic basis.")
    end 
    
    # Compute the  T^(mp) matrix element
    for  (ia, csfa)  in  enumerate(aLevel.basis.csfs)
        for  (ib, csfb)  in  enumerate(bLevel.basis.csfs)
            wb  = ComplexF64(0.)
            if  abs(aLevel.mc[ia] * bLevel.mc[ib]) > 1.0e-10
                if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity   error("stop b")    end 
                subshellList = aLevel.basis.subshells
                orbitals     = aLevel.basis.orbitals
                opa = SpinAngular.OneParticleOperator(mp.L, plus, true)
                wa  = SpinAngular.computeCoefficients(opa, aLevel.basis.csfs[ia], bLevel.basis.csfs[ib], subshellList)
                for  coeff in wa
                    ja   = Basics.subshell_2j(orbitals[coeff.a].subshell)
                    jb   = Basics.subshell_2j(orbitals[coeff.b].subshell)
                    if      mp == M1    tamp = InteractionStrength.hfs_tM1(orbitals[coeff.a], orbitals[coeff.b], grid)
                    elseif  mp == E2    tamp = InteractionStrength.hfs_tE2(orbitals[coeff.a], orbitals[coeff.b], grid)
                    elseif  mp == E1    tamp = InteractionStrength.hfs_tE1(orbitals[coeff.a], orbitals[coeff.b], grid)  
                    elseif  mp == E3    tamp = InteractionStrength.hfs_tE3(orbitals[coeff.a], orbitals[coeff.b], grid)
                    elseif  mp == M2    tamp = InteractionStrength.hfs_tM2(orbitals[coeff.a], orbitals[coeff.b], grid)   
                    elseif  mp == M3    tamp = InteractionStrength.hfs_tM3(orbitals[coeff.a], orbitals[coeff.b], grid)   
                    else    error("stop b")    
                    end 
                    wb = wb + coeff.T * tamp/ sqrt( ja + 1) * sqrt( (Basics.twice(aLevel.J) + 1))    
                end
            end 
            me = me + aLevel.mc[ia] * bLevel.mc[ib] * wb    
        end 
    end 

    return( me )
end 

"""
`Hfs.computeAmplitudesProperties(outcome::Hfs.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings, im::Hfs.InteractionMatrix) 
    ... to compute all amplitudes and properties of for a given level; an outcome::Hfs.Outcome is returned for which the 
        amplitudes and properties are now evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::Hfs.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings, im::Hfs.InteractionMatrix)
    AIoverMu = BoverQ = amplitudeT1 = amplitudeT2 = 0.;    J = AngularMomentum.oneJ(outcome.Jlevel.J)
    if  settings.calcT1  &&  outcome.Jlevel.J != AngularJ64(0)
        if  im.calcT1   amplitudeT1 = transpose(outcome.Jlevel.mc) * im.matrixT1 * outcome.Jlevel.mc
        else            amplitudeT1 = Hfs.amplitude("T^(1) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
        end
        
        AIoverMu = amplitudeT1 / sqrt(J * (J+1))
    end
    #
    if  settings.calcT2  &&  outcome.Jlevel.J != AngularJ64(0)
        if  im.calcT2   amplitudeT2 = transpose(outcome.Jlevel.mc) * im.matrixT2 * outcome.Jlevel.mc
        else            amplitudeT2 = Hfs.amplitude("T^(2) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
        end
        
        BoverQ   = 2 * amplitudeT2 * sqrt( (2J-1) / ((J+1)*(2J+3)) )  ## * sqrt(J)
    end
    #
    hfsMultiplet = Hfs.IJF_Multiplet()
    if  settings.calcIJFexpansion
        # Determine a IJF_Multiplet for the given Jlevel/outcome
        hfsMultiplet = Hfs.computeHyperfineMultiplet(outcome.Jlevel, nm, grid)
    end
    newOutcome = Hfs.Outcome( outcome.Jlevel, AIoverMu, BoverQ, amplitudeT1, amplitudeT2, 
                              nm.spinI, hfsMultiplet)
    return( newOutcome )
end


"""
`Hfs.computeHyperfineMultiplet(level::Level, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... to compute a hyperfine multiplet, i.e. a representation of hyperfine levels within a hyperfine-coupled basis as defined by the
        given (electronic) level; a hfsMultiplet::IJF_Multiplet is returned.
"""
function computeHyperfineMultiplet(level::Level, nm::Nuclear.Model, grid::Radial.Grid)
    #
    hfsBasis     = Hfs.defineHyperfineBasis(level, nm)
    hfsMultiplet = Hfs.computeHyperfineRepresentation(hfsBasis, nm, grid)

    return( hfsMultiplet )
end


"""
`Hfs.computeHyperfineMultiplet(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)`  
    ... to compute a hyperfine multiplet, i.e. a representation of hyperfine levels within a hyperfine-coupled basis as defined by the
        given (electronic) multiplet; a hfsMultiplet::IJF_Multiplet is returned.
"""
function computeHyperfineMultiplet(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)
    println("")
    printstyled("Hfs.computeHyperfineMultiplet(): The computation of the hyperfine multiplet starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------ \n", color=:light_green)
    #
    hfsBasis     = Hfs.defineHyperfineBasis(multiplet, nm)
    hfsMultiplet = Hfs.computeHyperfineRepresentation(hfsBasis, nm, grid)
    # Print all results to screen
    hfsMultiplet = Basics.sortByEnergy(hfsMultiplet)
    Basics.tabulate("multiplet: energies",                                      hfsMultiplet) 
    Basics.tabulate("multiplet: energy of each level relative to lowest level", hfsMultiplet) 
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary     
        Basics.tabulate("multiplet: energies",                                      hfsMultiplet, stream=iostream) 
        Basics.tabulate("multiplet: energy of each level relative to lowest level", hfsMultiplet, stream=iostream) 
    end
    #
    if    output    return( hfsMultiplet )
    else            return( nothing )
    end
end


"""
`Hfs.computeHyperfineRepresentation(hfsBasis::IJF_Basis, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... to set-up and diagonalized the Hamiltonian matrix of H^(DFB) + H^(hfs) within the atomic hyperfine (IJF-coupled) basis;
        a hfsMultiplet::IJF_Multiplet is returned.
"""
function computeHyperfineRepresentation(hfsBasis::IJF_Basis, nm::Nuclear.Model, grid::Radial.Grid)
    n = length(hfsBasis.vectors);   matrix = zeros(n,n)
    for  r = 1:n
        for  s = 1:n
            matrix[r,s] = 0.
            if  r == s    matrix[r,s] = matrix[r,s] + hfsBasis.vectors[r].levelJ.energy                  end
            if  hfsBasis.vectors[r].F              !=  hfsBasis.vectors[s].F                continue     end
            if  hfsBasis.vectors[r].levelJ.parity  !=  hfsBasis.vectors[s].levelJ.parity    continue     end
            # Now add the hyperfine matrix element for K = 1 (magnetic) and K = 2 (electric) contributions
            spinI = nm.spinI;   Jr = hfsBasis.vectors[r].levelJ.J;   Js = hfsBasis.vectors[s].levelJ.J;   F = hfsBasis.vectors[r].F
            Ix    = AngularMomentum.oneJ(spinI)
            #
            wa = AngularMomentum.phaseFactor([spinI, +1, Jr, +1, F])  * 
                    AngularMomentum.Wigner_6j(spinI, Jr, F, Js, spinI, AngularJ64(1))
            wb = Hfs.amplitude("T^(1) amplitude", hfsBasis.vectors[r].levelJ, hfsBasis.vectors[s].levelJ, grid::Radial.Grid)
            wc = nm.mu * sqrt( (Ix+1)/Ix )
            matrix[r,s] = matrix[r,s] + wa * wb * wc     * 1.0e-6 ## fudge-factor to keep HFS interaction small
            #
            if  spinI  in [AngularJ64(0), AngularJ64(1//2)]                                 continue     end
            wa = AngularMomentum.phaseFactor([spinI, +1, Jr, +1, F]) * 
                    AngularMomentum.Wigner_6j(spinI, Jr, F, Js, spinI, AngularJ64(2))
            wb = Hfs.amplitude("T^(2) amplitude", hfsBasis.vectors[r].levelJ, hfsBasis.vectors[s].levelJ, grid::Radial.Grid)
            wc = nm.Q / 2. * sqrt( (Ix+1)*(2Ix+3)/ (Ix*(2Ix-1)) )
            matrix[r,s] = matrix[r,s] + wa * wb * wc     * 1.0e-6 ## fudge-factor to keep HFS interaction small
        end
    end
    #
    # Diagonalize the matrix and set-up the representation
    eigen    = Basics.diagonalize("matrix: LinearAlgebra", matrix)
    levelFs  = Hfs.IJF_Level[]
    for  ev = 1:length(eigen.values)
        # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
        evector   = eigen.vectors[ev];    en = eigen.values[ev]
        parity    = Basics.plus;    F = AngularJ64(0);     MF = AngularM64(0)
        for  r = 1:length(hfsBasis.vectors)
            if  abs(evector[r]) > 1.0e-6    
                parity = hfsBasis.vectors[r].levelJ.parity
                F      = hfsBasis.vectors[r].F;     MF = AngularM64(hfsBasis.vectors[r].F);   break    end
        end
        newlevelF = Hfs.IJF_Level(nm.spinI, F, MF, parity, en, hfsBasis, evector) 
        push!( levelFs, newlevelF)
    end
    hfsMultiplet = Hfs.IJF_Multiplet("hyperfine", levelFs)
    
    return( hfsMultiplet )
end


#==
"""
`Hfs.computeInteractionAmplitudeM(mp::EmMultipole, leftIsomer::Nuclear.Isomer, rightIsomer::Nuclear.Isomer)` 
    ... to compute the hyperfine interaction amplitude (<leftIsomer || M^(mp)) || rightIsomer>) for the interaction of two
        nuclear levels; this ME is geometrically fixed if the left and right isomer are the same, and it depends
        on the nuclear ME otherwise. An amplitude::ComplexF64 is returned.
"""
function  computeInteractionAmplitudeM(mp::EmMultipole, leftIsomer::Nuclear.Isomer, rightIsomer::Nuclear.Isomer)
    amplitude = 1.
    # Calculate the geometrical factor if the left- and right-hand isomer is the same
    if  leftIsomer == rightIsomer
        floatI = Basics.twice(leftIsomer.spinI) / 2.
        if       mp == M1       amplitude = leftIsomer.mu * sqrt( (floatI + 1) / floatI)
        elseif   mp == E2       amplitude = leftIsomer.Q / 2 * sqrt( (floatI + 1) * (2*floatI + 3)/ (floatI * (2*floatI -1)) )
        else   error("stop a; mp = $mp")
        end
    else
        if       leftIsomer.multipoleM != rightIsomer.multipoleM     
               error("stop a; $(leftIsomer.multipoleM) != $(rightIsomer.multipoleM) ")    
        elseif   mp == leftIsomer.multipoleM     amplitude = (leftIsomer.elementM + rightIsomer.elementM) / 2
        else                                     amplitude = 0.
        end
    end
    
    return( amplitude )
end


"""
`Hfs.computeInteractionAmplitudeT(mp::EmMultipole, aLevel::Level, bLevel, grid::Radial.Grid)` 
    ... to compute the T^(mp) interaction matrices for the given basis, i.e. (<aLevel || T^(mp) || bLevel>).
        Both levels must refer to the same basis. A me::ComplexF64 is returned.
"""
function  computeInteractionAmplitudeT(mp::EmMultipole, aLevel::Level, bLevel, grid::Radial.Grid)
    #
    ncsf = length(aLevel.basis.csfs);  me = ComplexF64(0.)
    if  ncsf != length(bLevel.basis.csfs)  ||  aLevel.basis.subshells != bLevel.basis.subshells
        error("stop a: both levels must refer to the same electronic basis.")
    end 
    
    # Compute the  T^(mp) matrix element
    for  (ia, csfa)  in  enumerate(aLevel.basis.csfs)
        for  (ib, csfb)  in  enumerate(bLevel.basis.csfs)
            wb  = ComplexF64(0.)
            if  abs(aLevel.mc[ia] * bLevel.mc[ib]) > 1.0e-10
                if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity   error("stop b")    end 
                subshellList = aLevel.basis.subshells
                orbitals     = aLevel.basis.orbitals
                opa = SpinAngular.OneParticleOperator(mp.L, plus, true)
                wa  = SpinAngular.computeCoefficients(opa, aLevel.basis.csfs[ia], bLevel.basis.csfs[ib], subshellList)
                for  coeff in wa
                    ja   = Basics.subshell_2j(orbitals[coeff.a].subshell)
                    jb   = Basics.subshell_2j(orbitals[coeff.b].subshell)
                    if      mp == M1    tamp = InteractionStrength.hfs_tM1(orbitals[coeff.a], orbitals[coeff.b], grid)
                    elseif  mp == E2    tamp = InteractionStrength.hfs_tE2(orbitals[coeff.a], orbitals[coeff.b], grid)
                    else    error("stop b")    
                    end 
                    wb = wb + coeff.T * tamp
                end
            end 
            me = me + aLevel.mc[ia] * bLevel.mc[ib] * wb    
        end 
    end 

    return( me )
end 
==#



"""
`Hfs.computeInteractionMatrix(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)` 
    ... to compute the T^1 and/or T^2 interaction matrices for the given basis, i.e. (<csf_r || T^(n)) || csf_s>).
        An im::Hfs.InteractionMatrix is returned.
"""
function  computeInteractionMatrix(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)
    #
    ncsf = length(basis.csfs);    matrixT1 = zeros(ncsf,ncsf);    matrixT2 = zeros(ncsf,ncsf)
    #
    if  settings.calcT1
        calcT1 = true;    matrixT1 = zeros(ncsf,ncsf)
        for  r = 1:ncsf
            for  s = 1:ncsf
                if  basis.csfs[r].parity  != basis.csfs[s].parity   continue    end 
                #
                # Calculate the spin-angular coefficients
                if  Defaults.saRatip()
                    waR = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 1, basis.csfs[r], basis.csfs[s])
                    wa  = waR       
                end
                if  Defaults.saGG()
                    subshellList = basis.subshells
                    opa = SpinAngular.OneParticleOperator(1, plus, true)
                    waG = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList) 
                    wa  = waG
                end
                if  Defaults.saRatip() && Defaults.saGG() && true
                    if  length(waR) != 0     println("\n>> Angular coeffients from GRASP/MCT   = $waR ")    end
                    if  length(waG) != 0     println(  ">> Angular coeffients from SpinAngular = $waG ")    end
                end
                #
                for  coeff in wa
                    ja   = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                    jb   = Basics.subshell_2j(basis.orbitals[coeff.b].subshell)
                    tamp = InteractionStrength.hfs_tM1(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                    matrixT1[r,s] = matrixT1[r,s] + coeff.T * tamp  
                end
            end
        end
    else   
        calcT1 = false;    matrixT1 = zeros(2,2)
    end
    #
    if  settings.calcT2
        calcT2 = true;    matrixT2 = zeros(ncsf,ncsf)
        for  r = 1:ncsf
            for  s = 1:ncsf
                if  basis.csfs[r].parity  != basis.csfs[s].parity   continue    end 
                #
                # Calculate the spin-angular coefficients
                if  Defaults.saRatip()
                    waR = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 2, basis.csfs[r], basis.csfs[s])
                    wa  = waR       
                end
                if  Defaults.saGG()
                    subshellList = basis.subshells
                    opa = SpinAngular.OneParticleOperator(2, plus, true)
                    waG = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList) 
                    wa  = waG
                end
                if  Defaults.saRatip() && Defaults.saGG() && true
                    if  length(waR) != 0     println("\n>> Angular coeffients from GRASP/MCT   = $waR ")    end
                    if  length(waG) != 0     println(  ">> Angular coeffients from SpinAngular = $waG ")    end
                end
                #
                for  coeff in wa
                    ja   = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                    jb   = Basics.subshell_2j(basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.hfs_tE2(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                    matrixT2[r,s] = matrixT2[r,s] + coeff.T * tamp  
                end
            end
        end
    else   
        calcT2 = false;    matrixT2 = zeros(2,2)
    end
    #
    im = Hfs.InteractionMatrix(calcT1, calcT2, matrixT1, matrixT2)
    #
    return( im )
end


"""
`Hfs.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)`  
    ... to compute (as selected) the HFS A and B parameters as well as hyperfine energy splittings for the levels 
        of the given multiplet and as specified by the given settings. The results are printed in neat tables to 
        screen and, if requested, an arrays{Hfs.Outcome,1} with all the results are returned.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)
    println("")
    printstyled("Hfs.computeOutcomes(): The computation of the Hyperfine amplitudes and parameters starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------ \n", color=:light_green)
    println("")
    outcomes = Hfs.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    Hfs.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    im = Hfs.computeInteractionMatrix(multiplet.levels[1].basis, grid, settings)
    newOutcomes = Hfs.Outcome[]
    for  outcome in outcomes
        newOutcome = Hfs.computeAmplitudesProperties(outcome, nm, grid, settings, im) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    Hfs.displayResults(stdout, newOutcomes, nm, settings)
    # Compute and display the non-diagonal hyperfine amplitudes, if requested
    if  settings.calcNondiagonal    Hfs.displayNondiagonal(stdout, multiplet, grid, settings)   end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    
        Hfs.displayResults(iostream, newOutcomes, nm, settings) 
        if  settings.calcNondiagonal    Hfs.displayNondiagonal(iostream, multiplet, grid, settings)   end
    end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`Hfs.defineHyperfineBasis(level::Level, nm::Nuclear.Model)`  
    ... to define/set-up an atomic hyperfine (IJF-coupled) basis for the given electronic level; 
        a hfsBasis::IJF_Basis is returned.
"""
function  defineHyperfineBasis(level::Level, nm::Nuclear.Model)
    function  display_ijfVector(i::Int64, vector::Hfs.IJF_Vector) 
        si = string(i);   ni = length(si);    sa = repeat(" ", 5);    sa = sa[1:5-ni] * si * ")  "
        sa = sa * "[" * string( LevelSymmetry(vector.levelJ.J, vector.levelJ.parity) ) * "] " * string(vector.F) * repeat(" ", 4)
        return( sa )
    end

    vectors = IJF_Vector[]
    Flist = Basics.oplus(nm.spinI, level.J)
    for  F in Flist     push!(vectors,  IJF_Vector(nm.spinI, F, level) )      end
    #
    hfsBasis = IJF_Basis(vectors, level.basis)
    
    return( hfsBasis )
end


"""
`Hfs.defineHyperfineBasis(multiplet::Multiplet, nm::Nuclear.Model)`  
    ... to define/set-up an atomic hyperfine (IJF-coupled) basis for the given electronic multipet; 
        a hfsBasis::IJF_Basis is returned.
"""
function  defineHyperfineBasis(multiplet::Multiplet, nm::Nuclear.Model)
    function  display_ijfVector(i::Int64, vector::Hfs.IJF_Vector) 
        si = string(i);   ni = length(si);    sa = repeat(" ", 5);    sa = sa[1:5-ni] * si * ")  "
        sa = sa * "[" * string( LevelSymmetry(vector.levelJ.J, vector.levelJ.parity) ) * "] " * string(vector.F) * repeat(" ", 4)
        return( sa )
    end


    vectors = IJF_Vector[]
    #
    for i = 1:length(multiplet.levels)
        Flist = Basics.oplus(nm.spinI, multiplet.levels[i].J)
        for  F in Flist
            push!(vectors,  IJF_Vector(nm.spinI, F, multiplet.levels[i]) )
        end
    end
    #
    hfsBasis = IJF_Basis(vectors, multiplet.levels[1].basis)
    #
    println(" ")
    println("  Construction of a atomic hyperfine (IJF-coupled) basis of dimension $(length(hfsBasis.vectors)), based on " *
            "$(length(multiplet.levels)) electronic levels and nuclear spin $(hfsBasis.vectors[1].I).")
    println("  Basis vectors are defined as [J^P] F: \n")
    sa = "  "
    for  k = 1:length(hfsBasis.vectors)
        if  length(sa) > 100    println(sa);   sa = "  "   end
        sa = sa * display_ijfVector(k, hfsBasis.vectors[k])
    end
    if   length(sa) > 5    println(sa)   end
    println(" ")
    
    return( hfsBasis )
end


"""
`Hfs.computeModifiedEinsteinRates(upperOutcome::Outcome, lowerOutcome::Outcome, multipoles::Array{EmMultipole,1}, gauge::EmGauge,
                                  grid::Radial.Grid)`  
    ... to compute and tabulate the modified Einstein amplitudes and rates for the hyperfine-resolved 
        transitions between the upper and lower outcome. The procedures assumes that the two outcomes provide
        a proper IJF expansion (multiplet) of the hyperfine levels of interest.
        A neat table is printed but nothing is returned otherwise
"""
function  computeModifiedEinsteinRates(upperOutcome::Outcome, lowerOutcome::Outcome, multipoles::Array{EmMultipole,1}, gauge::EmGauge,
                                       grid::Radial.Grid)
    if  upperOutcome.nuclearI != lowerOutcome.nuclearI   
            error("Inconsistent nuclear spins; upper-I=$(upperOutcome.nuclearI)  !=  lower-I=$(lowerOutcome.nuclearI)")
    end
    
    stream      = stdout
    amplitudesJ = ComplexF64[]
    iJsym = LevelSymmetry(upperOutcome.Jlevel.J, upperOutcome.Jlevel.parity)
    fJsym = LevelSymmetry(lowerOutcome.Jlevel.J, lowerOutcome.Jlevel.parity)
    omega = upperOutcome.Jlevel.energy - lowerOutcome.Jlevel.energy
    
    # First compute the multipole transition amplitudes between the J-levels of the upper and lower outcome
    println(stream, " ")
    println(stream, "  Multipole amplitudes for J-levels with symmetry $iJsym --> $fJsym:")
    println(stream, " ")
    for  multipole in multipoles
        if  multipole  in  [M1, M2, M3]    gaugex = Basics.Magnetic    else    gaugex = gauge  end
        ampJ = PhotoEmission.amplitude("emission", multipole, gaugex, omega, upperOutcome.Jlevel, lowerOutcome.Jlevel, grid)
        push!(amplitudesJ, ampJ)
        println(stream, "  J-level amplitude for $multipole  = $ampJ  ")
    end
    
    # Print a table header
    nx = 150
    println(stream, " ")
    println(stream, "  HFS modified Einstein amplitudes and rates:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = " ";   sb = "  "
    sa = sa * TableStrings.center(14, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(14, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(12, "i--F--f";   na=2);                         sb = sb * TableStrings.hBlank(15)
    sa = sa * TableStrings.center(12, "Energy"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center( 9, "Multipole"; na=1);                         sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(11, "Gauge"    ; na=4);                         sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(26, "A--Einstein--B"; na=3);       
    sb = sb * TableStrings.center(26, TableStrings.inUnits("rate")*"           "*TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.center(26, "re-- <Ff |amplitude L| Fi> --im"; na=2);       
    
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  upperLevF in  upperOutcome.hfsMultiplet.levelFs
        for  lowerLevF in  lowerOutcome.hfsMultiplet.levelFs
            for  (m, mp)  in  enumerate(multipoles)
                Ji = upperOutcome.Jlevel.J;   Fi = upperLevF.F
                Jf = lowerOutcome.Jlevel.J;   Ff = lowerLevF.F
                
                ampF = ComplexF64(0.)
                for  umc in upperLevF.mc,   lmc in lowerLevF.mc
                    ampF = ampF  +  umc * lmc * sqrt( AngularMomentum.bracket([Fi, Ff]) ) *
                           AngularMomentum.phaseFactor([Ji, +1, upperOutcome.nuclearI, +1, Ff, +1, AngularJ64(mp.L)]) * 
                           AngularMomentum.AngularMomentum.Wigner_6j(Fi, Ff, mp.L, Jf, Ji, upperOutcome.nuclearI) * amplitudesJ[m]
                end
                sa = ""
                sa = sa * TableStrings.center(14, TableStrings.levels_if(upperOutcome.Jlevel.index, lowerOutcome.Jlevel.index); na=2)
                sa = sa * TableStrings.center(14, TableStrings.symmetries_if(iJsym, iJsym); na=0)
                sc = "         " * string(Fi) * "    " * string(Ff)
                sa = sa * TableStrings.center(12, sc[end-9:end]; na=3)
                en = upperOutcome.Jlevel.energy - lowerOutcome.Jlevel.energy
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en)) * "    "
                sa = sa * TableStrings.center(9,  string(mp); na=4)
                if  mp  in  [M1, M2, M3]    gaugex = Basics.Magnetic    else    gaugex = gauge  end
                sa = sa * TableStrings.flushleft(11, string(gaugex);  na=0)
                chRate =  8pi * Defaults.getDefaults("alpha") * en / (Basics.twice(Ji) + 1) * (abs(ampF)^2) * (Basics.twice(Jf) + 1)
                ## sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein A",  line, chRate)) * "  "
                ## sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein B",  line, chRate)) * "    "
                sa = sa * @sprintf("% .6e", chRate)  * "  " 
                sa = sa * @sprintf("% .6e", chRate)  * "    " 
                sa = sa * @sprintf("% .6e", ampF.re) * "  " 
                sa = sa * @sprintf("% .6e", ampF.im) * "  " 
                println(stream, sa)
            end
            println(stream, " ")
        end
    end
    println(stream, "  ", TableStrings.hLine(nx)) 
    
    return( nothing )
end


"""
`Hfs.determineOutcomes(multiplet::Multiplet, settings::Hfs.Settings)`  
    ... to determine a list of Outcomes's for the computation of HFS A- and B-parameters for the given multiplet. 
        It takes into account the particular selections and settings. An Array{Hfs.Outcome,1} is returned. Apart from the 
        level specification, all physical properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::Hfs.Settings) 
    outcomes = Hfs.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, Hfs.Outcome(level, 0., 0., 0., 0., AngularJ64(0), Hfs.IJF_Multiplet() ) )
        end
    end
    return( outcomes )
end


"""
`Hfs.displayNondiagonal(stream::IO, multiplet::Multiplet, grid::Radial.Grid, settings::Hfs.Settings)`  
    ... to compute and display all non-diagonal hyperfine amplitudes for the selected levels. A small neat table of 
        all (pairwise) hyperfine amplitudes is printed but nothing is returned otherwise.
"""
function  displayNondiagonal(stream::IO, multiplet::Multiplet, grid::Radial.Grid, settings::Hfs.Settings)
    # Determine pairs to be calculated
    pairs = Tuple{Int64,Int64}[]
    for  (f, fLevel)  in  enumerate(multiplet.levels)
        for  (i, iLevel)  in  enumerate(multiplet.levels)
            if  Basics.selectLevel(fLevel, settings.levelSelection)  &&   Basics.selectLevel(iLevel, settings.levelSelection)
                push!( pairs, (f,i) )
            end
        end
    end
    #
    nx = 107
    println(stream, " ")
    println(stream, "  Selected (non-) diagonal hyperfine amplitudes:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(10, "Level_f"; na=2)
    sa = sa * TableStrings.center(10, "Level_i"; na=2)
    sa = sa * TableStrings.center(10, "J^P_f";   na=3)
    sa = sa * TableStrings.center(57, "T1   --   Amplitudes   --   T2"; na=4);              
    sa = sa * TableStrings.center(10, "J^P_f";   na=4)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  (f,i) in pairs
        sa   = "  ";    
        sa   = sa * TableStrings.center(10, string(f); na=2)
        sa   = sa * TableStrings.center(10, string(i); na=2)
        symf = LevelSymmetry( multiplet.levels[f].J, multiplet.levels[f].parity)
        symi = LevelSymmetry( multiplet.levels[i].J, multiplet.levels[i].parity)
        sa   = sa * TableStrings.center(10, string(symf); na=4)
        T1   = Hfs.amplitude("T^(1) amplitude", multiplet.levels[f], multiplet.levels[i], grid, printout=false)
        T2   = Hfs.amplitude("T^(2) amplitude", multiplet.levels[f], multiplet.levels[i], grid, printout=false)
        sa   = sa * @sprintf("%.5e %s %.5e", T1.re, "  ", T1.im) * "    "
        sa   = sa * @sprintf("%.5e %s %.5e", T2.re, "  ", T2.im) * "    "
        sa   = sa * TableStrings.center(10, string(symi); na=4)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`Hfs.displayOutcomes(outcomes::Array{Hfs.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations A small neat table of all 
        selected levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{Hfs.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected HFS levels:")
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
    println(" ")
    #
    return( nothing )
end


"""
`Hfs.displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1}, nm::Nuclear.Model, settings::Hfs.Settings)`  
    ... to display the energies, A- and B-values, Delta E_F energy shifts, etc. for the selected levels. All nuclear 
        parameters are taken from the nuclear model. A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1}, nm::Nuclear.Model, settings::Hfs.Settings)
    nx = 117
    println(stream, " ")
    println(stream, "  HFS parameters and amplitudes:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(31, "A/mu [mu_N]  --  HFS  --  B/Q [barn]"; na=4);              
    sb = sb * TableStrings.center(31, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(32, "T1 -- Amplitudes -- T2"    ; na=4);        sb = sb * TableStrings.hBlank(36)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = outcome.Jlevel.energy
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))           * "      "
        wa = Defaults.convertUnits("energy: from atomic", outcome.AIoverMu) / nm.spinI.num * nm.spinI.den 
        wa = wa / (2 * 1836.15267 ) * Defaults.getDefaults("alpha")  # take 2m_p and alpha into account
        sa = sa * @sprintf("%.8e", wa)                                                   * "    " 
        wa = Defaults.convertUnits("energy: from atomic", outcome.BoverQ) 
        wa = wa / Defaults.convertUnits("cross section: from atomic to barn", 1.0)  # take Q [barn] into account
        sa = sa * @sprintf("%.8e", wa)                                                   * "    "
        sa = sa * @sprintf("%.8e %s %.8e", outcome.amplitudeT1.re, "  ", outcome.amplitudeT2.re) * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    # Printout the Delta E_F energy shifts of the hyperfine levels |alpha F> with regard to the (electronic) levels |alpha J>
    nx = 90
    println(stream, " ")
    println(stream, "  HFS Delta E_F energy shifts with regard to the (electronic) level energies E_J:")
    println(stream, " ")
    println(stream, "    Nuclear spin I:                          $(nm.spinI) ")
    println(stream, "    Nuclear magnetic-dipole moment      mu = $(nm.mu)    ")
    println(stream, "    Nuclear electric-quadrupole moment   Q = $(nm.Q)     ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);                            sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(10, "F^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Delta E_F"; na=4);                         
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "C factor"; na=4);                         
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = outcome.Jlevel.energy
        #
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))           * "    "
        J     = AngularMomentum.oneJ(outcome.Jlevel.J);   spinI = AngularMomentum.oneJ(nm.spinI)
        Flist = Basics.oplus(nm.spinI, outcome.Jlevel.J)
        first = true
        for  Fang in Flist
            Fsym    = LevelSymmetry(Fang, outcome.Jlevel.parity)
            F       = AngularMomentum.oneJ(Fang)
            Cfactor = F*(F+1) - J*(J+1) - spinI*(spinI+1)
            energy  = outcome.AIoverMu * nm.mu / spinI * Cfactor / 2.
            if  abs(outcome.BoverQ) > 1.0e-10
                energy  = energy +  outcome.BoverQ * nm.Q * 3/4 * (Cfactor*(Cfactor+1) - spinI*(spinI+1)*J*(J+1) ) /
                                    ( 2spinI*(2spinI-1)*J*(2J-1) )
            end
            sb = TableStrings.center(10, string(Fsym); na=2)
            sb = sb * TableStrings.flushright(16, @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))) * "    "
            sb = sb * TableStrings.flushright(12, @sprintf("%.5e", Cfactor))
            #
            if   first    println(stream,  sa*sb );   first = false
            else          println(stream,  TableStrings.hBlank( length(sa) ) * sb )
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    if  settings.calcIJFexpansion
        nx = 90
        println(stream, " ")
        println(stream, "  IJF-coupled hyperfine levels:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2)
        sa = sa * TableStrings.center(10, "J^P";   na=4)
        sa = sa * TableStrings.center( 6, "F^P";   na=4);   na = length(sa)   
        sa = sa * TableStrings.flushleft(34, "IJF basis (1:)"; na=4)
        println(stream, sa);    sa = repeat(" ", na);   
        sa = sa * TableStrings.flushleft(34, "Mixing coefficients (1:)"; na=4);                         
        println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=6);    na = length(sa)
            sa = sa * "        "
            if   length(outcome.hfsMultiplet.levelFs) == 0  sa = sa * "No IJF levels";   println(stream, sa);    continue    end
            nvecs = min( length(outcome.hfsMultiplet.levelFs[1].basis.vectors), 6)
            for  nvec = 1:nvecs
                sa = sa * "[" * string(outcome.hfsMultiplet.levelFs[1].basis.vectors[nvec].levelJ.J) * "] " *
                                string(outcome.hfsMultiplet.levelFs[1].basis.vectors[nvec].F) * "       "
            end
            println(stream, sa)
            for  levelF in outcome.hfsMultiplet.levelFs
                sa = repeat(" ", na)
                sa = sa * string(levelF.F) * "       "
                for  nvec = 1:nvecs
                    sa = sa * @sprintf("% .5e", levelF.mc[nvec] ) * "  "
                end
                println(stream, sa)
            end
            println(stream, " ")
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end

end # module


