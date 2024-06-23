

#################################################################################################################################
#################################################################################################################################


"""
`@enum   Basics.EmGauge`  
    ... defines a enumeration for the (allowed) gauges to deal with the radiation field; this differs from the UseGauge which can 
        be selected for explicit computations.
                    
        + NoGauge      ... No gauge define (yet).
        + Coulomb      ... Coulomb gauge (= velocity gauge in the non-relativistic limit; Grant, 1974).
        + Babushkin    ... Babushkin gauge (= length gauge in the non-relativistic limit; Grant, 1974). 
        + Magnetic     ... Magnetic gauge (transverse velocity gauge).
        + Velocity     ... Transverse velocity gauge (Johnson, 2007).
        + Length       ... Length gauge (Johnson, 2007).
"""
@enum   EmGauge    NoGauge    Coulomb   Babushkin   Magnetic   Velocity   Length
export  EmGauge,   NoGauge,   Coulomb,  Babushkin,  Magnetic,  Velocity,  Length


"""
`Basics.EmGauge(sa::String)`  ... constructor for a given String.
"""
function EmGauge(sa::String)
    if       sa in ["none", "NoGauge"]         wa = NoGauge
    elseif   sa == "Coulomb"                   wa = Coulomb
    elseif   sa == "Babushkin"                 wa = Babushkin
    elseif   sa in ["Magnetic", "magnetic"]    wa = Magnetic
    else     error("stop a")
    end

    return( wa )
end


# `Base.string(gauge::EmGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
function Base.string(gauge::EmGauge) 
    if      gauge == NoGauge     return( "no gauge" )
    elseif  gauge == Coulomb     return( "Coulomb" )   
    elseif  gauge == Babushkin   return( "Babushkin" )   
    elseif  gauge == Magnetic    return( "Magnetic" )   
    else    error("stop a")
    end
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.EmMultipole`  ... defines a struct for the em multipoles, E1, M1, E2, ...

    + L           ::Int64    ... multipolarity
    + electric    ::Bool     ... True if the multipole is an electric one, and false otherwise.
"""  
struct  EmMultipole
    L             ::Int64
    electric      ::Bool   
end

export  EmMultipole, E1, M1, E2, M2, E3, M3, E4, M4


const  E1 = EmMultipole(1, true);     const  M1 = EmMultipole(1, false)
const  E2 = EmMultipole(2, true);     const  M2 = EmMultipole(2, false)
const  E3 = EmMultipole(3, true);     const  M3 = EmMultipole(3, false)
const  E4 = EmMultipole(4, true);     const  M4 = EmMultipole(4, false)


"""
`Basics.EmMultipole(sa::String)`  ... constructor for a given String.
"""
function EmMultipole(sa::String)
    if       sa == "E1"    L = 1;   electric = true
    elseif   sa == "M1"    L = 1;   electric = false
    elseif   sa == "E2"    L = 2;   electric = true
    elseif   sa == "M2"    L = 2;   electric = false
    elseif   sa == "E3"    L = 3;   electric = true
    elseif   sa == "M3"    L = 3;   electric = false
    elseif   sa == "E4"    L = 4;   electric = true
    elseif   sa == "M4"    L = 4;   electric = false
    else     error("stop a")
    end

    EmMultipole(L, electric)
end


# `Base.show(io::IO, mp::EmMultipole)`  ... prepares a proper printout of the variable mp::EmMultipole.
function Base.show(io::IO, mp::EmMultipole) 
    if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
    print(io, string(mp) )
end


# `Base.string(mp::EmMultipole)`  ... provides a proper printout of the variable mp::EmMultipole.
function Base.string(mp::EmMultipole)
    mp.L < 1    &&    error("Improper multipolarity L = $(mp.L)") 
    if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
    return( sa )
end


"""
`Basics.multipole_p(mp::EmMultipole)`  ... returns the kind integer p=0 (magnetic) or p=1 (electric) of the given  mp::EmMultipole.
"""
function multipole_p(mp::EmMultipole)
    if  mp.electric   return(1)    else    return(0)   end
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  EmProperty`  ... defines a type to maintain two gauge forms of a computed result that depends on the radiation field.

    + Coulomb         ::Float64    ... Value for the Coulomb gauge of the radiation field.
    + Babushkin       ::Float64    ... Value for the Coulomb gauge of the radiation field.
"""
struct  EmProperty 
    Coulomb           ::Float64
    Babushkin         ::Float64
end

export EmProperty


"""
`Basics.EmProperty(wa::Float64)`  ... constructor for an `constant` instance of EmProperty.
"""
function EmProperty(wa::Float64)
    EmProperty(wa, wa)
end


Base.:+(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb + b.Coulomb, a.Babushkin + b.Babushkin)
Base.:+(a::EmProperty, b) = EmProperty(a.Coulomb + b, a.Babushkin + b)
Base.:+(a, b::EmProperty) = b + a
Base.:-(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb - b.Coulomb, a.Babushkin - b.Babushkin)
Base.:*(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb * b.Coulomb, a.Babushkin * b.Babushkin)
Base.:*(a, b::EmProperty) = EmProperty(a * b.Coulomb, a * b.Babushkin)
Base.:/(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb / b.Coulomb, a.Babushkin / b.Babushkin)
Base.:/(a, b::EmProperty) = EmProperty(a / b.Coulomb, a / b.Babushkin)


# `Base.show(io::IO, property::EmProperty)`  ... prepares a proper printout of the variable property::EmProperty.
function Base.show(io::IO, property::EmProperty) 
    sa = Base.string(property);                print(io, sa)
end


# `Base.string(property::EmProperty)`  ... provides a String notation for the variable property::EmProperty.
function Base.string(property::EmProperty)
    sa = "$(property.Coulomb) [Coulomb],  $(property.Babushkin) [Babushkin]"
    return( sa )
end


#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.IsotopicFraction`  
    ... defines a type to deal with different -- individual or averaged -- isotopic fractions in a computation.

    + Z        ::Float64    ... nuclear charge
    + A        ::Float64    ... (mean) mass number of the isotope/isotopic mixture.
    + x        ::Float64    ... fraction 0 <= x <=1 of the given isotope/isotopic mixture.
"""
struct  IsotopicFraction 
    Z          ::Float64
    A          ::Float64
    x          ::Float64 
end

export IsotopicFraction


"""
`Basics.IsotopicFraction()`  ... constructor for an `empty` instance of IsotopicFraction.
"""
function IsotopicFraction()
    IsotopicFraction(0., 0., 0.)
end


# `Base.show(io::IO, frac::IsotopicFraction)`  ... prepares a proper printout of the variable frac::IsotopicFraction.
function Base.show(io::IO, frac::IsotopicFraction) 
    sa = Base.string(frac);                print(io, sa)
end


# `Base.string(frac::IsotopicFraction)`  ... provides a String notation for the variable frac::IsotopicFraction.
function Base.string(frac::IsotopicFraction)
    sa = "Isotopic fraction (Z=$(frac.Z), A=$(frac.A), x=$(frac.x))"
    return( sa )
end

#################################################################################################################################
#################################################################################################################################


"""
`@enum   Parity`  ... defines a enumeration for the allowed values of parity.

    + plus, minus               ... with obvious meaning
"""
@enum   Parity    plus   minus
export  Parity,   plus,  minus


"""
`Basics.Parity(sa::String)`  ... constructor for a given String.
"""
function Parity(sa::String)
    if       sa == "+"    wa = plus
    elseif   sa == "-"    wa = minus
    else     error("stop a")
    end

    return( wa )
end


"""
`Basics.invertParity(p::Parity)`  ... inverts the given parity plus <--> minus.
"""
function invertParity(p::Parity)
    if       p == plus    return( minus )
    elseif   p == minus   return( plus )
    else     error("stop a")
    end
end


# `Base.string(p::Parity)`  ... provides a proper printout of the variable p::Parity.
function Base.string(p::Parity) 
    if      p == plus     return( "+" )
    elseif  p == minus    return( "-" )  
    else    error("stop a")
    end
end


function  Base.:*(p1::Parity, p2::Parity)
    if      p1 == plus   &&  p2 == plus    return( plus ) 
    elseif  p1 == minus  &&  p2 == minus   return( plus ) 
    else                                   return( minus ) 
    end
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.Shell`  ... defines a enumeration for the allowed values of a non-relativistic shell.  

    + n     ::Int64  ... principal quantum number 
    + l     ::Int64  ... orbital angular quantum number
"""
struct  Shell
    n       ::Int64 
    l       ::Int64 
end

export  Shell


"""
`Basics.Shell(sa::Union{String,SubString{String}})`  ... constructor for a given String, such as 1s, 2s, 2p, ... .
"""
function Shell(sa::Union{String,SubString{String}}) 
    wa = strip( Base.convert(String, sa) )
    n  = parse(Int64, wa[1:end-1]);    l = shellNotation( string(wa[end]) )
    Shell(n, l)   
end


# `Base.show(io::IO, sh::Shell)`  ... prepares a proper printout of the variable sh::Shell.
function Base.show(io::IO, sh::Shell) 
    print(io, string(sh) )
end


# `Base.string(sh::Shell)`  ... provides a proper printout of the variable sh::Shell.
function Base.string(sh::Shell) 
    return( string(sh.n) * shellNotation(sh.l) )
end


"""
`Basics.shellNotation(l::Int64)`  ... returns the corresponding spectroscopy letter for a given orbital angular quantum number l.
"""
function shellNotation(l::Int64)
    !(0 <= l <= 21)      &&   return (string(l))
    wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
    return( wa[l+1] )  
    #
    ##x !(0 <= l <= 25)      &&   error("Orbital QN 0 <= l <= 25; l = $l")
    ##x wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
    ##x return( wa[l+1] )   
end


"""
`Basics.shellNotation(sa::String)`  ... returns for a given spectroscopy letter the orbital angular quantum number l.
"""
function shellNotation(sa::String)
    wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
    wb = findall(in([sa]), wa)
    !(length(wb)==1)      &&   error("Unrecognized orbital letter = $sa")
    return( wb[1]-1)   
end


"""
`Basics.shellSplitIntoSubshells(sh::Shell)'  ... to convert a non-relativistic shell into a list of Subshell[..]
"""
function shellSplitIntoSubshells(sh::Shell)
    subshells = Subshell[]
    if     sh.l == 0   push!(subshells, Subshell(sh.n, -1) )
    else   j2 = 2sh.l - 1;    kappa = Int64(  (j2+1)/2 );   push!(subshells, Subshell(sh.n, kappa) )
        j2 = 2sh.l + 1;    kappa = Int64( -(j2+1)/2 );   push!(subshells, Subshell(sh.n, kappa) )
    end
    return( subshells )
end


"""
`Basics.shellSplitOccupation(sh::Shell, occ::Int64)'  
    ... to split the occupation of a shell into proper subshell occupations; an  Array{Dict{Basics.Subshell,Int64},1} is returned.
"""
function shellSplitOccupation(sh::Shell, occ::Int64)

    subshells = shellSplitIntoSubshells(sh)
    wa        = Dict{Subshell,Int64}[]

    for  k1 in 0:occ
        for  k2 in 0:occ
            if  k1 + k2 == occ
                if       length( subshells ) == 1  
                    !(0 <= occ <= Basics.subshell_2j( subshells[1] ) + 1)   &&    error("Wrong pair of subshell = $(subshells[1]) and occupation = $occ") 
                    wb = Dict( subshells[1] => occ )
                    push!(wa, wb)
                    return( wa )

                    elseif  length( subshells ) == 2
                        j2_1 = Basics.subshell_2j(subshells[1])
                    j2_2 = Basics.subshell_2j(subshells[2])
                    if  k1 <= j2_1 + 1   &&    k2 <= j2_2 + 1    
                    wb = Dict( subshells[1] => k1, subshells[2] => k2 )
                    # Delete contributions with zero occupation
                    ## if  k1 == 0    delete!(wb, subshells[1])   end
                    ## if  k2 == 0    delete!(wb, subshells[2])   end
                    push!(wa, wb)
                    end
                else    error("stop a") 
                end
            end
        end
    end

    return( wa )
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.Subshell`  ... defines a type for the allowed values of a relativistic subhell.  

    + n        ::Int64  ... principal quantum number 
    + kappa    ::Int64  ... relativistic angular quantum number
"""
struct  Subshell
    n          ::Int64 
    kappa      ::Int64 
end

export  Subshell, subshell_2j



"""
`Basics.Subshell(sa::String)`  ... constructor for a given String, such as 1s_1/2, 2s_1/2, 2p_3/2, ... .
"""
function Subshell(sa::String) 
    wa = strip( sa );   
    wb = findnext("_", wa, 1);    if  wb == nothing   error("No underscore in sa = $sa")           else   wb = wb[1]   end
    wc = findnext("/", wa, wb);   if  wc == nothing   error("No / in wa[] = $(wa[wb[1]+1:end])")   else   wc = wc[1]   end
    wd = string( wa[wb+1:wc-1] )
    n  = parse(Int64, wa[1:wb-2]);                  l  = shellNotation( string(wa[ wb-1 ]) )
    j2 = parse(Int64, wd );                         !(string(wa[end]) == "2")    &&  error("Unrecognized subshell string sa = $sa") 
    if  l - j2/2 > 0     kappa = Int64( (j2+1)/2 )    else    kappa = Int64( -(j2+1)/2 )   end
    Subshell(n, kappa)   
end


# `Base.show(io::IO, sh::Subshell)`  ... prepares a proper printout of the variable sh::Subshell.
function Base.show(io::IO, sh::Subshell) 
    print(io, string(sh) )
end


# `Base.string(sh::Subshell)`  ... provides a proper printout of the variable sh::Subshell.
function Base.string(sh::Subshell) 
    sa = string(sh.n) * shellNotation( Basics.subshell_l(sh) ) * "_" * string( Basics.subshell_j(sh) )
    return( sa )
end


"""
`Basics.subshellGrasp(sa::Stringl)`  ... returns a sh::Subshell if sa denotes a (relativistic) Subshell in Grasp-notation.
"""
function subshellGrasp(sa::String)
        if   sa == "1s"   return( Subshell("1s_1/2")  )   
    elseif   sa == "2s"   return( Subshell("2s_1/2")  )   elseif   sa == "2p-"   return( Subshell("2p_1/2") )
    elseif   sa == "2p"   return( Subshell("2p_3/2")  ) 
    elseif   sa == "3s"   return( Subshell("3s_1/2")  )   elseif   sa == "3p-"   return( Subshell("3p_1/2") )
    elseif   sa == "3p"   return( Subshell("3p_3/2")  )   elseif   sa == "3d-"   return( Subshell("3d_3/2") )
    elseif   sa == "3d"   return( Subshell("3d_5/2")  ) 
    elseif   sa == "4s"   return( Subshell("4s_1/2")  )   elseif   sa == "4p-"   return( Subshell("4p_1/2") )
    elseif   sa == "4p"   return( Subshell("4p_3/2")  )   elseif   sa == "4d-"   return( Subshell("4d_3/2") )
    elseif   sa == "4d"   return( Subshell("4d_5/2")  )   elseif   sa == "4f-"   return( Subshell("4f_5/2") )
    elseif   sa == "4f"   return( Subshell("4f_7/2")  ) 
    elseif   sa == "5s"   return( Subshell("5s_1/2")  )   elseif   sa == "5p-"   return( Subshell("5p_1/2") )
    elseif   sa == "5p"   return( Subshell("5p_3/2")  )   elseif   sa == "5d-"   return( Subshell("5d_3/2") )
    elseif   sa == "5d"   return( Subshell("5d_5/2")  )   elseif   sa == "5f-"   return( Subshell("5f_5/2") )
    elseif   sa == "5f"   return( Subshell("5f_7/2")  )   elseif   sa == "5g-"   return( Subshell("5g_7/2") )
    elseif   sa == "5g"   return( Subshell("5g_9/2")  ) 
    elseif   sa == "6s"   return( Subshell("6s_1/2")  )   elseif   sa == "6p-"   return( Subshell("6p_1/2") )
    elseif   sa == "6p"   return( Subshell("6p_3/2")  )   elseif   sa == "6d-"   return( Subshell("6d_3/2") )
    elseif   sa == "6d"   return( Subshell("6d_5/2")  )   elseif   sa == "6f-"   return( Subshell("6f_5/2") )
    elseif   sa == "6f"   return( Subshell("6f_7/2")  )   elseif   sa == "6g-"   return( Subshell("6g_7/2") )
    elseif   sa == "6g"   return( Subshell("6g_9/2")  ) 
    elseif   sa == "7s"   return( Subshell("7s_1/2")  )   elseif   sa == "7p-"   return( Subshell("7p_1/2") )
    elseif   sa == "7p"   return( Subshell("7p_3/2")  )   elseif   sa == "7d-"   return( Subshell("7d_3/2") )
    else     error("Grasp notation ... not yet implemented; sa = $sa")
    end
end


"""
`Basics.subshell_l(sh::Subshell)`  ... returns the orbital angular quantum number l for a given sh::Subshell; an Int64 is returned
"""
function subshell_l(sh::Subshell)
    if   sh.kappa < 0   l = abs(sh.kappa) -1   else   l = sh.kappa   end
    return( l )   
end


"""
`Basics.subshell_j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
"""
function subshell_j(sh::Subshell)
    return( AngularJ64(abs(sh.kappa) - 1//2) )   
end


"""
`Basics.subshell_2j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
"""
function subshell_2j(sh::Subshell)
    return( 2* abs(sh.kappa) - 1 )   
end


"""
`Basics.subshellsFromClosedShellConfiguration("[Ne]")`  
    ... to provide a list of (relativistic) subshells for the given closed-shell configuration.
"""
function subshellsFromClosedShellConfiguration(sa::String)
    if       sa == "[He]"    wa = [ Subshell("1s_1/2")]    
    elseif   sa == "[Ne]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]   
    elseif   sa == "[Mg]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), Subshell("3s_1/2")]   
    elseif   sa == "[Ar]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                        Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2")]   
    elseif   sa == "[Kr]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                        Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"),   
                                                        Subshell("3d_3/2"), Subshell("3d_5/2"), 
                                    Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2")]   
    elseif   sa == "[Xe]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                        Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"),   
                                                        Subshell("3d_3/2"), Subshell("3d_5/2"), 
                                    Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2"), 
                                                        Subshell("4d_3/2"), Subshell("4d_5/2"), 
                                    Subshell("5s_1/2"), Subshell("5p_1/2"), Subshell("5p_3/2")]   
    else
        error("Unsupported keystring = $sa ")
    end
    
    return( wa )
end

#################################################################################################################################
#################################################################################################################################


"""
`struct   Basics.SubshellStateR`  
    ... defines a struct for the relativistic antisymmetric subshell states in the seniority scheme; this struct can be accessed 
        only internally and, therefore, the only the standard constructor is supported.

    + subshell     ::Subshell        ... to refer to a particular subshell. 
    + occ          ::Int64           ... occupation of the subshell
    + nu           ::Int64           ... seniority of the subshell state, often called nu
    + Jsub2        ::Int64           ... 2*Jsub
"""
struct  SubshellStateR
    subshell       ::Subshell
    occ            ::Int64
    nu             ::Int64
    Jsub2          ::Int64
end     

export  SubshellStateR


# `Base.show(io::IO, state::SubshellStateR)`  ... prepares a proper printout of the variable state::SubshellStateR.
function Base.show(io::IO, state::SubshellStateR) 
    print(io, string(state) )
end


# `Base.string(state::SubshellStateR)`  ... provides a proper printout of the variable state::SubshellStateR.
function Base.string(state::SubshellStateR)
    sa = "SubshellState: [" * string(state.subshell) * "^$(state.occ), nu=$(state.nu), 2*J_sub=$(state.Jsub2)]"
end

#################################################################################################################################
#################################################################################################################################


"""
`@enum   Basics.UseGauge`  ... defines a enumeration for the (allowed) gauges to be selected in explicit computations
"""
@enum   UseGauge    UseCoulomb   UseBabushkin
export  UseGauge,   UseCoulomb,  UseBabushkin


"""
`Basics.UseGauge(sa::String)`  ... constructor for a given String.
"""
function UseGauge(sa::String)
    if       sa == "Coulomb"                   wa = UseCoulomb
    elseif   sa == "Babushkin"                 wa = UseBabushkin
    else     error("stop a")
    end

    UseGauge(wa)
end


# `Base.string(gauge::UseGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
function Base.string(gauge::UseGauge) 
    if      gauge == UseCoulomb     return( "Coulomb" )   
    elseif  gauge == UseBabushkin   return( "Babushkin" )   
    else    error("stop a")
    end
end


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.CorePolarization`  
    ... defines a type to describe the influence of core-polarization upon the potential or interaction operator.

    + doApply          ::Bool                ... True if core-polarization is to be applied, and false otherwise.
    + coreAlpha        ::Float64             ... Polarizibility of the ionic core.
    + coreRadius       ::Float64             ... Radius of the ionic core.
    + coreShells       ::Array{Shell,1}      ... Shells of the ionic core
"""
struct  CorePolarization
    doApply            ::Bool     
    coreAlpha          ::Float64 
    coreRadius         ::Float64
    coreShells         ::Array{Shell,1} 
end 

export  CorePolarization


"""
`Basics.CorePolarization()`  ... constructor for the default values of core-polarization contributions.
"""
function CorePolarization()
    CorePolarization(false, 0., 0., Shell[])
end


# `Base.show(io::IO, cp::CorePolarization)`  ... prepares a proper printout of the variable cp::CorePolarization.
function Base.show(io::IO, cp::CorePolarization) 
    sa = Base.string(cp);                print(io, sa)
end


# `Base.string(cp::CorePolarization)`  ... provides a String notation for the variable cp::CorePolarization.
function Base.string(cp::CorePolarization)
    if  cp.doApply  sa = "Apply "   else    sa = "Do not apply "     end
    sa = sa * "core-polarization with alpha_c = $(cp.coreAlpha) a.u., r_c = $(cp.coreRadius) a.u. " *
                "core shells = $(cp.coreShells)."
    return( sa )
end

