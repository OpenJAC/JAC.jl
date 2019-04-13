#=
Here, we define all basic types/struct that are specific to the JAC module; we shall use proper docstrings in order to allow easy access to the 
individual subcomponents of each type/struct definition
=#

export  AngularJ64, AngularM64, HalfInt, HalfInteger, SolidAngle, Parity, LevelSymmetry, Shell, Subshell,
        EmMultipole, E1, M1, E2, M2, E3, M3, E4, M4, EmProperty, ExpStokes, EmStokes, TensorComp


include("inc-halfintegers.jl")


abstract type  AngularJ  end

"""
`struct  AngularJ64 <: AngularJ`  ... defines a type for angular momenta j.

    + num  ::Int64              ... numerator
    + den  ::Int64              ... denominator, must be 1 or 2
"""
struct  AngularJ64 <: AngularJ
    num      ::Int64
    den      ::Int64
end


"""
`JAC.AngularJ64(m::Integer)`  ... constructor for a given integer (numerator).
"""
function AngularJ64(j::Integer)
    j < 0   &&   error("j must be positive.")
    AngularJ64(j, 1)
end


"""
`JAC.AngularJ64(rational::Rational{Int64})`  ... constructor for a given  Rational{Int64}.
"""
function AngularJ64(rational::Rational{Int64})
    !(rational.den in [1,2])   &&   error("Denominator must be 1 or 2.")
    rational.num < 0           &&   error("j must be positive.")
    AngularJ64(rational.num, rational.den)
end


"""
`Base.show(io::IO, j::AngularJ64)`  ... prepares a proper printout of the variable j::JAC.AngularJ64.
"""
function Base.show(io::IO, j::JAC.AngularJ64) 
    if      j.den == 1    print(io, j.num)
    elseif  j.den == 2    print(io, j.num, "/2")
    else    error("stop a")
    end
end


"""
`JAC.oplus(ja::AngularJ64, jb::AngularJ64)`  ... adds the angular momenta ja `oplus` jb and returns a list::Array{AngularJ64,1} of j-values
                                                 in the interval |ja - jb| <= j <= ja + jb.
"""
function  oplus(ja::AngularJ64, jb::AngularJ64)
    if  ja.den == 1   ja2 = 2ja.num   else   ja2 = ja.num   end
    if  jb.den == 1   jb2 = 2jb.num   else   jb2 = jb.num   end
    jList = AngularJ64[]
    for  j = abs(ja2 - jb2):2:ja2+jb2    push!(jList, AngularJ64(j//2) )    end
    return( jList )
end

oplus(ja, jb) = oplus(HalfInt(ja), HalfInt(jb))
oplus(ja::Union{HalfInt,Integer}, jb::Union{HalfInt,Integer}) = abs(ja-jb):ja+jb

const ⊕ = oplus


"""
`JAC.projections(ja::AngularJ64)`  ... returns all allowed projections of a given angular momenta ja as a list::Array{AngularM64,1} 
                                      of m-values, i.e. -ja, -ja+1, ..., j.
"""
function  projections(ja::AngularJ64)
    if  ja.den == 1   ja2 = 2ja.num   else   ja2 = ja.num   end
    mList = AngularM64[]
    for  m = -ja2:2:ja2    push!(mList, AngularM64(m//2) )    end
    return( mList )
end

projections(j) = projections(HalfInt(j))
projections(j::Union{HalfInt,Integer}) =
    j≥0 ? (-j:j) : throw(DomainError(j, "angular momentum j must be non-negative."))


abstract type  AngularM  end

"""
`struct  AngularM64 <: AngularM`  ... defines a type for angular momentum projections m.

    + num  ::Int64              ... numerator
    + den  ::Int64              ... denominator, must be 1 or 2
"""
struct  AngularM64 <: AngularM
    num      ::Int64
    den      ::Int64
end


"""
`JAC.AngularM64(m::Integer)`  ... constructor for a given integer (numerator).
"""
function AngularM64(m::Integer)
    AngularM64(m, 1)
end


"""
`JAC.AngularM64(m::Integer, j::AngularJ64)`  ... constructor for a given integer (numerator) that is consistent with a given j-value.
"""
function AngularM64(m::Integer, j::AngularJ64)
    !(j.den == 1)      &&  error("m must be integer for j = $(j).")
      j.num < abs(m)   &&  error("abs(m) must be <= j = $(j).")
    AngularM64(m, 1)
end


"""
`JAC.AngularM64(rational::Rational{Int64})`  ... constructor for a given  Rational{Int64}.
"""
function AngularM64(rational::Rational{Int64})
    !(rational.den in [1,2])  &&  error("Denominator must be 1 or 2.")
    AngularM64(rational.num, rational.den)
end


"""
`JAC.AngularM64(rational::Rational{Int64}, j::AngularJ64)`  ... constructor for a given  Rational{Int64} that is consistent with a given m-value.
"""
function AngularM64(rational::Rational{Int64}, j::AngularJ64)
    !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
    !(j.den == rational.den)      &&   error("j,m must be both integer or half-integer.")
      j.num < abs(rational.num)   &&   error("abs(m) must be <= j = $(j).")
    AngularM64(rational.num, rational.den)
end


"""
`JAC.AngularM64(j::AngularJ64)`  ... constructor for a given j::AngularJ64 to define m = j.
"""
function AngularM64(j::AngularJ64)
    AngularM64(j.num, j.den)
end


"""
`Base.show(io::IO, m::AngularM64)`  ... prepares a proper printout of the variable m::JAC.AngularM64.
"""
function Base.show(io::IO, m::JAC.AngularM64) 
    if      m.den == 1    print(io, m.num)
    elseif  m.den == 2    print(io, m.num, "/2")
    else    error("stop a")
    end
end


"""
`JAC.add(ma::AngularM64, mb::AngularM64)`  
    ... adds the projections of the angular momenta ma + mb and returns a mc::AngularM64.
"""
function  add(ma::AngularM64, mb::AngularM64)
    if  ma.den == 1   ma2 = 2ma.num   else   ma2 = ma.num   end
    if  mb.den == 1   mb2 = 2mb.num   else   mb2 = mb.num   end
    return( AngularM64( (ma2+mb2)//2 ) )
end



# Conversion between HalfInt and AngularJ64, AngularM64

twice(x::Union{AngularJ64,AngularM64}) = ifelse(isone(x.den), twice(x.num), x.num)

AngularJ64(x::HalfInt) = isinteger(x) ? AngularJ64(Integer(x)) : AngularJ64(twice(x), 2)
AngularM64(x::HalfInt) = isinteger(x) ? AngularM64(Integer(x)) : AngularM64(twice(x), 2)



"""
`struct  Eigen`  ... defines a simple struct to communicate eigenvalues and eigenvectors if different diagonalization procedures are used.

    + values   ::Array{Float64,1}            ... List of eigenvalues.
    + vectors  ::Array{Vector{Float64},1}    ... List of eigenvectors.
"""
struct  Eigen
    values     ::Array{Float64,1}
    vectors    ::Array{Vector{Float64},1}
end


"""
`JAC.isEqual(a::Eigen, b::Eigen, accuracy::Float64)`  ... determines whether a == b with the requested accuracy; a value::Bool is returned.
                                                          The test is made element-wise on abs(a.number - b.number) > accuracy.
"""
function isEqual(a::Eigen, b::Eigen, accuracy::Float64)
    # Check the dimensions
    if  accuracy <= 0.                          error("Requires accuracy > 0.")                                    end
    if  size(a.values) != size(b.values)   ||   size(a.vectors) != size(b.vectors)                return( false )  end
    # Check the eigenvalues and vectors
    for  i = 1:length(a.values)    
        if  abs(a.values[i] - b.values[i]) > accuracy                                             return( false )  end    
        for  j = 1:length(a.vectors[i])   if  abs(a.vectors[i][j] - b.vectors[i][j]) > accuracy   return( false )  end    end     
    end
    return( true )
end





"""
`struct  SolidAngle`         ... defines a type for a solid angle Omega = (theta, phi).

    + theta  ::Float64       ... polar angle theta with regard to the z-axis; 0 <= theta <= pi.
    + phi    ::Float64       ... azimuthal angle phi (with regard to the x-z plane in mathematical positive direction); 0 <= phi <= 2 pi.
"""
struct  SolidAngle
    theta    ::Float64
    phi      ::Float64
end


"""
`Base.show(io::IO, Omega::SolidAngle)`  ... prepares a proper printout of the variable Omega::SolidAngle.
"""
function Base.show(io::IO, Omega::SolidAngle)
    sa = "Omega == (" * string(Omega.theta) * ", " * string(Omega.phi) * ")"
    print(io, sa)
end


"""
`@enum   Parity`  ... defines a enumeration for the allowed values of parity.

    + plus, minus               ... with obvious meaning
"""
@enum   Parity    plus   minus


"""
`JAC.Parity(sa::String)`  ... constructor for a given String.
"""
function Parity(sa::String)
    if       sa == "+"    wa = plus
    elseif   sa == "-"    wa = minus
    else     error("stop a")
    end

    return( wa )
end


"""
`JAC.invertParity(p::Parity)`  ... inverts the given parity plus <--> minus.
"""
function invertParity(p::Parity)
    if       p == plus    return( minus )
    elseif   p == minus   return( plus )
    else     error("stop a")
    end
end


"""
`Base.show(io::IO, p::Parity)`  ... prepares a proper printout of the variable p::Parity.
"""
function Base.show(io::IO, p::Parity) 
    print(io, string(p) )
end


"""
`Base.string(p::Parity)`  ... provides a proper printout of the variable p::Parity.
"""
function Base.string(p::Parity) 
    if      p == plus     return( "+" )
    elseif  p == minus    return( "-" )  
    else    error("stop a")
    end
end


"""
`struct   LevelSymmetry`  ... defines a struct for defining the overall J^P symmetry of a level.

    + J          ::AngularJ64  ... total angular momentum of a level
    + parity     ::Parity      ... total parity of the level
"""
struct  LevelSymmetry
    J            ::AngularJ64
    parity       ::Parity  
end


"""
`JAC.LevelSymmetry(rational::Rational{Int64}, parity::Parity)`  ... constructor for a given (Rational,Parity).
"""
function  LevelSymmetry(rational::Rational{Int64}, parity::Parity)
    !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
    LevelSymmetry( AngularJ64(rational.num, rational.den), parity )    
end


"""
`JAC.LevelSymmetry(i::Int64, parity::Parity)`  ... constructor for a given (Int64,Parity).
"""
function  LevelSymmetry(i::Int64, parity::Parity)
    LevelSymmetry( AngularJ64(i), parity )    
end


"""
`JAC.LevelSymmetry(rational::Rational{Int64}, sa::String)`  ... constructor for a given (Rational,String).
"""
function  LevelSymmetry(rational::Rational{Int64},sa::String)
    !(rational.den in [1,2])      &&   error("Denominator must be 1 or 2.")
    LevelSymmetry( AngularJ64(rational.num, rational.den), Parity(sa) )    
end


"""
`JAC.LevelSymmetry(i::Int64, sa::String)`  ... constructor for a given (Int64,String).
"""
function  LevelSymmetry(i::Int64,sa::String)
    LevelSymmetry( AngularJ64(i), Parity(sa) )    
end


"""
`Base.show(io::IO, sym::LevelSymmetry)`  ... prepares a proper printout of the variable sym::LevelSymmetry.
"""
function Base.show(io::IO, sym::LevelSymmetry) 
    print(io, string(sym) )
end


"""
`Base.string(sym::LevelSymmetry)`  ... provides a proper printout of the variable sym::LevelSymmetry.
"""
function Base.string(sym::LevelSymmetry) 
    if      sym.parity == plus     return( "$(sym.J) +" )
    elseif  sym.parity == minus    return( "$(sym.J) -" )  
    else    error("stop a")
    end
end


"""
`struct  Shell`  ... defines a enumeration for the allowed values of a non-relativistic shell.  

    + n     ::Int64  ... principal quantum number 
    + l     ::Int64  ... orbital angular quantum number
"""
struct  Shell
    n       ::Int64 
    l       ::Int64 
end


"""
`JAC.Shell(sa::Union{String,SubString{String}})`  ... constructor for a given String, such as 1s, 2s, 2p, ... .
"""
function Shell(sa::Union{String,SubString{String}}) 
    wa = strip( Base.convert(String, sa) )
    n  = parse(Int64, wa[1:end-1]);    l = shellNotation( string(wa[end]) )
    Shell(n, l)   
end


"""
`Base.show(io::IO, sh::Shell)`  ... prepares a proper printout of the variable sh::Shell.
"""
function Base.show(io::IO, sh::Shell) 
    print(io, string(sh) )
end


"""
`Base.string(sh::Shell)`  ... provides a proper printout of the variable sh::Shell.
"""
function Base.string(sh::Shell) 
    return( string(sh.n) * shellNotation(sh.l) )
end


"""
`JAC.shellNotation(l::Int64)`  ... returns the corresponding spectroscopy letter for a given orbital angular quantum number l.
"""
function shellNotation(l::Int64)
    !(0 <= l <= 25)      &&   error("Orbital QN 0 <= l <= 25; l = $l")
    wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
    return( wa[l+1] )   
end


"""
`JAC.shellNotation(sa::String)`  ... returns for a given spectroscopy letter the orbital angular quantum number l.
"""
function shellNotation(sa::String)
    wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
    wb = findall(in([sa]), wa)
    !(length(wb)==1)      &&   error("Unrecognized orbital letter = $sa")
    return( wb[1]-1)   
end


"""
`JAC.shellSplitIntoSubshells(sh::Shell)'  ... to convert a non-relativistic shell into a list of Subshell[..]
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
`JAC.shellSplitOccupation(sh::Shell, occ::Int64)'  ... to split the occupation of a shell into proper subshell occupations; 
                                                       an  Array{Dict{JAC.Subshell,Int64},1} is returned.
"""
function shellSplitOccupation(sh::Shell, occ::Int64)

    subshells = shellSplitIntoSubshells(sh)
    wa        = Dict{JAC.Subshell,Int64}[]

    for  k1 in 0:occ
        for  k2 in 0:occ
            if  k1 + k2 == occ
                if       length( subshells ) == 1  
                    !(0 <= occ <= subshell_2j( subshells[1] ) + 1)   &&    error("Wrong pair of subshell = $(subshells[1]) and occupation = $occ") 
                    wb = Dict( subshells[1] => occ )
                    push!(wa, wb)
                    return( wa )

                elseif  length( subshells ) == 2
                    j2_1 = subshell_2j(subshells[1])
                    j2_2 = subshell_2j(subshells[2])
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


"""
`struct  Subshell`  ... defines a type for the allowed values of a relativistic subhell.  

    + n        ::Int64  ... principal quantum number 
    + kappa    ::Int64  ... relativistic angular quantum number
"""
struct  Subshell
    n          ::Int64 
    kappa      ::Int64 
end


"""
`JAC.Subshell(sa::String)`  ... constructor for a given String, such as 1s_1/2, 2s_1/2, 2p_3/2, ... .
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


"""
`Base.show(io::IO, sh::Subshell)`  ... prepares a proper printout of the variable sh::Subshell.
"""
function Base.show(io::IO, sh::Subshell) 
    print(io, string(sh) )
end


"""
`Base.string(sh::Subshell)`  ... provides a proper printout of the variable sh::Subshell.
"""
function Base.string(sh::Subshell) 
    sa = string(sh.n) * shellNotation( subshell_l(sh) ) * "_" * string( subshell_j(sh) )
    return( sa )
end


"""
`JAC.subshellGrasp(sa::Stringl)`  ... returns a sh::Subshell if sa denotes a (relativistic) Subshell in Grasp-notation.
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
`JAC.subshell_l(sh::Subshell)`  ... returns the orbital angular quantum number l for a given sh::Subshell; an Int64 is returned
"""
function subshell_l(sh::Subshell)
    if   sh.kappa < 0   l = abs(sh.kappa) -1   else   l = sh.kappa   end
    return( l )   
end


"""
`JAC.subshell_j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
"""
function subshell_j(sh::Subshell)
    return( AngularJ64(abs(sh.kappa) - 1//2) )   
end


"""
`JAC.subshell_2j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
"""
function subshell_2j(sh::Subshell)
    return( 2* abs(sh.kappa) - 1 )   
end


"""
`JAC.subshellsFromClosedShellConfiguration("[Ne]")`  ... to provide a list of (relativistic) subshells for the given closed-shell configuration.
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


"""
`struct   SubshellStateR`  ... defines a struct for the relativistic antisymmetric subshell states in the seniority scheme; this struct 
                               can be accessed only internally and, therefore, the only the standard constructor is supported.

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


"""
`Base.show(io::IO, state::SubshellStateR)`  ... prepares a proper printout of the variable state::SubshellStateR.
"""
function Base.show(io::IO, state::SubshellStateR) 
    print(io, string(state) )
end


"""
`Base.string(state::SubshellStateR)`  ... provides a proper printout of the variable state::SubshellStateR.
"""
function Base.string(state::SubshellStateR)
    sa = "SubshellState: [" * string(state.subshell) * "^$(state.occ), nu=$(state.nu), 2*J_sub=$(state.Jsub2)]"
end


"""
`JAC.subshellStateString(subshell::String, occ::Int64, seniority::Int64, Jsub::AngularJ64, X::AngularJ64)`  ... to provide a 
                         string of a given subshell state in the form '[2p_1/2^occ]_(seniority, J_sub), X=Xo'
"""
function subshellStateString(subshell::String, occ::Int64, seniority::Int64, Jsub::AngularJ64, X::AngularJ64)
    sa = "[" * subshell * "^$occ]_($seniority, " * string(Jsub) * ") X=" * string(X)
    return( sa )
end


"""
`@enum   AtomicProcess`  ... defines a enumeration for the atomic processes that are supported in the JAC program in the computation 
                             of transition arrays or whole excitation/decay cascades.

    + Auger         ... Auger transitions, i.e. single autoionization or the emission of a single free electron into the continuum.
    + AugerInPlasma ... Auger transitions but calculated for a specified plasma model.
    + Compton       ... Rayleigh-Compton scattering cross sections.
    + Coulex        ... Coulomb-excitation of target or projeticle electrons by fast, heavy ions.
    + Coulion       ... Coulomb-ionization of target or projeticle electrons by fast, heavy ions.
    + Dierec        ... di-electronic recombination, i.e. the dielectronic capture of a free electron and the subsequent emission of a photon.
    + MultiPhotonDE ... multi-photon excitation and decay rates, including 2-photon, etc. processes.
    + MultiPI       ... multi-photon (single-electron) ionization.
    + MultiPDI      ... multi-photon (single-electron) double ionization.
    + Photo         ... Photoionization processes, i.e. the emission of a single free electron into the continuum due to an external light field.
    + PhotoExc      ... Photoexcitation rates.
    + PhotoExcFluor ... photoexcitation fluorescence rates and cross sections.
    + PhotoExcAuto  ... photoexcitation autoionization cross sections and collision strengths.
    + PhotoInPlasma ... Photoionization processes but calculated for a specified plasma model.
    + PhotoIonFluor ... photoionization fluorescence rates and cross sections.
    + PhotoIonAuto  ... photoionization autoionization cross sections and collision strengths.
    + Radiative    ... Radiative (multipole) transitions between bound-state levels of the same charge state.
    + Rec           ... radiative electron capture, i.e. the capture of a free electron with the simultaneous emission of a photon.
    + Eimex         ... electron-impact excitation cross sections and collision strengths.
    + RAuger        ... Radiative Auger rates.
"""
@enum   AtomicProcess  NoProcess  Auger  AugerInPlasma  Compton  Coulex  Dierec  Eimex  ImpactExcAuto  InternalConv  MultiPhotonDE  MultiPI  MultiPDI=
                   20  Photo  PhotoExc  PhotoExcAuto  PhotoExcFluor  PhotoInPlasma  PhotoIonAuto  PhotoIonFluor   Radiative  RAuger=
                   40  Rec  PairA1P  Coulion
export  AtomicProcess, NoProcess, Auger, AugerInPlasma, Compton, Coulex, Dierec, Eimex, ImpactExcAuto, InternalConv, MultiPhotonDE, MultiPI, MultiPDI,
                       Photo, PhotoExc, PhotoExcAuto, PhotoExcFluor, PhotoInPlasma, PhotoIonAuto, PhotoIonFluor,  Radiative, RAuger, 
                       Rec, PairA1P, Coulion


"""
`JAC.AtomicProcess(sa::String)`  ... constructor for a given String.
"""
function AtomicProcess(sa::String)
    if       sa in [ "none", "NoProcess"]                       wa = NoProcess
    elseif   sa in [ "Auger"]                                   wa = Auger
    elseif   sa in [ "Auger in plasma"]                         wa = AugerInPlasma
    elseif   sa in [ "Compton", "Rayleigh"]                     wa = Compton
    elseif   sa in [ "Coulex", "Coulomb excitation"]            wa = Coulex
    elseif   sa in [ "Coulion", "Coulomb ionization"]           wa = Coulion
    elseif   sa in [ "dierec", "Dierec"]                        wa = Dierec
    elseif   sa in [ "eimex", "Eimex"]                          wa = Eimex
    elseif   sa in [ "internal conversion"]                     wa = InternalConv
    elseif   sa in [ "multi-DE"]                                wa = MultiPhotonDE
    elseif   sa in [ "multi-PI"]                                wa = MultiPI
    elseif   sa in [ "multi-PDI"]                               wa = MultiPDI
    elseif   sa in [ "photo", "Photo"]                          wa = Photo
    elseif   sa in [ "photo in plasma"]                         wa = PhotoInPlasma
    elseif   sa in [ "PhotoExc", "photoexcitation"]             wa = PhotoExc
    elseif   sa in [ "photo-EA"]                                wa = PhotoExcAuto
    elseif   sa in [ "photo-IA"]                                wa = PhotoIonAuto
    elseif   sa in [ "photo-IF"]                                wa = PhotoIonFluor
    elseif   sa in [ "photo-EF"]                                wa = PhotoExcFluor
    elseif   sa in [ "radiative", "Radiative"]                  wa = Radiative
    elseif   sa in [ "radiativeAuger", "RAuger"]                wa = RAuger
    elseif   sa in [ "rec", "Rec"]                              wa = Rec
    elseif   sa in [ "pair-annihilation-1-photon", "PA1P"]      wa = PairA1P
    elseif   sa in [ "impact-EA"]                               wa = ImpactExcAuto
    elseif   sa in [ "rec", "Rec"]                              wa = Rec
    else     error("stop a")
    end

    AtomicProcess(wa)
end


"""
`Base.show(io::IO, process::AtomicProcess)`  ... prepares a proper printout of the variable process::AtomicProcess.
"""
function Base.show(io::IO, process::AtomicProcess) 
    print(io, string(process) )
end


"""
`Base.string(process::AtomicProcess)`  ... provides a proper printout of the variable process::AtomicProcess.
"""
function Base.string(process::AtomicProcess) 
    if      process == NoProcess       return( "no process" )
    elseif  process == Auger           return( "Auger" )  
    elseif  process == AugerInPlasma   return( "Auger in plasma" )  
    elseif  process == Compton         return( "Rayleigh-Compton" )  
    elseif  process == Dierec          return( "Dielectronic recombination" )  
    elseif  process == Eimex           return( "Eimex" )  
    elseif  process == MultiPhotonDE   return( "multi-photon excitation & decay" )  
    elseif  process == Photo           return( "Photo" )  
    elseif  process == PhotoExcFluor   return( "Photo-Excitation-Fluoresence" )  
    elseif  process == PhotoExcAuto    return( "Photo-Excitation-Autoionization" )  
    elseif  process == PhotoInPlasma   return( "Photo in plasma" )  
    elseif  process == PhotoIonFluor   return( "Photo-Ionization-Fluoresence" )  
    elseif  process == PhotoIonAuto    return( "Photo-Ionization-Autoionization" )  
    elseif  process == PhotoExc        return( "Photo-Excitation" )  
    elseif  process == Radiative      return( "Radiative" )  
    elseif  process == RAuger          return( "Radiative Auger" )  
    elseif  process == Rec             return( "Rec" )  
    else    error("stop a")
    end
end



"""
`struct  EmMultipole`  ... defines a struct for the em multipoles, E1, M1, E2, ...

    + L           ::Int64    ... multipolarity
    + electric    ::Bool     ... True if the multipole is an electric one, and false otherwise.
"""  
struct  EmMultipole
    L             ::Int64
    electric      ::Bool   
end


const  E1 = EmMultipole(1, true);     const  M1 = EmMultipole(1, false)
const  E2 = EmMultipole(2, true);     const  M2 = EmMultipole(2, false)
const  E3 = EmMultipole(3, true);     const  M3 = EmMultipole(3, false)
const  E4 = EmMultipole(4, true);     const  M4 = EmMultipole(4, false)


"""
`JAC.EmMultipole(sa::String)`  ... constructor for a given String.
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


"""
`Base.show(io::IO, mp::EmMultipole)`  ... prepares a proper printout of the variable mp::EmMultipole.
"""
function Base.show(io::IO, mp::EmMultipole) 
    if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
    print(io, string(mp) )
end


"""
`Base.string(mp::EmMultipole)`  ... provides a proper printout of the variable mp::EmMultipole.
"""
function Base.string(mp::EmMultipole)
    mp.L < 1    &&    error("Improper multipolarity L = $(mp.L)") 
    if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
    return( sa )
end


"""
`JAC.multipole_p(mp::EmMultipole)`  ... returns the kind integer p=0 (magnetic) or p=1 (electric) of the given  mp::EmMultipole.
"""
function multipole_p(mp::EmMultipole)
    if  mp.electric   return(1)    else    return(0)   end
end


"""
`@enum   EmGauge`  ... defines a enumeration for the (allowed) gauges to deal with the radiation field; this differs from the UseGauge
                       which can be selected for explicit computations.
                       
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
`JAC.EmGauge(sa::String)`  ... constructor for a given String.
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


"""
`Base.show(io::IO, gauge::EmGauge)`  ... prepares a proper printout of the variable gauge::EmGauge.
"""
function Base.show(io::IO, gauge::EmGauge) 
    print(io, string(gauge) )
end


"""
`Base.string(gauge::EmGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
"""
function Base.string(gauge::EmGauge) 
    if      gauge == NoGauge     return( "no gauge" )
    elseif  gauge == Coulomb     return( "Coulomb" )   
    elseif  gauge == Babushkin   return( "Babushkin" )   
    elseif  gauge == Magnetic    return( "Magnetic" )   
    else    error("stop a")
    end
end


"""
`@enum   UseGauge`  ... defines a enumeration for the (allowed) gauges to be selected in explicit computations
"""
@enum   UseGauge    UseCoulomb   UseBabushkin
export  UseGauge,   UseCoulomb,  UseBabushkin


"""
`JAC.UseGauge(sa::String)`  ... constructor for a given String.
"""
function UseGauge(sa::String)
    if       sa == "Coulomb"                   wa = UseCoulomb
    elseif   sa == "Babushkin"                 wa = UseBabushkin
    else     error("stop a")
    end

    UseGauge(wa)
end


"""
`Base.show(io::IO, gauge::UseGauge)`  ... prepares a proper printout of the variable gauge::UseGauge.
"""
function Base.show(io::IO, gauge::UseGauge) 
    print(io, string(gauge) )
end


"""
`Base.string(gauge::UseGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
"""
function Base.string(gauge::UseGauge) 
    if      gauge == UseCoulomb     return( "Coulomb" )   
    elseif  gauge == UseBabushkin   return( "Babushkin" )   
    else    error("stop a")
    end
end



"""
`struct  EmProperty`  ... defines a type to maintain two gauge forms of a computed result that depends on the radiation field.

    + Coulomb         ::Float64    ... Value for the Coulomb gauge of the radiation field.
    + Babushkin       ::Float64    ... Value for the Coulomb gauge of the radiation field.
"""
struct  EmProperty 
    Coulomb           ::Float64
    Babushkin         ::Float64
end 

Base.:+(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb + b.Coulomb, a.Babushkin + b.Babushkin)
Base.:+(a::EmProperty, b) = EmProperty(a.Coulomb + b, a.Babushkin + b)
Base.:+(a, b::EmProperty) = b + a


"""
`Base.show(io::IO, property::EmProperty)`  ... prepares a proper printout of the variable property::EmProperty.
"""
function Base.show(io::IO, property::EmProperty) 
    sa = Base.string(property);                print(io, sa)
end


"""
`Base.string(property::EmProperty)`  ... provides a String notation for the variable property::EmProperty.
"""
function Base.string(property::EmProperty)
    sa = "$(property.Coulomb) [Coulomb],  $(property.Babushkin) [Babushkin]"
    return( sa )
end


"""
`struct  ExpStokes`  ... defines a type to maintain the (experimentally) given Stokes parameter for the polarization of incoming
                         photons or electrons.

    + P1              ::Float64    ... Stokes P1 parameter.
    + P2              ::Float64    ... Stokes P2 parameter.
    + P3              ::Float64    ... Stokes P3 parameter.
 """
struct  ExpStokes 
    P1                ::Float64
    P2                ::Float64
    P3                ::Float64
end 


"""
`JAC.ExpStokes()`  ... constructor for an (empty or unpolarized) Stokes vector.
"""
function ExpStokes()
    ExpStokes(0., 0., 0.)
end


"""
`Base.show(io::IO, ps::ExpStokes)`  ... prepares a proper printout of the variable ps::ExpStokes.
"""
function Base.show(io::IO, ps::ExpStokes) 
    sa = Base.string(ps);                print(io, sa)
end


"""
`Base.string(ps::ExpStokes)`  ... provides a String notation for the variable ps::ExpStokes.
"""
function Base.string(ps::ExpStokes)
    sa = "exp. Stokes parameters P1 = $(ps.P1),  P2 = $(ps.P2),  P3 = $(ps.P3)"
    return( sa )
end


"""
`struct  EmStokes`  ... defines a type to maintain the (computed) given Stokes parameter for the polarization of emitted
                        photons or electrons.

    + P1              ::EmProperty    ... Stokes P1 parameter in Coulomb and Babushkin gauge
    + P2              ::EmProperty    ... Stokes P2 parameter.
    + P3              ::EmProperty    ... Stokes P3 parameter.
 """
struct  EmStokes 
    P1                ::EmProperty
    P2                ::EmProperty
    P3                ::EmProperty
end 



"""
`Base.show(io::IO, ps::EmStokes)`  ... prepares a proper printout of the variable ps::EmStokes.
"""
function Base.show(io::IO, ps::EmStokes) 
    sa = Base.string(ps);                print(io, sa)
end



"""
`Base.string(ps::EmStokes)`  ... provides a String notation for the variable ps::EmStokes.
"""
function Base.string(ps::EmStokes)
    sa = "Stokes parameter P1 = $(ps.P1),  P2 = $(ps.P2),  P3 = $(ps.P3)  "
    return( sa )
end


"""
`struct  TensorComp`  ... defines a type for a component of the statistical tensor as associated with an atomic or ionic level.

    + k               ::Int64          ... rank of the tensor component
    + q               ::Int64          ... q-component
    + comp            ::ComplexF64     ... value of the tensor component
"""
struct  TensorComp 
    k                 ::Int64
    q                 ::Int64 
    comp              ::ComplexF64
end 



"""
`Base.show(io::IO, st::TensorComp)`  ... prepares a proper printout of the variable st::TensorComp.
"""
function Base.show(io::IO, st::TensorComp) 
    sa = Base.string(st);                print(io, sa)
end



"""
`Base.string(st::TensorComp)`  ... provides a String notation for the variable st::TensorComp.
"""
function Base.string(st::TensorComp)
    sa = "Statistical tensor component rho_($(st.k),$(st.q)) = $(st.comp)  "
    return( sa )
end


"""
`JAC.tensorComp(k::Int64, q::Int64, tcomps::Array{TensorComp,1};  withZeros::Bool = false)`  ... returns the value of the statistical 
     tensor component rho_kq if it is contained in tcomps.
"""
function tensorComp(k::Int64, q::Int64, tcomps::Array{TensorComp,1};  withZeros::Bool = false)
    for  tc in  tcomps
        if   k == tc.k    &&    q == tc.q    return( tc.comp )    end
    end
    
    if   withZeros    return( ComplexF64(0.) )    else    error("Statistical tensor component not found for k = $k and q = $q ")    end
end


"""
`@enum   AtomicLevelProperty`  ... defines a enumeration for the atomic level properties that are supported in the JAC program.

    + NoProperty      ... No level property defined.
    + Einstein        ... Einstein A, B coefficients and oscillator strength; although not a 'level property', `Einstein` is treated within
                          a single basis and does not support to include relaxation effects, etc. It can be used for a quick overwiew to
                          transition probabilities or in cascade computations, however,
    + HFS             ... Hyperfine A and B parameters.
    ##x + IJF_Expansion   ... Expansion of atomic states in a IJF-coupled basis.
    + Isotope         ... Isotope shift M and F parameters.
    + LandeJ          ... Lande g_J factors.
    + LandeF          ... Lande g_F factors.
    + Polarizibility  ... static and dynamic polarizibilities.
    + Plasma          ... CI computations including interactions from various plasma models.
    + Zeeman          ... Zeeman splitting of fine-structure levels.
"""
@enum   AtomicLevelProperty   NoProperty   AlphaX   EinsteinX   FormF   Green  HFS   Isotope   LandeJ   LandeF  Plasma=
                         20   Polarity   Yields  Zeeman
export  AtomicLevelProperty,  NoProperty,  AlphaX,  EinsteinX,  FormF,  Green, HFS,  Isotope,  LandeJ,  LandeF, Plasma,  
                              Polarity,  Yields,  Zeeman


"""
`JAC.AtomicLevelProperty(sa::Union{String,SubString{String}})`  ... constructor for a given String.
"""
function AtomicLevelProperty(sa::Union{String,SubString{String}})
    sa = strip( Base.convert(String, sa) )
    if       sa in ["none", "NoProperty"]                wa = NoProperty
    elseif   sa in ["alpha"]                             wa = AlphaX
    elseif   sa in ["Einstein"]                          wa = EinsteinX
    elseif   sa in ["form factor"]                       wa = FormF
    elseif   sa in ["Einstein"]                          wa = EinsteinX
    elseif   sa in ["Green", "Greens function"]          wa = Green
    ##x elseif   sa in ["IJF_expansion", "IJF"]              wa = IJF_expansion
    elseif   sa == ["isotope shift", "IS"]               wa = Isotope
    elseif   sa == ["LandeJ", "g_J"]                     wa = LandeJ
    elseif   sa == ["LandeF", "g_F"]                     wa = LandeF
    elseif   sa == ["Plasma", "plasma"]                  wa = Plasma
    elseif   sa == ["polarity", "polarizibility"]        wa = Polarity
    elseif   sa == ["yields", "fluorescence yield"]      wa = Yields
    elseif   sa == ["Zeeman"]                            wa = Zeeman
    else     error("stop a")
    end

    AtomicLevelProperty(wa)
end


"""
`Base.show(io::IO, property::AtomicLevelProperty)`  ... prepares a proper printout of the variable property::AtomicLevelProperty.
"""
function Base.show(io::IO, property::AtomicLevelProperty) 
    print(io, string(property) )
end


"""
`Base.string(property::AtomicLevelProperty)`  ... provides a proper printout of the variable property::AtomicLevelProperty.
"""
function Base.string(property::AtomicLevelProperty) 
    if      property == NoProperty    return( "no property" )
    elseif  property == AlphaX        return( "alpha variation" )  
    elseif  property == EinsteinX     return( "Einstein A, B coefficients" )  
    elseif  property == FormF         return( "form & scattering factors" )  
    elseif  property == Green         return( "approximate Greens function" )  
    elseif  property == HFS           return( "HFS A, B parameters" )  
    elseif  property == Isotope       return( "isotope shift M and F parameters" )  
    elseif  property == LandeJ        return( "Lande g_J factors" )  
    elseif  property == LandeF        return( "Lande g_F factors" )  
    elseif  property == Plasma        return( "CI computations with plasma-type interactions" )  
    elseif  property == Polarity      return( "multipole polarizibility" )  
    elseif  property == Yields        return( "fluorescence & Auger yields" )  
    elseif  property == Zeeman        return( "Zeeman splitting in an external magnetic field" )  
    else    error("stop a")
    end
end



"""
`@enum   RadialMesh`  ... defines a enumeration for dealing with the radial mesh points and integration scheme.

    + MeshGrasp       ... use the exponential or exponential-linearized grid from Grasp92.
    + MeshGL          ... use a Gauss-Legendre grid with grid points due to the orderGL.
"""
@enum   RadialMesh    MeshGrasp    MeshGL
export  RadialMesh,   MeshGrasp,   MeshGL



"""
`@enum   QedModel`  ... defines a enumeration for dealing with the radiative (QED) contribution to the Hamiltonian matrix.

    + QedSydney       ... to use the local QED potentials by Flambaum & Ginges PRA 72, 052115 (2005).
    + QedPetersburg   ... to use the local QED potential similar as suggested by Shabaev et al., PRA 88, 012513 (2013). 
"""
@enum   QedModel    QedSydney    QedPetersburg
export  QedModel,   QedSydney,   QedPetersburg



"""
`@enum   ContinuumNormalization`  ... defines a enumeration for dealing with the normalization of continuum orbitals.

    + PureSine       ... normalize with regard to an (asymtotic) pure sine funtion, sin(kr).
    + CoulombSine    ... normalize with regard to an (asymtotic) Coulombic-sine funtion, sin(kr + ...).
    + OngRussek      ... normalize by following Ong & Russek (1973).
"""
@enum   ContinuumNormalization    PureSine    CoulombSine    OngRussek
export  ContinuumNormalization,   PureSine,   CoulombSine,   OngRussek



"""
`@enum   ContinuumSolutions`  ... defines a enumeration for solving the continuum orbitals in a given potential.

    + ContBessel              ... generate a pure Bessel function for the large component together with kinetic balance.
    + ContSine                ... generate a pure Sine function for the large component together with kinetic balance.
    + NonrelativisticCoulomb  ... generate a non-relativistic Coulomb function for the large component together with kinetic balance.
    + AsymptoticCoulomb       ... generate a pure (asymptotic) Coulombic function for both components.
    + BsplineGalerkin         ... generate a continuum orbital with the Galerkin method.
"""
@enum   ContinuumSolutions    ContBessel    ContSine     AsymptoticCoulomb    NonrelativisticCoulomb    BsplineGalerkin
export  ContinuumSolutions,   ContBessel,   ContSine,    AsymptoticCoulomb,   NonrelativisticCoulomb,   BsplineGalerkin



"""
`@enum   Warnings`  ... defines a enumeration for dealing with warnings that are made during a run or REPL session.
                        Cf. JAC.warn().

    + AddWarning    ... add a Warning to a warningList.
    + PrintWarnings ... print all warnings into a jac-war.report file.
    + ResetWarnings ... reset (empty) the warningList, usually at the beginning of a new run.
"""
@enum   Warnings    AddWarning    PrintWarnings    ResetWarnings 
export  Warnings,   AddWarning,   PrintWarnings,   ResetWarnings 


"""
`@enum   Guint`  ... defines a enumeration for dealing with graphical user interfaces (GUI).

    + Gui          ... use a graphical (interactive) user interface.
    + GuiSettings  ... use a graphical (interactive) user interface to define some settings.
    + Cascade      ... use a graphical (interactive) user interface to define some settings.
"""
@enum   Guint    Gui    GuiSettings   GuiCascade
export  Guint,   Gui,   GuiSettings,  GuiCascade



"""
`JAC.yesno(question::String, sa::String)`  ... Returns true if the answer 'yes' or 'y' is found, and false otherwise; 
                                               sa = {"Y", "N"} determines how the zero-String "" is interpreted. The given question is repeated 
                                               until a proper answer is obtained.
"""
function yesno(question::String, sa::String)
    while  true
        print(question * "  ");    reply = strip( readline(STDIN) )
        if      sa == "Y"
            if      reply in ["", "Y", "Yes", "y", "yes"]   return( true )
            elseif  reply in [    "N", "No",  "n", "no" ]   return( false )
                    println("... answer not recognized ... redo:")
            end
        elseif  sa == "N"
            if      reply in [    "Y", "Yes", "y", "yes"]   return( true )
            elseif  reply in ["", "N", "No",  "n", "no" ]   return( false )
                    println("... answer not recognized ... redo:")
            end
        end
    end
end


