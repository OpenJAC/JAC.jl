
"""
`module  JAC.PeriodicTable`  
    ... a submodel of JAC that contains methods and data from the periodic table of elements. Although we here collect 
        data for various elements, no attempt will be made to set-up a full compilation. Instead, useful information
        and some semi-empirical data are compiled to simplify the use of the JAC package in the future.
"""
module PeriodicTable

    using  Printf

    
    """
    `struct  PeriodicTable.Isotope`  ... defines a type for collecting data for an individual isotope.

        + something     ::Int64                 ... Something
    """
    struct  Isotope 
    end


    """
    `PeriodicTable.Isotope()`  ... simple constructur for an empty instance of an isotope.
    """
    function Isotope(i::Int64)
        Isotope()
    end
    
   
    """
    `struct  PeriodicTable.Element`  ... defines a type for collecting data for each element.

        + nuclearCharge     ::Int64                 ... Nuclear charge
        + symbol            ::Symbol                ... Symbol of the element.
        + name              ::String                ... Name of the element.
        + mass              ::Float64               ... Mean mass number of the element.
        + configuration     ::String                ... Ground configuration of the element.
        + isotopes          ::Array{Isotope,1}        ... List of known isotopes (need still to be defined)
    """
    struct Element
        nuclearCharge       ::Int64 
        symbol              ::Symbol
        name                ::String
        mass                ::Float64
        configuration       ::String
        isotopes            ::Array{Isotope,1}
    end 


    """
    `PeriodicTable.Element(charge::Int64, symbol::Symbol)`  ... simple constructur for defining an element
    """
    function Element(charge::Int64, symbol::Symbol)
        Element(charge, symbol, "not defined", 999.99, Isotope[])
    end


    # `Base.show(io::IO, element::PeriodicTable.Element)`  ... prepares a proper printout of the variable element::PeriodicTable.Element.
    function Base.show(io::IO, element::PeriodicTable.Element) 
        println(io, "Element  $(element.name) ($(element.symbol), Z=$(element.nuclearCharge), mass=$(element.mass))  ")
        println(io, "configuration:            $(element.configuration)  ")
    end

    
    """
    `PeriodicTable.data::Dict{Int64,PeriodicTable.Element}`
        ... provides all available data for an element with given nuclear charge Z::Int64 by PeriodicTable.data[Z].
    """
    const data = Dict{Int64,PeriodicTable.Element}(
                1 => Element(  1, :H,  "hydrogen",     999.99,  "1s",                           Isotope[]),
                2 => Element(  2, :He, "helium",       999.99,  "1s^2",                         Isotope[]),
                3 => Element(  3, :Li, "lithum",       999.99,  "[He] 2s",                      Isotope[]),
                4 => Element(  4, :Be, "beryllium",    999.99,  "[He] 2s^2",                    Isotope[]),
                5 => Element(  5, :B,  "boron",        999.99,  "[He] 2s^2 2p",                 Isotope[]),
                6 => Element(  6, :C,  "carbon",       999.99,  "[He] 2s^2 2p^2",               Isotope[]),
                7 => Element(  7, :N,  "nitrogen",     999.99,  "[He] 2s^2 2p^3",               Isotope[]),
                8 => Element(  8, :O,  "oxygen",       999.99,  "[He] 2s^2 2p^4",               Isotope[]),
                9 => Element(  9, :F,  "fluorine",     999.99,  "[He] 2s^2 2p^5",               Isotope[]),
               10 => Element( 10, :Ne, "neon",         999.99,  "[Ne]",                         Isotope[]),
               11 => Element( 11, :Na, "sodium",       999.99,  "[Ne] 3s",                      Isotope[]),
               12 => Element( 12, :Mg, "magnesium",    999.99,  "[Ne] 3s^2",                    Isotope[]),
               13 => Element( 13, :Al, "aluminium",    999.99,  "[Ne] 3s^2 3p",                 Isotope[]),
               14 => Element( 14, :Si, "silicon",      999.99,  "[Ne] 3s^2 3p^2",               Isotope[]),
               15 => Element( 15, :P,  "phosphorus",   999.99,  "[Ne] 3s^2 3p^3",               Isotope[]),
               16 => Element( 16, :S,  "sulfur",       999.99,  "[Ne] 3s^2 3p^4",               Isotope[]),
               17 => Element( 17, :Cl, "chlorine",     999.99,  "[Ne] 3s^2 3p^5",               Isotope[]),
               18 => Element( 18, :Ar, "argon",        999.99,  "[Ar]",                         Isotope[]),
               19 => Element( 19, :K,  "potassium",    999.99,  "[Ar] 4s",                      Isotope[]),
               20 => Element( 20, :Ca, "calcium",      999.99,  "[Ar] 4s^2",                    Isotope[]),
               21 => Element( 21, :Sc, "scandium",     999.99,  "[Ar] 3d 4s^2",                 Isotope[]),
               22 => Element( 22, :Ti, "titanium",     999.99,  "[Ar] 3d^2 4s^2",               Isotope[]),
               23 => Element( 23, :V,  "vanadium",     999.99,  "[Ar] 3d^3 4s^2",               Isotope[]),
               24 => Element( 24, :Cr, "chromium",     999.99,  "[Ar] 3d^5 4s",                 Isotope[]),
               25 => Element( 25, :Mn, "manganese",    999.99,  "[Ar] 3d^5 4s^2",               Isotope[]),
               26 => Element( 26, :Fe, "iron",         999.99,  "[Ar] 3d^6 4s^2",               Isotope[]),
               27 => Element( 27, :Co, "cobalt",       999.99,  "[Ar] 3d^7 4s^2",               Isotope[]),
               28 => Element( 28, :Ni, "nickel",       999.99,  "[Ar] 3d^8 4s^2",               Isotope[]),
               29 => Element( 29, :Cu, "copper",       999.99,  "[Ar] 3d^10 4s",                Isotope[]),
               30 => Element( 30, :Zn, "zinc",         999.99,  "[Ar] 3d^10 4s^2",              Isotope[]),
               31 => Element( 31, :Ga, "gallium",      999.99,  "[Ar] 3d^10 4s^2 4p",           Isotope[]),
               32 => Element( 32, :Ge, "germanium",    999.99,  "[Ar] 3d^10 4s^2 4p^2",         Isotope[]),
               33 => Element( 33, :As, "arsenic",      999.99,  "[Ar] 3d^10 4s^2 4p^3",         Isotope[]),
               34 => Element( 34, :Se, "selenium",     999.99,  "[Ar] 3d^10 4s^2 4p^4",         Isotope[]),
               35 => Element( 35, :Br, "bromine",      999.99,  "[Ar] 3d^10 4s^2 4p^5",         Isotope[]),
               36 => Element( 36, :Kr, "krypton",      999.99,  "[Kr]",                         Isotope[]),
               37 => Element( 37, :Rb, "rubidium",     999.99,  "[Kr] 5s",                      Isotope[]),
               38 => Element( 38, :Sr, "strontium",    999.99,  "[Kr] 5s^2",                    Isotope[]),
               39 => Element( 39, :Y,  "yttrium",      999.99,  "[Kr] 4d 5s^2",                 Isotope[]),
               40 => Element( 40, :Zr, "zirconium",    999.99,  "[Kr] 4d^2 5s^2",               Isotope[]),
               41 => Element( 41, :Nb, "niobium",      999.99,  "[Kr] 4d^4 5s",                 Isotope[]),
               42 => Element( 42, :Mo, "molybdenum",   999.99,  "[Kr] 4d^5 5s",                 Isotope[]),
               43 => Element( 43, :Tc, "technetium",   999.99,  "[Kr] 4d^5 5s^2",               Isotope[]),
               44 => Element( 44, :Ru, "ruthenium",    999.99,  "[Kr] 4d^7 5s",                 Isotope[]),
               45 => Element( 45, :Rh, "rhodium",      999.99,  "[Kr] 4d^8 5s",                 Isotope[]),
               46 => Element( 46, :Pd, "palladium",    999.99,  "[Kr] 4d^10",                   Isotope[]),
               47 => Element( 47, :Ag, "silver",       999.99,  "[Kr] 4d^10 5s",                Isotope[]),
               48 => Element( 48, :Cd, "cadmium",      999.99,  "[Kr] 4d^10 5s^2",              Isotope[]),
               49 => Element( 49, :In, "indium",       999.99,  "[Kr] 4d^10 5s^2 5p",           Isotope[]),
               50 => Element( 50, :Sn, "tin",          999.99,  "[Kr] 4d^10 5s^2 5p^2",         Isotope[]), 
               51 => Element( 51, :Sb, "antimony",     999.99,  "[Kr] 4d^10 5s^2 5p^3",         Isotope[]),
               52 => Element( 52, :Te, "tellurium",    999.99,  "[Kr] 4d^10 5s^2 5p^4",         Isotope[]),
               53 => Element( 53, :I,  "iodine",       999.99,  "[Kr] 4d^10 5s^2 5p^5",         Isotope[]),
               54 => Element( 54, :Xe, "xenon",        999.99,  "[Xe]",                         Isotope[]),
               55 => Element( 55, :Cs, "caesium",      999.99,  "[Xe] 6s",                      Isotope[]),
               56 => Element( 56, :Ba, "barium",       999.99,  "[Xe] 6s^2",                    Isotope[]),
               57 => Element( 57, :La, "lanthanum",    999.99,  "[Xe] 5d 6s^2",                 Isotope[]),
               58 => Element( 58, :Ce, "cerium",       999.99,  "[Xe] 4f 5d 6s^2",              Isotope[]),
               59 => Element( 59, :Pr, "praseodymium", 999.99,  "[Xe] 4f^3 6s^2",               Isotope[]),
               60 => Element( 60, :Nd, "neodymium",    999.99,  "[Xe] 4f^4 6s^2",               Isotope[]),
               61 => Element( 61, :Pm, "promethium",   999.99,  "[Xe] 4f^5 6s^2",               Isotope[]),
               62 => Element( 62, :Sm, "samarium",     999.99,  "[Xe] 4f^6 6s^2",               Isotope[]),
               63 => Element( 63, :Eu, "europium",     999.99,  "[Xe] 4f^7 6s^2",               Isotope[]),
               64 => Element( 64, :Gd, "gadolinium",   999.99,  "[Xe] 4f^7 5d 6s^2",            Isotope[]),
               65 => Element( 65, :Tb, "terbium",      999.99,  "[Xe] 4f^9 6s^2",               Isotope[]),
               66 => Element( 66, :Dy, "dysprosium",   999.99,  "[Xe] 4f^10 6s^2",              Isotope[]),
               67 => Element( 67, :Ho, "holmium",      999.99,  "[Xe] 4f^11 6s^2",              Isotope[]),
               68 => Element( 68, :Er, "erbium",       999.99,  "[Xe] 4f^12 6s^2",              Isotope[]),
               69 => Element( 69, :Tm, "thulium",      999.99,  "[Xe] 4f^13 6s^2",              Isotope[]),
               70 => Element( 70, :Yb, "ytterbium",    999.99,  "[Xe] 4f^14 6s^2",              Isotope[]),
               71 => Element( 71, :Lu, "lutetium",     999.99,  "[Xe] 4f^14 5d 6s^2",           Isotope[]),
               72 => Element( 72, :Hf, "hafnium",      999.99,  "[Xe] 4f^14 5d^2 6s^2",         Isotope[]),
               73 => Element( 73, :Ta, "tantalum",     999.99,  "[Xe] 4f^14 5d^3 6s^2",         Isotope[]),
               74 => Element( 74, :W,  "tungsten",     999.99,  "[Xe] 4f^14 5d^4 6s^2",         Isotope[]),
               75 => Element( 75, :Re, "rhenium",      999.99,  "[Xe] 4f^14 5d^5 6s^2",         Isotope[]),
               76 => Element( 76, :Os, "osmium",       999.99,  "[Xe] 4f^14 5d^6 6s^2",         Isotope[]),
               77 => Element( 77, :Ir, "iridium",      999.99,  "[Xe] 4f^14 5d^7 6s^2",         Isotope[]),
               78 => Element( 78, :Pt, "platinum",     999.99,  "[Xe] 4f^14 5d^9 6s",           Isotope[]),
               79 => Element( 79, :Au, "gold",         999.99,  "[Xe] 4f^14 5d^10 6s",          Isotope[]),
               80 => Element( 80, :Hg, "mercury",      999.99,  "[Xe] 4f^14 5d^10 6s^2",        Isotope[]),
               81 => Element( 81, :Tl, "thallium",     999.99,  "[Xe] 4f^14 5d^10 6s^2 6p",     Isotope[]),
               82 => Element( 82, :Pb, "lead",         999.99,  "[Xe] 4f^14 5d^10 6s^2 6p^2",   Isotope[]),
               83 => Element( 83, :Bi, "bismuth",      999.99,  "[Xe] 4f^14 5d^10 6s^2 6p^3",   Isotope[]),
               84 => Element( 84, :Po, "polonium",     999.99,  "[Xe] 4f^14 5d^10 6s^2 6p^4",   Isotope[]),
               85 => Element( 85, :At, "astatine",     999.99,  "[Xe] 4f^14 5d^10 6s^2 6p^5",   Isotope[]),
               86 => Element( 86, :Rn, "radon",        999.99,  "[Rn]",                         Isotope[]),
               87 => Element( 87, :Fr, "francium",     999.99,  "[Rn] 7s",                      Isotope[]),
               88 => Element( 88, :Ra, "radium",       999.99,  "[Rn] 7s^2",                    Isotope[]),
               89 => Element( 89, :Ac, "actinium",     999.99,  "[Rn] 6d 7s^2",                 Isotope[]),
               90 => Element( 90, :Th, "thorium",      999.99,  "[Rn] 6d^2 7s^2",               Isotope[]),
               91 => Element( 91, :Pa, "protactinium", 999.99,  "[Rn] 5f^2 6d 7s^2",            Isotope[]),
               92 => Element( 92, :U,  "uranium",      999.99,  "[Rn] 5f^3 6d 7s^2",            Isotope[]),
               93 => Element( 93, :Np, "neptunium",    999.99,  "[Rn] 5f^4 6d 7s^2",            Isotope[]),
               94 => Element( 94, :Pu, "plutonium",    999.99,  "[Rn] 5f^6 7s^2",               Isotope[]),
               95 => Element( 95, :Am, "americium",    999.99,  "[Rn] 5f^7 7s^2",               Isotope[]),
               96 => Element( 96, :Cm, "curium",       999.99,  "[Rn] 5f^7 6d 7s^2",            Isotope[]),
               97 => Element( 97, :Bk, "berkelium",    999.99,  "[Rn] 5f^9 7s^2",               Isotope[]),
               98 => Element( 98, :Cf, "californium",  999.99,  "[Rn] 5f^10 7s^2",              Isotope[]),
               99 => Element( 99, :Es, "einsteinium",  999.99,  "[Rn] 5f^11 7s^2",              Isotope[]),
              100 => Element(100, :Fm, "fermium",      999.99,  "[Rn] 5f^12 7s^2",              Isotope[]),
              101 => Element(101, :Md, "mendelevium",  999.99,  "[Rn] 5f^13 7s^2",              Isotope[]),
              102 => Element(102, :No, "nobelium",     999.99,  "[Rn] 5f^14 7s^2",              Isotope[]),
              103 => Element(103, :Lr, "lawrencium",   999.99,  "[Rn] 5f^14 7s^2 7p",           Isotope[]),
              104 => Element(104, :Rf, "rutherfordium",999.99,  "[Rn] 5f^14 6d^2 7s^2",         Isotope[]),
              105 => Element(105, :Db, "dubnium",      999.99,  "[Rn] 5f^14 6d^3 7s^2",         Isotope[]),
              106 => Element(106, :Sg, "seaborgium",   999.99,  "[Rn] 5f^14 6d^4 7s^2",         Isotope[]),
              107 => Element(107, :Bh, "bohrium",      999.99,  "[Rn] 5f^14 6d^5 7s^2",         Isotope[]),
              108 => Element(108, :Hs, "hassium",      999.99,  "[Rn] 5f^14 6d^6 7s^2",         Isotope[]),
              109 => Element(109, :Mt, "meitnerium",   999.99,  "[Rn] 5f^14 6d^7 7s^2",         Isotope[]),
              110 => Element(110, :Ds, "darmstadtium", 999.99,  "[Rn] 5f^14 6d^8 7s^2",         Isotope[]),
              111 => Element(111, :Rg, "roentgenium",  999.99,  "[Rn] 5f^14 6d^9 7s^2",         Isotope[]),
              112 => Element(112, :Cn, "copernicium",  999.99,  "[Rn] 5f^14 6d^10 7s^2",        Isotope[]),
              113 => Element(113, :Nh, "nihonium",     999.99,  "[Rn] 5f^14 6d^8 7s^2 7p",      Isotope[]),
              114 => Element(114, :Fl, "flerovium",    999.99,  "[Rn] 5f^14 6d^2 7s^2 7p^2",    Isotope[]),
              115 => Element(115, :Mc, "moscovium",    999.99,  "[Rn] 5f^14 6d^3 7s^2 7p^3",    Isotope[]),
              116 => Element(116, :Lv, "livermorium",  999.99,  "[Rn] 5f^14 6d^4 7s^2 7p^4",    Isotope[]),
              117 => Element(117, :Ts, "tennessine",   999.99,  "[Rn] 5f^14 6d^5 7s^2 7p^5",    Isotope[]),
              118 => Element(118, :Og, "oganesson",    999.99,  "[Rn] 5f^14 6d^6 7s^2 7p^6",    Isotope[]) )


    """
    `PeriodicTable.bindingEnergies_Williams2000(Z::Int64)`  ... former  JAC.store_Williams2000(Z::Int64)
        ... to return the `stored values' of Williams et al. (2000), ... for the element with nuclear charge Z.
    """
    function bindingEnergies_Williams2000(Z::Int64)
        # This procedure 'stores' the binding energies of the inner-shell electrons 1s, ..., 3d for all elements from Z = 1, ..., 54;
        # all binding energies are given in eV
        # 
        #                        1s_1/2,   2s_1/2,   2p_1/2,   2p_3/2,   3s_1/2,   3p_1/2,   3p_3/2,   3d_3/2,   3d_5/2,   4s_1/2,   4p_1/2,   4p_3/2 
        if     Z ==  1   wa = [    13.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z ==  2   wa = [    24.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z ==  3   wa = [    54.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z ==  4   wa = [   111.5,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
        elseif Z ==  5   wa = [   188.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z ==  6   wa = [   284.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z ==  7   wa = [   409.9,     37.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z ==  8   wa = [   543.1,     41.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z ==  9   wa = [   696.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 10   wa = [   870.2,     48.5,     21.7,     21.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        #
        elseif Z == 11   wa = [  1070.8,     63.5,    30.65,    30.81,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
        elseif Z == 12   wa = [  1303.0,     88.7,    49.78,    49.50,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] 
        elseif Z == 13   wa = [  1559.6,    117.8,    72.95,    72.55,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 14   wa = [  1839.0,    149.7,    99.82,    99.42,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 15   wa = [  2145.5,    189.0,   136.0,    135.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 16   wa = [  2472.0,    230.9,   163.6,    162.5,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 17   wa = [  2822.4,    270.0,   202.0,    200.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 18   wa = [  3205.9,    326.3,   250.6,    248.4,      29.3,     15.9,     15.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 19   wa = [  3608.4,    378.6,   297.3,    294.6,      34.8,     18.3,     18.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 20   wa = [  4038.5,    438.4,   349.7,    346.2,      44.3,     25.4,     25.4,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        #
        elseif Z == 21   wa = [  4492.0,    498.0,   403.6,    398.7,      51.1,     28.3,     28.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 22   wa = [  4966.0,    560.9,   460.2,    453.8,      58.7,     32.6,     32.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 23   wa = [  5465.0,    626.7,   519.8,    512.1,      66.3,     37.2,     37.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 24   wa = [  5989.0,    696.0,   583.8,    574.1,      74.1,     42.2,     42.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 25   wa = [  6539.0,    769.1,   649.9,    638.7,      82.3,     47.2,     47.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 26   wa = [  7112.0,    844.6,   719.9,    706.8,      91.3,     52.7,     52.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 27   wa = [  7709.0,    925.1,   793.2,    778.1,     101.0,     58.9,     59.9,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 28   wa = [  8333.0,   1008.6,   870.0,    852.7,     110.8,     68.0,     66.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 29   wa = [  8979.0,   1096.7,   952.3,    932.7,     122.5,     77.3,     75.1,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 30   wa = [  9659.0,   1196.2,  1044.9,   1021.8,     139.8,     91.4,     88.6,     10.2,     10.1,     -1.0,     -1.0,     -1.0]
        #
        elseif Z == 31   wa = [ 10367.0,   1299.0,  1143.2,   1116.4,     159.5,    103.5,    100.0,     18.7,     18.7,     -1.0,     -1.0,     -1.0]
        elseif Z == 32   wa = [ 11103.0,   1414.6,  1248.1,   1217.0,     180.1,    124.9,    120.8,     29.8,     29.2,     -1.0,     -1.0,     -1.0]
        elseif Z == 33   wa = [ 11867.0,   1527.0,  1359.1,   1323.6,     204.7,    146.2,    141.2,     41.7,     41.7,     -1.0,     -1.0,     -1.0]
        elseif Z == 34   wa = [ 12658.0,   1652.0,  1474.3,   1433.9,     229.6,    166.5,    160.7,     55.5,     54.6,     -1.0,     -1.0,     -1.0]
        elseif Z == 35   wa = [ 13474.0,   1782.0,  1596.0,   1550.0,     257.0,    189.0,    182.0,     70.0,     69.0,     -1.0,     -1.0,     -1.0]
        elseif Z == 36   wa = [ 14326.0,   1921.0,  1730.9,   1678.4,     292.8,    222.2,    214.4,     95.0,     93.8,     27.5,     14.1,     14.1]
        #
        else  error("Data not available for Z = $Z")
        end
        return( wa )
    end

end # module
