
"""
`module  JAC.PeriodicTable`  
... a submodel of JAC that contains methods and data from the periodic table of elements. Although we here collect 
    data for various elements, no attempt will be made to set-up a full compilation. Instead, useful information
    and some semi-empirical data are compiled to simplify the use of the JAC package in the future.
"""
module PeriodicTable


using  Printf, ..Basics, ..ManyElectron


"""
`PeriodicTable.getAtomicNumber(symbol::Symbol)`  
    ... returns the atomic number (nuclear charge Z) for the given element symbol.
"""
function getAtomicNumber(symbol::Symbol)
    for  (k,v) in PeriodicTable.data
        if  v.symbol == symbol    return( v.nuclearCharge )    end
    end
    
    error("Undefined element symbol in the periodic table.")
end


"""
`PeriodicTable.getData()`  
    ... provides various data from the periodic table of elements. 

+ `("mass", Z::Int64)`  or  `("mass", symbol::Symbol)`   
    ... to get mean mass value of the element with nuclear charge Z or symbol.

+ `("density", Z::Int64)`  or  `("density", symbol::Symbol)`   
    ... to get the density [g/cm^3] of the element with nuclear charge Z or symbol.

+ `("1st IP", Z::Int64)`  or  `("1st IP", symbol::Symbol)`   
    ... to get the first ionization potential [eV] of the element with nuclear charge Z or symbol.

+ `("ground configuration", Z::Int64)`  or  `("ground configuration", symbol::Symbol)`   
    ... to get the ground configuration of the neutral atom with nuclear charge Z or symbol.

+ `("atomic radius", Z::Int64)`  or  `("atomic radius", symbol::Symbol)`   
    ... to get the atomic radius of the neutral atom with nuclear charge Z or symbol.

+ `("polarizibility", Z::Int64)`  or  `("polarizibility", symbol::Symbol)`   
    ... to get the polarizibility of the neutral atom with nuclear charge Z or symbol.

+ `("ground configuration", Z::Int64)`  or  `("ground configuration", symbol::Symbol)`   
    ... to get the ground configuration of the neutral atom with nuclear charge Z or symbol.

+ `("chemical hardness", Z::Int64)`  or  `("chemical hardness", symbol::Symbol)`   
    ... to get the chemical hardness of the neutral atom with nuclear charge Z or symbol.
"""
function getData(sa::String, zsymbol::Union{Int64,Symbol})
    if  typeof(zsymbol) == Symbol   Z = PeriodicTable.getAtomicNumber(zsymbol)    else    Z = zsymbol     end

    if        sa == "mass"              return ( PeriodicTable.data[Z].mass )   
    elseif    sa == "density"           println("\nDensities are given [g/cm^3]: ...")
                                        return ( PeriodicTable.data[Z].density )   
    elseif    sa == "1st IP"            println("\nIonization potentials are given [eV]: ...")
                                        return ( PeriodicTable.data[Z].firstIP )   
    elseif    sa == "ground configuration"        
                                        return ( Configuration(PeriodicTable.data[Z].configuration) )   
    elseif    sa == "atomic radius"     println("\nAtomic radius is given [Angstroem]: ...")
                                        return ( PeriodicTable.data[Z].atomicRadius )   
    elseif    sa == "polarizibility"    println("\nAtomic polarizibility is given [10^-24 c.c]: ...")
                                        return ( PeriodicTable.data[Z].polarizibility )   
    elseif    sa == "chemical hardness" println("\nChemical hardness is given [eV]: ...")
                                        return ( PeriodicTable.data[Z].chemHardness )   
    else      error("Unsupported keystring:: $sa")
    end
end


"""
`PeriodicTable.getIsotopeData()`  
    ... provides various data from the periodic table of elements. Not yet implemented.

+ `("mass", Z::Int64)`  or  `("mass", symbol::Symbol)`   
    ... to get mean mass value of the element with nuclear charge Z or symbol

"""
function getIsotopeData(sa::String, zsymbol::Union{Int64,Symbol})

    if        sa in ["alpha", "fine-structure constant alpha"]          return (FINE_STRUCTURE_CONSTANT)   
    else      error("Unsupported keystring:: $sa")
    end
end



########################################################################################################################
########################################################################################################################
########################################################################################################################



"""
`struct  PeriodicTable.Isotope`  ... defines a type for collecting data for an individual isotope.

    + massNumber         ::Float64           ... Mass number of the isotope.
    + isotopeMass        ::Float64           ... Mass of the isotope.
    + spin               ::AngularJ64        ... Nuclear spin of the isotope
    + mu                 ::Float64           ... Magnetic dipole moment of the isotope in nuclear magnetons.
    + Q                  ::Float64           ... Electric quadrupole moment of the isotope.
"""
struct  Isotope 
    massNumber           ::Float64     
    isotopeMass          ::Float64  
    spin                 ::AngularJ64
    mu                   ::Float64   
    Q                    ::Float64 
end


"""
`PeriodicTable.Isotope()`  ... simple constructur for an empty instance of an isotope.
"""
function Isotope(i::Int64)
    Isotope( i, NaN, AngularJ64(0), NaN, NaN)
end


"""
`struct  PeriodicTable.Element`  ... defines a type for collecting data for each element.

    + nuclearCharge     ::Int64                 ... Nuclear charge
    + symbol            ::Symbol                ... Symbol of the element.
    + name              ::String                ... Name of the element.
    + mass              ::Float64               ... Mean mass number of the element.
    + density           ::Float64               ... Density of the element [g/cm^3].
    + firstIP           ::Float64               ... First ionization potential [eV].
    + configuration     ::String                ... Ground configuration of the element.
    + atomicRadius      ::Float64               ... Mean atomic radius of the neutral atom (Ghosh & Biswas, 2002).
    + polarizibility    ::Float64               ... Mean polarizibility of the neutral atom (Ghosh & Biswas, 2002).
    + chemHardness      ::Float64               ... Mean chemical hardness of the element (Ghosh & Biswas, 2002).
"""
struct Element
    nuclearCharge       ::Int64 
    symbol              ::Symbol
    name                ::String
    mass                ::Float64
    density             ::Float64 
    firstIP             ::Float64  
    configuration       ::String
    atomicRadius        ::Float64  
    polarizibility      ::Float64 
    chemHardness        ::Float64 
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
    ... provides all available data for an element with given nuclear charge Z::Int64 by PeriodicTable.data[Z];
    cf. PeriodicTable.Element
"""
const data = Dict{Int64,PeriodicTable.Element}(
            1 => Element(  1, :H,  "hydrogen",         1.008,   0.09,  13.60,  "1s",			            0.5292,     0.6669,   13.588    ),
            2 => Element(  2, :He, "helium",           4.003,   0.18,  24.59,  "1s^2",			            0.3113,     0.1358,   22.383    ),
            3 => Element(  3, :Li, "lithum",           6.941,   0.53,   5.39,  "[He] 2s",			        1.6282,    19.4239,    4.4164   ),
            4 => Element(  4, :Be, "beryllium",        9.012,   1.85,   9.32,  "[He] 2s^2", 		        1.0855,     5.7558,    6.6244   ),
            5 => Element(  5, :B,  "boron",           10.811,   2.34,   8.30,  "[He] 2s^2 2p",		        0.8141,     2.428 ,    8.8328   ),
            6 => Element(  6, :C,  "carbon",          12.011,   2.26,  11.26,  "[He] 2s^2 2p^2",		    0.6513,     1.2432,   11.0407   ),
            7 => Element(  7, :N,  "nitrogen",        14.007,   1.25,  14.53,  "[He] 2s^2 2p^3",		    0.5427,     0.7193,   13.25     ),
            8 => Element(  8, :O,  "oxygen",          15.999,   1.43,  13.62,  "[He] 2s^2 2p^4",		    0.4652,     0.453 ,   15.4574   ),
            9 => Element(  9, :F,  "fluorine",        18.998,   1.70,  17.42,  "[He] 2s^2 2p^5",		    0.4071,     0.3034,   17.6634   ),
            10 => Element( 10, :Ne, "neon",            20.180,   0.90,  21.56,  "[Ne]",			            0.3618,     0.2131,   19.875    ),
            11 => Element( 11, :Na, "sodium",          22.990,   0.97,   5.14,  "[Ne] 3s",			        2.1649,    45.659 ,    3.3215   ),
            12 => Element( 12, :Mg, "magnesium",       24.305,   1.74,   7.65,  "[Ne] 3s^2", 		        1.6711,    21.0   ,    4.303    ),
            13 => Element( 13, :Al, "aluminium",       26.982,   2.70,   5.99,  "[Ne] 3s^2 3p",		        1.3607,    11.337 ,    5.2846   ),
            14 => Element( 14, :Si, "silicon",         28.086,   2.33,   8.15,  "[Ne] 3s^2 3p^2",		    1.1476,     6.8011,    6.2659   ),
            15 => Element( 15, :P,  "phosphorus",      30.974,   1.82,  10.49,  "[Ne] 3s^2 3p^3",		    0.9922,     4.3955,    7.2473   ),
            16 => Element( 16, :S,  "sulfur",          32.065,   2.07,  10.36,  "[Ne] 3s^2 3p^4",		    0.8738,     3.0023,    8.2293   ),
            17 => Element( 17, :Cl, "chlorine",        35.453,   3.21,  12.97,  "[Ne] 3s^2 3p^5",		    0.7807,     2.1412,    9.2107   ),
            18 => Element( 18, :Ar, "argon",           39.948,   1.78,  15.76,  "[Ar]",			            0.7056,     1.5808,    10.191   ),
            19 => Element( 19, :K,  "potassium",       39.098,   0.86,   4.34,  "[Ar] 4s",			        3.5598,   202.9969,     2.02    ),
            20 => Element( 20, :Ca, "calcium",         40.078,   1.55,   6.11,  "[Ar] 4s^2", 		        2.7479,    93.3717,     2.6162  ),
            21 => Element( 21, :Sc, "scandium",        44.956,   2.99,   6.56,  "[Ar] 3d 4s^2",		        2.6106,    80.0633,     2.7545  ),
            22 => Element( 22, :Ti, "titanium",        47.867,   4.54,   6.83,  "[Ar] 3d^2 4s^2",		    2.4861,    69.1712,     2.8924  ),
            23 => Element( 23, :V,  "vanadium",        50.942,   6.11,   6.75,  "[Ar] 3d^3 4s^2",		    2.3732,    60.1472,     3.03    ),
            24 => Element( 24, :Cr, "chromium",        51.996,   7.19,   6.77,  "[Ar] 3d^5 4s",		        2.2701,    52.6438,     3.1676  ),
            25 => Element( 25, :Mn, "manganese",       54.938,   7.43,   7.43,  "[Ar] 3d^5 4s^2",		    2.1754,    46.3265,     3.3055  ),
            26 => Element( 26, :Fe, "iron",            55.845,   7.87,   7.90,  "[Ar] 3d^6 4s^2",		    2.0885,    40.9939,     3.443   ),
            27 => Element( 27, :Co, "cobalt",          58.933,   8.90,   7.88,  "[Ar] 3d^7 4s^2",		    2.008 ,    36.4555,     3.5811  ),
            28 => Element( 28, :Ni, "nickel",          58.693,   8.90,   7.64,  "[Ar] 3d^8 4s^2",		    1.9337,    32.5372,     3.7187  ),
            29 => Element( 29, :Cu, "copper",          63.546,   8.96,   7.73,  "[Ar] 3d^10 4s",		        1.8648,    29.1769,     3.8561  ),
            30 => Element( 30, :Zn, "zinc",            65.390,   7.13,   9.39,  "[Ar] 3d^10 4s^2",		    1.8004,    26.2615,     3.994   ),
            31 => Element( 31, :Ga, "gallium",         69.723,   5.91,   6.00,  "[Ar] 3d^10 4s^2 4p",	    1.5663,    17.2917,     4.5909  ),
            32 => Element( 32, :Ge, "germanium",       72.640,   5.32,   7.90,  "[Ar] 3d^10 4s^2 4p^2",	    1.3862,    11.9864,     5.1874  ),
            33 => Element( 33, :As, "arsenic",         74.922,   5.72,   9.79,  "[Ar] 3d^10 4s^2 4p^3",	    1.2431,     8.6443,     5.7846  ),
            34 => Element( 34, :Se, "selenium",        78.960,   4.79,   9.75,  "[Ar] 3d^10 4s^2 4p^4",	    1.1269,     6.4397,     6.381   ),
            35 => Element( 35, :Br, "bromine",         79.904,   3.12,  11.81,  "[Ar] 3d^10 4s^2 4p^5",	    1.0305,     4.9244,     6.978   ),
            36 => Element( 36, :Kr, "krypton",         83.800,   3.75,  14.00,  "[Kr]",			            0.9493,     3.8497,     7.5748  ),
            37 => Element( 37, :Rb, "rubidium",        85.468,   1.63,   4.18,  "[Kr] 5s",			        4.8106,   500.9683,     1.4948  ),
            38 => Element( 38, :Sr, "strontium",       87.620,   2.54,   5.69,  "[Kr] 5s^2", 		        3.7135,   230.4462,     1.9364  ),
            39 => Element( 39, :Y,  "yttrium",         88.906,   4.47,   6.22,  "[Kr] 4d 5s^2",		        3.5278,   197.5715,     2.0383  ),
            40 => Element( 40, :Zr, "zirconium",       91.224,   6.51,   6.63,  "[Kr] 4d^2 5s^2",		    3.3598,   170.6683,     2.1402  ),
            41 => Element( 41, :Nb, "niobium",         92.906,   8.57,   6.76,  "[Kr] 4d^4 5s",		        3.2071,   148.4397,     2.2422  ),
            42 => Element( 42, :Mo, "molybdenum",      95.940,  10.22,   7.09,  "[Kr] 4d^5 5s",		        3.0677,   129.9126,     2.344   ),
            43 => Element( 43, :Tc, "technetium",      98.000,  11.50,   7.28,  "[Kr] 4d^5 5s^2",		    2.9398,   114.3315,     2.446   ),
            44 => Element( 44, :Ru, "ruthenium",      101.070,  12.37,   7.36,  "[Kr] 4d^7 5s",		        2.8222,   101.1523,     2.5479  ),
            45 => Element( 45, :Rh, "rhodium",        102.906,  12.41,   7.46,  "[Kr] 4d^8 5s",		        2.7137,    89.9286,     2.6498  ),
            46 => Element( 46, :Pd, "palladium",      106.420,  12.02,   8.34,  "[Kr] 4d^10",		        2.6132,    80.3028,     2.7517  ),
            47 => Element( 47, :Ag, "silver",         107.868,  10.50,   7.58,  "[Kr] 4d^10 5s",		        2.5199,    72.005 ,     2.8536  ),
            48 => Element( 48, :Cd, "cadmium",        112.411,   8.65,   8.99,  "[Kr] 4d^10 5s^2",		    2.433 ,    64.8095,     2.9555  ),
            49 => Element( 49, :In, "indium",         114.818,   7.31,   5.79,  "[Kr] 4d^10 5s^2 5p",	    2.1167,    42.6767,     3.3972  ),
            50 => Element( 50, :Sn, "tin",            118.710,   7.31,   7.34,  "[Kr] 4d^10 5s^2 5p^2",	    1.8732,    29.5777,     3.8388  ),
            51 => Element( 51, :Sb, "antimony",       121.760,   6.68,   8.61,  "[Kr] 4d^10 5s^2 5p^3",	    1.6799,    21.3335,     4.2805  ),
            52 => Element( 52, :Te, "tellurium",      127.600,   6.24,   9.01,  "[Kr] 4d^10 5s^2 5p^4",	    1.5228,    15.8906,     4.7221  ),
            53 => Element( 53, :I,  "iodine",         126.905,   4.93,  10.45,  "[Kr] 4d^10 5s^2 5p^5",	    1.3926,    12.1532,     5.0636  ),
            54 => Element( 54, :Xe, "xenon",          131.293,   5.90,  12.13,  "[Xe]",			            1.2828,     9.4993,     5.6056  ),
            55 => Element( 55, :Cs, "caesium",        132.906,   1.87,   3.89,  "[Xe] 6s",			        6.0615,  1002.1964,     1.1863  ),
            56 => Element( 56, :Ba, "barium",         137.327,   3.59,   5.21,  "[Xe] 6s^2", 		        4.6788,   460.9098,     1.5369  ),
            57 => Element( 57, :La, "lanthanum",      138.906,   6.15,   5.58,  "[Xe] 5d 6s^2",		        3.8102,   248.9177,     1.8873  ),
            58 => Element( 58, :Ce, "cerium",         140.116,   6.77,   5.54,  "[Xe] 4f 5d 6s^2",		    3.2133,   149.2883,     2.2378  ),
            59 => Element( 59, :Pr, "praseodymium",   140.908,   6.77,   5.47,  "[Xe] 4f^3 6s^2",		    2.778 ,    96.4738,     2.5885  ),
            60 => Element( 60, :Nd, "neodymium",      144.240,   7.01,   5.53,  "[Xe] 4f^4 6s^2",		    2.4468,    65.9186,     2.9389  ),
            61 => Element( 61, :Pm, "promethium",     145.000,   7.30,   5.58,  "[Xe] 4f^5 6s^2",		    2.1861,    47.0135,     3.2893  ),
            62 => Element( 62, :Sm, "samarium",       150.360,   7.52,   5.64,  "[Xe] 4f^6 6s^2",		    1.9756,    34.6984,     3.6398  ),
            63 => Element( 63, :Eu, "europium",       151.964,   5.24,   5.67,  "[Xe] 4f^7 6s^2",		    1.802 ,    26.3316,     3.9905  ),
            64 => Element( 64, :Gd, "gadolinium",     157.250,   7.90,   6.15,  "[Xe] 4f^7 5d 6s^2", 	    1.6565,    20.4544,     4.341   ),
            65 => Element( 65, :Tb, "terbium",        158.925,   8.23,   5.86,  "[Xe] 4f^9 6s^2",		    1.5328,    16.2057,     4.6913  ),
            66 => Element( 66, :Dy, "dysprosium",     162.500,   8.55,   5.94,  "[Xe] 4f^10 6s^2",		    1.4262,    13.0543,     5.0419  ),
            67 => Element( 67, :Ho, "holmium",        164.930,   8.80,   6.02,  "[Xe] 4f^11 6s^2",		    1.3335,    10.6707,     5.3924  ),
            68 => Element( 68, :Er, "erbium",         167.259,   9.07,   6.11,  "[Xe] 4f^12 6s^2",		    1.2521,     8.8334,     5.743   ),
            69 => Element( 69, :Tm, "thulium",        168.934,   9.32,   6.18,  "[Xe] 4f^13 6s^2",		    1.1801,     7.3955,     6.0934  ),
            70 => Element( 70, :Yb, "ytterbium",      173.040,   6.90,   6.25,  "[Xe] 4f^14 6s^2",		    1.1159,     6.253 ,     6.4439  ),
            71 => Element( 71, :Lu, "lutetium",       174.967,   9.84,   5.43,  "[Xe] 4f^14 5d 6s^2",	    1.0583,     5.3338,     6.7947  ),
            72 => Element( 72, :Hf, "hafnium",        178.490,  13.31,   6.83,  "[Xe] 4f^14 5d^2 6s^2",	    1.0079,     4.6075,     7.1344  ),
            73 => Element( 73, :Ta, "tantalum",       180.948,  16.65,   7.55,  "[Xe] 4f^14 5d^3 6s^2",	    0.9594,     3.9739,     7.4951  ),
            74 => Element( 74, :W,  "tungsten",       183.840,  19.35,   7.86,  "[Xe] 4f^14 5d^4 6s^2",	    0.9165,     3.4643,     7.8459  ),
            75 => Element( 75, :Re, "rhenium",        186.207,  21.04,   7.83,  "[Xe] 4f^14 5d^5 6s^2",	    0.8773,     3.0385,     8.1965  ),
            76 => Element( 76, :Os, "osmium",         190.230,  22.60,   8.44,  "[Xe] 4f^14 5d^6 6s^2",	    0.8413,     2.6796,     8.5472  ),
            77 => Element( 77, :Ir, "iridium",        192.217,  22.40,   8.97,  "[Xe] 4f^14 5d^7 6s^2",	    0.8182,     2.3756,     8.8973  ),
            78 => Element( 78, :Pt, "platinum",       195.078,  21.45,   8.96,  "[Xe] 4f^14 5d^9 6s",	    0.7776,     2.1175,     9.2474  ),
            79 => Element( 79, :Au, "gold",           196.967,  19.32,   9.23,  "[Xe] 4f^14 5d^10 6s",	    0.7492,     1.8924,     9.598   ),
            80 => Element( 80, :Hg, "mercury",        200.590,  13.55,  10.44,  "[Xe] 4f^14 5d^10 6s^2",	    3.0636,   129.646 ,     2.3456  ),
            81 => Element( 81, :Tl, "thallium",       204.383,  11.85,   6.11,  "[Xe] 4f^14 5d^10 6s^2 6p",  2.667 ,    85.3653,     2.6962  ),
            82 => Element( 82, :Pb, "lead",           207.200,  11.35,   7.42,  "[Xe] 4f^14 5d^10 6s^2 6p^2",2.3603,    59.1717,     3.0466  ),
            83 => Element( 83, :Bi, "bismuth",        208.980,   9.75,   7.29,  "[Xe] 4f^14 5d^10 6s^2 6p^3",2.1167,    42.6767,     3.3972  ),
            84 => Element( 84, :Po, "polonium",       209.000,   9.30,   8.42,  "[Xe] 4f^14 5d^10 6s^2 6p^4",1.9187,    31.7858,     3.7477  ),
            85 => Element( 85, :At, "astatine",       210.000,    NaN,   9.30,  "[Xe] 4f^14 5d^10 6s^2 6p^5",1.7546,    24.3079,     4.0983  ),
            86 => Element( 86, :Rn, "radon",          222.000,   9.73,  10.75,  "[Rn]",			            1.6164,    19.0046,     4.4487  ),
            87 => Element( 87, :Fr, "francium",       223.000,    NaN,   4.07,  "[Rn] 7s",			        7.2404,  1705.0   ,     0.9913  ),
            88 => Element( 88, :Ra, "radium",         226.000,   5.50,   5.28,  "[Rn] 7s^2", 		        5.5887,   785.4977,     1.2867  ),
            89 => Element( 89, :Ac, "actinium",       227.000,  10.07,   5.17,  "[Rn] 6d 7s^2",		        5.3091,   673.4033,     1.3544  ),
            90 => Element( 90, :Th, "thorium",        232.038,  11.72,   6.31,  "[Rn] 6d^2 7s^2",		    5.0569,   581.6815,     1.4222  ),
            91 => Element( 91, :Pa, "protactinium",   231.036,  15.40,   5.89,  "[Rn] 5f^2 6d 7s^2", 	    3.7042,   228.7156,     1.9413  ),
            92 => Element( 92, :U,  "uranium",        238.029,  18.95,   6.19,  "[Rn] 5f^3 6d 7s^2", 	    3.2177,   149.9164,     2.2348  ),
            93 => Element( 93, :Np, "neptunium",      237.000,  20.20,   6.27,  "[Rn] 5f^4 6d 7s^2", 	    2.8443,   103.5473,     2.5281  ),
            94 => Element( 94, :Pu, "plutonium",      244.000,  19.84,   6.03,  "[Rn] 5f^6 7s^2",		    2.3596,    59.1266,     3.0473  ),
            95 => Element( 95, :Am, "americium",      243.000,  13.67,   5.97,  "[Rn] 5f^7 7s^2",		    2.1525,    44.8789,     3.3407  ),
            96 => Element( 96, :Cm, "curium",         247.000,  13.50,   5.99,  "[Rn] 5f^7 6d 7s^2", 	    2.1097,    42.2547,     3.4084  ),
            97 => Element( 97, :Bk, "berkelium",      247.000,  14.78,   6.20,  "[Rn] 5f^9 7s^2",		    1.8308,    27.5918,     3.9277  ),
            98 => Element( 98, :Cf, "californium",    251.000,  15.10,   6.28,  "[Rn] 5f^10 7s^2",		    1.7035,    22.2453,     4.2212  ),
            99 => Element( 99, :Es, "einsteinium",    252.000,    NaN,   6.42,  "[Rn] 5f^11 7s^2",		    1.5928,    18.1843,     4.5146  ),
            100 => Element(100, :Fm, "fermium",        257.000,    NaN,   6.50,  "[Rn] 5f^12 7s^2",		    1.4956,    15.0542,     4.808   ),
            101 => Element(101, :Md, "mendelevium",    258.000,    NaN,   6.58,  "[Rn] 5f^13 7s^2",		    1.4096,    12.6038,     5.1013  ),
            102 => Element(102, :No, "nobelium",       259.000,    NaN,   6.65,  "[Rn] 5f^14 7s^2",		    1.3329,    10.6563,     5.3949  ),
            103 => Element(103, :Lr, "lawrencium",     262.000,    NaN,   4.90,  "[Rn] 5f^14 7s^2 7p",	    1.3164,    10.2631,     5.4629  ),
            104 => Element(104, :Rf, "rutherfordium",  261.000,    NaN,    NaN,  "[Rn] 5f^14 6d^2 7s^2",	       NaN,        NaN,        NaN  ),
            105 => Element(105, :Db, "dubnium",        262.000,    NaN,    NaN,  "[Rn] 5f^14 6d^3 7s^2",	       NaN,        NaN,        NaN  ),
            106 => Element(106, :Sg, "seaborgium",     266.000,    NaN,    NaN,  "[Rn] 5f^14 6d^4 7s^2",	       NaN,        NaN,        NaN  ),
            107 => Element(107, :Bh, "bohrium",        264.000,    NaN,    NaN,  "[Rn] 5f^14 6d^5 7s^2",	       NaN,        NaN,        NaN  ),
            108 => Element(108, :Hs, "hassium",        277.000,    NaN,    NaN,  "[Rn] 5f^14 6d^6 7s^2",	       NaN,        NaN,        NaN  ),
            109 => Element(109, :Mt, "meitnerium",     268.000,    NaN,    NaN,  "[Rn] 5f^14 6d^7 7s^2",	       NaN,        NaN,        NaN  ),
            110 => Element(110, :Ds, "darmstadtium",       NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^8 7s^2",	       NaN,        NaN,        NaN  ),
            111 => Element(111, :Rg, "roentgenium",        NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^9 7s^2",	       NaN,        NaN,        NaN  ),
            112 => Element(112, :Cn, "copernicium",        NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^10 7s^2",	       NaN,        NaN,        NaN  ),
            113 => Element(113, :Nh, "nihonium",           NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^8 7s^2 7p",	   NaN,        NaN,        NaN  ),
            114 => Element(114, :Fl, "flerovium",          NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^2 7s^2 7p^2",    NaN,        NaN,        NaN  ),
            115 => Element(115, :Mc, "moscovium",          NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^3 7s^2 7p^3",    NaN,        NaN,        NaN  ),
            116 => Element(116, :Lv, "livermorium",        NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^4 7s^2 7p^4",    NaN,        NaN,        NaN  ),
            117 => Element(117, :Ts, "tennessine",         NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^5 7s^2 7p^5",    NaN,        NaN,        NaN  ),
            118 => Element(118, :Og, "oganesson",          NaN,    NaN,    NaN,  "[Rn] 5f^14 6d^6 7s^2 7p^6",    NaN,        NaN,        NaN  )	 )
                                                        
                                                            
"""
`PeriodicTable.isotopes::Dict{Int64,Array{PeriodicTable.Isotope,1}}`
    ... provides isotope data for an element with given nuclear charge Z::Int64 from PeriodicTable.isotopes[Z];
        cf. PeriodicTable.Isotope
"""
const isotopes = Dict{Int64,Array{PeriodicTable.Isotope,1}}(
            1 => Isotope[],
            2 => Isotope[ Isotope( 2, 2.015894,      AngularJ64(0), NaN, NaN),    Isotope( 2, 3.01602932265, AngularJ64(1//2), NaN, NaN), 
                            Isotope( 4, 4.00260325413, AngularJ64(0), NaN, NaN)],
            3 => Isotope[],
            4 => Isotope[],
            5 => Isotope[],
            6 => Isotope[],
            7 => Isotope[],
            8 => Isotope[],
            9 => Isotope[],
            10 => Isotope[],
            11 => Isotope[],
            12 => Isotope[],
            13 => Isotope[],
            14 => Isotope[],
            15 => Isotope[],
            16 => Isotope[],
            17 => Isotope[],
            18 => Isotope[],
            19 => Isotope[],
            20 => Isotope[],
            21 => Isotope[],
            22 => Isotope[],
            23 => Isotope[],
            24 => Isotope[],
            25 => Isotope[],
            26 => Isotope[],
            27 => Isotope[],
            28 => Isotope[],
            29 => Isotope[],
            30 => Isotope[],
            31 => Isotope[],
            32 => Isotope[],
            33 => Isotope[],
            34 => Isotope[],
            35 => Isotope[],
            36 => Isotope[],
            37 => Isotope[],
            38 => Isotope[],
            39 => Isotope[],
            40 => Isotope[],
            41 => Isotope[],
            42 => Isotope[],
            43 => Isotope[],
            44 => Isotope[],
            45 => Isotope[],
            46 => Isotope[],
            47 => Isotope[],
            48 => Isotope[],
            49 => Isotope[],
            50 => Isotope[], 
            51 => Isotope[],
            52 => Isotope[],
            53 => Isotope[],
            54 => Isotope[],
            55 => Isotope[],
            56 => Isotope[],
            57 => Isotope[],
            58 => Isotope[],
            59 => Isotope[],
            60 => Isotope[],
            61 => Isotope[],
            62 => Isotope[],
            63 => Isotope[],
            64 => Isotope[],
            65 => Isotope[],
            66 => Isotope[],
            67 => Isotope[],
            68 => Isotope[],
            69 => Isotope[],
            70 => Isotope[],
            71 => Isotope[],
            72 => Isotope[],
            73 => Isotope[],
            74 => Isotope[],
            75 => Isotope[],
            76 => Isotope[],
            77 => Isotope[],
            78 => Isotope[],
            79 => Isotope[],
            80 => Isotope[],
            81 => Isotope[],
            82 => Isotope[],
            83 => Isotope[],
            84 => Isotope[],
            85 => Isotope[],
            86 => Isotope[],
            87 => Isotope[],
            88 => Isotope[],
            89 => Isotope[],
            90 => Isotope[],
            91 => Isotope[],
            92 => Isotope[],
            93 => Isotope[],
            94 => Isotope[],
            95 => Isotope[],
            96 => Isotope[],
            97 => Isotope[],
            98 => Isotope[],
            99 => Isotope[],
            100 => Isotope[],
            101 => Isotope[],
            102 => Isotope[],
            103 => Isotope[],
            104 => Isotope[],
            105 => Isotope[],
            106 => Isotope[],
            107 => Isotope[],
            108 => Isotope[],
            109 => Isotope[],
            110 => Isotope[],
            111 => Isotope[],
            112 => Isotope[],
            113 => Isotope[],
            114 => Isotope[],
            115 => Isotope[],
            116 => Isotope[],
            117 => Isotope[],
            118 => Isotope[] )


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
    elseif Z ==  6 # wa = [   284.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
                        wa = [   284.2,     19.4,      7.0,      7.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] #ZS
    elseif Z ==  7   wa = [   409.9,     37.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  8   wa = [   543.1,     41.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
    elseif Z ==  9 # wa = [   696.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0]
                        wa = [   696.7,     31.0,      9.0,      9.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0] #ZS
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


"""
`PeriodicTable.bindingEnergies_XrayDataBooklet(Z::Int64)` 
    ... to return the `stored values' of Bearden et al. (1967), Cardona et al. (1978) and Fuggle et al. (1980) 
        for the element with nuclear charge Z; all binding energies are given in eV.
"""
function bindingEnergies_XrayDataBooklet(Z::Int64)
    # 
    #                        1s_1/2,   2s_1/2,   2p_1/2,   2p_3/2,   3s_1/2,   3p_1/2,   3p_3/2,   3d_3/2,   3d_5/2,   4s_1/2,   4p_1/2,   4p_3/2, 
    #                              4d_3/2,   4d_5/2,   4f_5/2,  4f_7/2,   5s_1/2,   5p_1/2,   5p_3/2,   5d_3/2,   5d_5/2,   6s_1/2,   6p_1/2,    6p_3/2
    if     Z ==  1   wa = [    13.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  2   wa = [    24.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  3   wa = [    54.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  4   wa = [   111.5,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  5   wa = [   188.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  6   wa = [   284.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  7   wa = [   409.9,     37.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  8   wa = [   543.1,     41.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z ==  9   wa = [   696.7,     31.0,      9.0,      9.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 10   wa = [   870.2,     48.5,     21.7,     21.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 11   wa = [  1070.8,     63.5,    30.65,    30.81,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 12   wa = [  1303.0,     88.7,    49.78,    49.50,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0, 
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 13   wa = [  1559.6,    117.8,    72.95,    72.55,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 14   wa = [  1839.0,    149.7,    99.82,    99.42,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 15   wa = [  2145.5,    189.0,   136.0,    135.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 16   wa = [  2472.0,    230.9,   163.6,    162.5,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 17   wa = [  2822.4,    270.0,   202.0,    200.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 18   wa = [  3205.9,    326.3,   250.6,    248.4,      29.3,     15.9,     15.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 19   wa = [  3608.4,    378.6,   297.3,    294.6,      34.8,     18.3,     18.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 20   wa = [  4038.5,    438.4,   349.7,    346.2,      44.3,     25.4,     25.4,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 21   wa = [  4492.0,    498.0,   403.6,    398.7,      51.1,     28.3,     28.3,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 22   wa = [  4966.0,    560.9,   460.2,    453.8,      58.7,     32.6,     32.6,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 23   wa = [  5465.0,    626.7,   519.8,    512.1,      66.3,     37.2,     37.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 24   wa = [  5989.0,    696.0,   583.8,    574.1,      74.1,     42.2,     42.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 25   wa = [  6539.0,    769.1,   649.9,    638.7,      82.3,     47.2,     47.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 26   wa = [  7112.0,    844.6,   719.9,    706.8,      91.3,     52.7,     52.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 27   wa = [  7709.0,    925.1,   793.2,    778.1,     101.0,     58.9,     59.9,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 28   wa = [  8333.0,   1008.6,   870.0,    852.7,     110.8,     68.0,     66.2,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 29   wa = [  8979.0,   1096.7,   952.3,    932.7,     122.5,     77.3,     75.1,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 30   wa = [  9659.0,   1196.2,  1044.9,   1021.8,     139.8,     91.4,     88.6,     10.2,     10.1,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 31   wa = [ 10367.0,   1299.0,  1143.2,   1116.4,     159.5,    103.5,    100.0,     18.7,     18.7,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 32   wa = [ 11103.0,   1414.6,  1248.1,   1217.0,     180.1,    124.9,    120.8,     29.8,     29.2,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 33   wa = [ 11867.0,   1527.0,  1359.1,   1323.6,     204.7,    146.2,    141.2,     41.7,     41.7,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 34   wa = [ 12658.0,   1652.0,  1474.3,   1433.9,     229.6,    166.5,    160.7,     55.5,     54.6,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 35   wa = [ 13474.0,   1782.0,  1596.0,   1550.0,     257.0,    189.0,    182.0,     70.0,     69.0,     -1.0,     -1.0,     -1.0,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 36   wa = [ 14326.0,   1921.0,  1730.9,   1678.4,     292.8,    222.2,    214.4,     95.0,     93.8,     27.5,     14.1,     14.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    #
    elseif Z == 37   wa = [ 15200.0,   2065.0,  1864.0,   1804.0,     326.7,    248.7,    239.1,    113.0,    112.0,     30.5,     16.3,     15.3,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 38   wa = [ 16105.0,   2216.0,  2007.0,   1940.0,     358.7,    280.3,    270.0,    136.0,    134.2,     38.9,     21.3,     20.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 39   wa = [ 17038.0,   2373.0,  2156.0,   2080.0,     392.0,    310.6,    298.8,    157.7,    155.8,     43.8,     24.4,     23.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 40   wa = [ 17998.0,   2532.0,  2307.0,   2223.0,     430.3,    343.5,    329.8,    181.1,    178.8,     50.6,     28.5,     27.1,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 41   wa = [ 18986.0,   2698.0,  2465.0,   2371.0,     466.6,    376.1,    360.6,    205.0,    202.3,     56.4,     32.6,     30.8,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 42   wa = [ 20000.0,   2866.0,  2625.0,   2520.0,     506.3,    411.6,    394.0,    231.1,    227.9,     63.2,     37.6,     35.5,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 43   wa = [ 21044.0,   3043.0,  2793.0,   2677.0,     544.0,    447.6,    417.7,    257.6,    253.9,     69.5,     42.3,     39.9,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 44   wa = [ 22117.0,   3224.0,  2967.0,   2838.0,     586.1,    483.5,    461.4,    284.2,    280.0,     75.0,     46.3,     43.2,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 45   wa = [ 23220.0,   3412.0,  3146.0,   3004.0,     628.1,    521.3,    496.5,    311.9,    307.2,     81.4,     50.5,     47.3,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 46   wa = [ 24350.0,   3604.0,  3330.0,   3173.0,     671.6,    599.9,    532.3,    340.5,    335.2,     87.1,     55.7,     50.9,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]
    elseif Z == 47   wa = [ 25514.0,   3806.0,  3524.0,   3351.0,     719.0,    603.8,    573.0,    374.0,    368.3,     97.0,     63.7,     58.3,
                                    -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0 ]  
    #
    elseif Z == 48   wa = [26711.0,    4018.0,  3727.0,   3538.0,     772.0,    652.6,    618.4,    411.9,    405.2,    109.8,     63.9,     63.9,
                                    11.7,     10.7,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 49   wa = [27940.0,    4238.0,  3938.0,   3730.0,     827.2,    703.2,    665.3,    451.4,    443.9,    122.9,     73.5,     73.5,
                                    17.7,     16.9,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 50   wa = [29200.0,    4465.0,  4156.0,   3929.0,     884.7,    756.5,    714.6,    493.2,    484.9,    137.1,     83.6,     83.6,
                                    24.9,     23.9,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 51   wa = [30491.0,    4698.0,  4380.0,   4132.0,     946.0,    812.7,    766.4,    537.5,    528.2,    153.2,     95.6,     95.6,
                                    33.3,     32.1,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 52   wa = [31814.0,    4939.0,  4612.0,   4341.0,    1006.0,    870.8,    820.0,    583.4,    573.0,    169.4,    103.3,    103.3,
                                    41.9,     40.4,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 53   wa = [33169.0,    5188.0,  4852.0,   4557.0,    1072.0,    931.0,    875.0,    630.8,    619.3,    186.0,    123.0,    123.0,
                                    50.6,     48.9,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 54   wa = [34561.0,    5453.0,  5107.0,   4786.0,    1148.7,   1002.1,    940.6,    689.0,    676.4,    213.2,    146.7,    145.5,
                                    69.5,     67.5,    -1.0,     -1.0,      23.3,     13.4,     12.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 55   wa = [35985.0,    5714.0,  5359.0,   5012.0,    1211.0,   1071.0,   1003.0,    740.5,    726.6,    232.3,    172.4,    161.3,
                                    79.8,     77.5,    -1.0,     -1.0,      22.7,     14.2,     12.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 56   wa = [37441.0,    5989.0,  5624.0,   5247.0,    1293.0,   1137.0,   1063.0,    795.7,    780.5,    253.5,    192.0,    178.6,
                                    92.6,     89.9,    -1.0,     -1.0,      30.3,     17.0,     14.8,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0] 
    elseif Z == 57   wa = [38925.0,    6266.0,  5891.0,   5483.0,    1362.0,   1209.0,   1128.0,    853.0,    836.0,    274.7,    205.8,    196.0,
                                    105.3,    102.5,    -1.0,     -1.0,      34.3,     19.3,     16.8,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 58   wa = [40443.0,    6549.0,  6164.0,   5723.0,    1436.0,   1274.0,   1187.0,    902.4,    883.8,    291.0,    223.2,    206.5,
                                    109.0,     -1.0,     0.1,      0.1,      37.8,     19.8,     17.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 59   wa = [41991.0,    6835.0,  6440.0,   5964.0,    1511.0,   1337.0,   1242.0,    948.3,    928.8,    304.5,    236.3,    217.6,
                                    115.1,    115.1,     2.0,      2.0,      37.4,     22.3,     22.3,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 60   wa = [43569.0,    7126.0,  6722.0,   6208.0,    1575.0,   1403.0,   1297.0,   1003.3,    980.4,    319.2,    243.3,    224.6,
                                    120.5,    120.5,     1.5,      1.5,      37.5,     21.1,     21.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 61   wa = [45184.0,    7428.0,  7013.0,   6459.0,      -1.0,   1471.0,   1357.0,   1052.0,   1027.0,     -1.0,    242.0,    242.0,
                                    120.0,    120.0,    -1.0,     -1.0,      -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 62   wa = [46834.0,    7737.0,  7312.0,   6716.0,    1723.0,   1541.0,   1420.0,   1110.9,   1083.4,    347.2,    265.6,    247.4,
                                    129.0,    129.0,     5.2,      5.2,      37.4,     21.3,     21.3,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 63   wa = [48519.0,    8052.0,  7617.0,   6977.0,    1800.0,   1614.0,   1481.0,   1158.6,   1127.5,    360.0,    284.0,    257.0,
                                    133.0,    127.7,     0.0,      0.0,      32.0,     22.0,     22.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 64   wa = [50239.0,    8376.0,  7930.0,   7243.0,    1881.0,   1688.0,   1544.0,   1221.9,   1189.6,    378.6,    286.0,    271.0,
                                    -1.0,    142.6,     8.6,      8.6,      36.0,     28.0,     21.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 65   wa = [51996.0,    8708.0,  8252.0,   7514.0,    1968.0,   1768.0,   1611.0,   1276.9,   1241.1,    396.0,    322.4,    284.1,
                                    150.5,    150.5,     7.7,      2.4,      45.6,     28.7,     22.6,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 66   wa = [53789.0,    9046.0,  8581.0,   7790.0,    2047.0,   1842.0,   1676.0,   1333.0,   1292.6,    414.2,    333.5,    293.2,
                                    153.6,    153.6,     8.0,      4.3,      49.9,     26.3,     26.3,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 67   wa = [55618.0,    9394.0,  8918.0,   8071.0,    2128.0,   1923.0,   1741.0,   1392.0,   1351.0,    432.4,    343.5,    308.2,
                                    160.0,    160.0,     8.6,      5.2,      49.3,     30.8,     24.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 68   wa = [57486.0,    9751.0,  9264.0,   8358.0,    2207.0,   2006.0,   1812.0,   1453.0,   1409.0,    449.8,    366.2,    320.2,
                                    167.6,    167.6,    -1.0,      4.7,      50.6,     31.4,     24.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 69   wa = [59390.0,   10116.0,  9617.0,   8648.0,    2307.0,   2090.0,   1885.0,   1515.0,   1468.0,    470.9,    385.9,    332.6,
                                    175.5,    175.5,    -1.0,      4.6,      54.7,     31.8,     25.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 70   wa = [61332.0,   10486.0,  9978.0,   8944.0,    2398.0,   2173.0,   1950.0,   1576.0,   1528.0,    480.5,    388.7,    339.7,
                                    191.2,    182.4,     2.5,      1.3,      52.0,     30.3,     24.1,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    #
    elseif Z == 71   wa = [63314.0,   10870.0, 10349.0,   9244.0,    2491.0,   2264.0,   2024.0,   1639.0,   1589.0,    506.8,    412.4,    359.2,
                                    206.1,    196.3,    89.0,      7.5,      57.3,     33.6,     26.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 72   wa = [65351.0,   11271.0, 10739.0,   9561.0,    2601.0,   2365.0,   2108.0,   1716.0,   1662.0,    538.0,    438.2,    380.7,
                                    220.0,    211.5,    15.9,     14.2,      64.2,     38.0,     29.9,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 73   wa = [67416.0,   11682.0, 11136.0,   9881.0,    2708.0,   2469.0,   2194.0,   1793.0,   1735.0,    563.4,    463.4,    400.9,
                                    237.9,    226.4,    23.5,     21.6,      69.7,     42.2,     32.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]		
    elseif Z == 74   wa = [69525.0,   12100.0, 11544.0,  10207.0,    2820.0,   2575.0,   2281.0,   1872.0,   1809.0,    594.1,    490.4,    423.6,
                                    255.9,    243.5,    33.6,     31.4,      75.6,     45.3,     36.8,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 75   wa = [71676.0,   12527.0, 11959.0,  10535.0,    2932.0,   2682.0,   2367.0,   1949.0,   1883.0,    625.4,    518.7,    446.8,
                                    273.9,    260.5,    42.9,     40.5,      83.0,     45.6,     34.6,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 76   wa = [73871.0,   12968.0, 12385.0,  10871.0,    3049.0,   2792.0,   2457.0,   2031.0,   1960.0,    658.2,    549,1,    470.7,
                                    293.1,    278.5,    53.4,     50.7,      84.0,     58.0,     44.5,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 77   wa = [76111.0,   13419.0, 12824.0,  11215.0,    3174.0,   2909.0,   2551.0,   2116.0,   2040.0,    691.1,    577.8,    495.8,
                                    311.9,    296.3,    63.8,     60.8,      95.2,     63.0,     48.0,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 78   wa = [78395.0,   13880.0, 13273.0,  11564.0,    3296.0,   3027.0,   2645.0,   2202.0,   2122.0,    725.4,    609.1,    519.4,
                                    331.6,    314.6,    74.5,     71.2,     101.7,     65.3,     51.7,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 79   wa = [80725.0,   14353.0, 13734.0,  11919.0,    3425.0,   3148.0,   2743.0,   2291.0,   2206.0,    762.1,    642.7,    546.3,
                                    353.2,    335.1,    87.6,     84.0,     107.2,     74.2,     57.2,     -1.0,     -1.0,     -1.0,     -1.0,      -1.0]	
    elseif Z == 80   wa = [83102.0,   14839.0, 14209.0,  12284.0,    3562.0,   3279.0,   2847.0,   2385.0,   2295.0,    802.2,    680.2,    576.6,
                                    378.2,    358.8,   104.0,     99.9,     127.0,     83.1,     64.5,     96.0,      7.8,     -1.0,     -1.0,      -1.0]	
    elseif Z == 81   wa = [85530.0,   15347.0, 14698.0,  12658.0,    3704.0,   3416.0,   2957.0,   2485.0,   2389.0,    846.2,    720.5,    609.5,
                                    405.7,    385.0,   122.2,    117.8,     136.0,     94.6,     73.5,     14.7,     12.5,     -1.0,     -1.0,      -1.0]
    elseif Z == 82   wa = [88005.0,   15861.0, 15200.0,  13035.0,    3851.0,   3554.0,   3066.0,   2586.0,   2484.0,    891.8,    761.9,    643.5,
                                    434.3,    412.2,   141.7,    136.9,     147.0,    106.4,     83.3,     20.7,     18.1,     -1.0,     -1.0,      -1.0]
    elseif Z == 83   wa = [90524.0,   16388.0, 15711.0,  13419.0,    3999.0,   3696.0,   3177.0,   2688.0,   2580.0,    939.0,    805.2,    678.8,
                                    464.0,    440.1,   162.3,    157.0,     159.3,    119.0,     92.6,     26.9,     23.8,     -1.0,     -1.0,      -1.0]
    elseif Z == 84   wa = [93105.0,   16939.0, 16244.0,  13814.0,    4149.0,   3854.0,   3302.0,   2798.0,   2683.0,    995.0,    851.0,    705.0,
                                    500.0,    473.0,   184.0,    184.0,     177.0,    132.0,    104.0,     31.0,     31.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 85   wa = [95730.0,   17493.0, 16785.0,  14214.0,    4317.0,   4008.0,   3426.0,   2909.0,   2787.0,    1042.0,   886.0,    740.0,
                                    533.0,    507.0,   210.0,    210.0,     195.0,    148.0,    115.0,     40.0,     40.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 86   wa = [98404.0,   18049.0, 17337.0,  14619.0,    4482.0,   4159.0,   3538.0,   3022.0,   2892.0,    1097.0,   929.0,    768.0,
                                    567.0,    541.0,   238.0,    238.0,     214.0,    164.0,    127.0,     48.0,     48.0,     26.0,     -1.0,      -1.0]	
    elseif Z == 87   wa = [101137.0,  18639.0, 17907.0,  15031.0,    4652.0,   4327.0,   3663.0,   3136.0,   3000.0,    1153.0,   980.0,    810.0,
                                    603.0,    577.0,   268.0,    268.0,     234.0,    182.0,    140.0,     58.0,     58.0,     34.0,     15.0,      15.0]
    elseif Z == 88   wa = [103922.0,  19237.0, 18484.0,  15444.0,    4822.0,   4490.0,   3792.0,   3248.0,   3105.0,    1208.0,  1058.0,    879.0,
                                    636.0,    603.0,   299.0,    299.0,     254.0,    200.0,    153.0,     68.0,     68.0,     44.0,     19.0,      19.0]
    elseif Z == 89   wa = [106755.0,  19840.0, 19083.0,  15871.0,    5002.0,   4656.0,   3909.0,   3370.0,   3219.0,    1269.0,  1080.0,    890.0,
                                    675.0,    639.0,   319.0,    319.0,     272.0,    215.0,    167.0,     80.0,     80.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 90   wa = [109651.0,  20472.0, 19693.0,  16300.0,    5182.0,   4830.0,   4046.0,   3491.0,   3332.0,    1330.0,  1168.0,    966.4,
                                    712.1,    675.2,   342.4,    333.1,     290.0,    229.0,    182.0,     92.5,     85.4,     41.4,     24.5,      16.6]
    elseif Z == 91   wa = [112601.0,  21105.0, 20314.0,  16733.0,    5367.0,   5001.0,   4174.0,   3611.0,   3442.0,    1387.0,  1224.0,   1007.0,
                                    743.0,    708.0,   371.0,    360.0,     310.0,    232.0,    232.0,     94.0,     94.0,     -1.0,     -1.0,      -1.0]
    elseif Z == 92   wa = [115606.0,  21757.0, 20948.0,  17166.0,    5548.0,   5182.0,   4303.0,   3728.0,   3552.0,    1439.0,  1271.0,   1043.0, 
                                    778.3,    736.2,   388.2,    377.4,     321.0,    257.0,    192.0,    102.8,     94.2,     43.9,     26.8,      16.8] 
    #       
    else  error("Data not available for Z = $Z")
    end
    return( wa )
end

end # module

