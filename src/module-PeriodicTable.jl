
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

end # module

