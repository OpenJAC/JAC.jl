
"""
`module JAC.AngularCoefficientsRatip2013`  
    ... a submodel of JAC that provides data types for angular (T and V) coefficients and methods to calculate the angular 
        coefficients between pairs of relativistic configuration state functions via an external Fortran library.
"""
module AngularCoefficientsRatip2013

    ##x using ..JAC: AngularJ64, twice, CsfR, getsubshells, have_same_subshells, Parity, plus, minus, Subshell, SubshellQuantumNumbers
    ##x using ..JAC: AngularJ64, CsfR, Parity, plus, minus, Subshell # , SubshellQuantumNumbers
    using  ..Basics,  ..Defaults,  ..ManyElectron
    using  Libdl: dlopen, dlsym

    export AngularTcoeff, AngularVcoeff

    const ANCOPATH = joinpath(@__DIR__, "..", "deps", "bin", "libanco-ratip2013.so")


    """
    `AngularCoefficientsRatip2013.twice(j::AngularJ64)`  ... returns `2*j` for a given angular momentum quantum number `j`.
    """
    twice(j::AngularJ64) = Int(2*j.num//j.den)


    """
    `AngularCoefficientsRatip2013.getsubshells(csf::CsfR)` 
        ... returns the array of subshells of a relativistic configuration state function.
    """
    getsubshells(csf::CsfR) = csf.useStandardSubshells ?   Defaults.GBL_STANDARD_SUBSHELL_LIST : csf.subshells


    """
    `AngularCoefficientsRatip2013.have_same_subshells(csfa::CsfR, csfb::CsfR)`  
        ... checks whether two relativistic configurations state functions have the same subshells; true is returned if 
            both CSF have the same subshell list, and false otherwise.
    """
    function have_same_subshells(csfa::CsfR, csfb::CsfR)
        return csfa.useStandardSubshells && csfb.useStandardSubshells || getsubshells(csfa)==getsubshells(csfb)
    end


    """
    `struct AngularCoefficientsRatip2013.Fnkappa`  
        ... type that mirrors `nkappa` from the Fortran library. This type is used to communicate with the external Fortran 
            library. It can be converted from and to a `Subshell` using the appropriate constructors.
    """
    struct Fnkappa      # nkappa
         n    ::Int32   # integer
         kappa::Int32   # integer
    end

    ##x Fnkappa(s::Subshell) = begin qn = SubshellQuantumNumbers(string(s)); return Fnkappa(qn[1], qn[2]) end
    Fnkappa(s::Subshell) = Fnkappa(s.n, s.kappa)


    """
    `AngularCoefficientsRatip2013.Subshell(nk::Fnkappa)`  
        ... constructor to define a Subshell (instance) from the quantum numbers n, kappa.
    """
    Subshell(nk::Fnkappa) = Subshell(nk.n, nk.kappa)


    """
    `AngularCoefficientsRatip2013.maxocc(nk::Fnkappa)`  ... return the maximal occupation number for a given relativistic subshell.
    """
    maxocc(nk::Fnkappa) = 2abs(nk.kappa)


    """
    `struct AngularCoefficientsRatip2013.FancoTcoeff`  
        ... type that mirrors `anco_T_coeff` from the Fortran library. This type is used to communicate with the external 
            Fortran library. It can be converted to a `AngularTcoeff` using the appropriate constructor.
    """
    struct FancoTcoeff # anco_T_coeff
        nu::Int32      # integer
        a::Fnkappa     # nkappa
        b::Fnkappa     # nkappa
        T::Float64     # real(kind=dp)
    end


    """
    `struct AngularCoefficientsRatip2013.FancoVcoeff`  
        ... type that mirrors `anco_V_coeff` from the Fortran library. This type is used to communicate with the external 
            Fortran library. It can be converted to an `AngularVcoeff` using the appropriate constructor.
    """
    struct FancoVcoeff # anco_V_coeff
        nu::Int32      # integer
        a::Fnkappa     # nkappa
        b::Fnkappa     # nkappa
        c::Fnkappa     # nkappa
        d::Fnkappa     # nkappa
        V::Float64     # real(kind=dp)
    end


    """
    `struct AngularCoefficientsRatip2013.Fmtcoefficient`  
        ... type that mirrors `mct_coefficient` from the Fortran library. This type is used to communicate with the
            external Fortran library. It can be converted to an `AngularTcoeff` using the appropriate constructor.
    """
    struct Fmctcoefficient # mct_coefficient
        r::Int32           # integer
        s::Int32           # integer
        a::Int16           # integer(kind=i2b)
        b::Int16           # integer(kind=i2b)
        ν::Int16           # integer(kind=i2b)
        T::Float64         # real(kind=dp)
    end


    """
    `struct AngularTcoeff`  ... type for representing angular T coefficients.
    """
    struct AngularTcoeff
        nu::Int
        a::Subshell
        b::Subshell
        T::Float64
    end


    # `Base.show(io::IO, coeff::AngularTcoeff)`  ... prepares a proper printout of the variable coeff::AngularTcoeff.
    function Base.show(io::IO, coeff::AngularTcoeff) 
        sa = " T^$(coeff.nu)[" * string(coeff.a) * "," * string(coeff.b) * "] = $(coeff.T)"
        print(io, sa )
    end


    """
    `AngularCoefficientsRatip2013.AngularTcoeff(t::FancoTcoeff)`  
        ... constructor to define an AngularTcoeff (instance) from an internal FancoTcoeff coefficient.
    """
    AngularTcoeff(t::FancoTcoeff) = AngularTcoeff(t.nu, Subshell(t.a), Subshell(t.b), t.T)

    """
    `AngularCoefficientsRatip2013.AngularTcoeff(t::Fmctcoefficient)`  
        ... constructor to define an AngularTcoeff (instance) from an internal Fmctcoefficient and the list of relativistic 
            subshells the CSFs are based on.
    """
    AngularTcoeff(t::Fmctcoefficient, subshells::AbstractArray{Subshell}) = AngularTcoeff(t.ν, subshells[t.a], subshells[t.b], t.T)


    """
    `struct AngularCoefficientsRatip2013.AngularTcoeff`  ... type for representing angular V coefficients.
    """
    struct AngularVcoeff
        nu::Int
        a::Subshell
        b::Subshell
        c::Subshell
        d::Subshell
        V::Float64
    end


    # `Base.show(io::IO, coeff::AngularVcoeff)`  ... prepares a proper printout of the variable coeff::AngularVcoeff.
    function Base.show(io::IO, coeff::AngularVcoeff) 
        sa = " V^$(coeff.nu)[" * string(coeff.a) * "," * string(coeff.b) * ";" * string(coeff.c) * "," * string(coeff.d) * "] = $(coeff.V)"
        print(io, sa )
    end

    """
    `AngularCoefficientsRatip2013.AngularVcoeff(t::FancoVcoeff)`  
        ... constructor to define an AngularVcoeff (instance) from an internal FancoVcoeff coefficient.
    """
    AngularVcoeff(t::FancoVcoeff) = AngularVcoeff(t.nu, Subshell(t.a), Subshell(t.b), Subshell(t.c), Subshell(t.d), t.V)


    const NUMBER_OF_MCT_COEFFICIENTS_MAX = 1000

    # The following `Ref`s need to be initialized at runtime (via the __init__() function)
    # in order to enable precompilation

    const ANCOLIB  = Ref{Ptr{Nothing}}() # dlopen(ANCOPATH)

    const no_anco_t_list             = Ref{Ptr{Int32}}() # cglobal(dlsym(ANCOLIB[], :no_anco_t_list), Int32)
    const no_anco_v_list             = Ref{Ptr{Int32}}() # cglobal(dlsym(ANCOLIB[], :no_anco_v_list), Int32)
    const number_of_mct_coefficients = Ref{Ptr{Int32}}() # cglobal(dlsym(ANCOLIB[], :number_of_mct_coefficients), Int32)

    const anco_t_list = Ref{Ptr{FancoTcoeff}}()     # cglobal(dlsym(ANCOLIB[], :anco_t_list), FancoTcoeff)
    const anco_v_list = Ref{Ptr{FancoVcoeff}}()     # cglobal(dlsym(ANCOLIB[], :anco_v_list), FancoVcoeff)
    const mct_list    = Ref{Ptr{Fmctcoefficient}}() # cglobal(dlsym(ANCOLIB[], :mct_list), Fmctcoefficient)

    const t_vector   = Ref{Vector{FancoTcoeff}}()     # unsafe_wrap(Vector{FancoTcoeff}, anco_t_list, 1000)
    const v_vector   = Ref{Vector{FancoVcoeff}}()     # unsafe_wrap(Vector{FancoVcoeff}, anco_v_list, 2000)
    const mct_vector = Ref{Vector{Fmctcoefficient}}() # unsafe_wrap(Vector{Fmctcoefficient}, mct_list, NUMBER_OF_MCT_COEFFICIENTS_MAX)

    function __init__()
        ANCOLIB[] = dlopen(ANCOPATH)

        no_anco_t_list[]             = cglobal(dlsym(ANCOLIB[], :no_anco_t_list), Int32)
        no_anco_v_list[]             = cglobal(dlsym(ANCOLIB[], :no_anco_v_list), Int32)
        number_of_mct_coefficients[] = cglobal(dlsym(ANCOLIB[], :number_of_mct_coefficients), Int32)

        anco_t_list[] = cglobal(dlsym(ANCOLIB[], :anco_t_list), FancoTcoeff)
        anco_v_list[] = cglobal(dlsym(ANCOLIB[], :anco_v_list), FancoVcoeff)
        mct_list[]    = cglobal(dlsym(ANCOLIB[], :mct_list), Fmctcoefficient)

        t_vector[]   = unsafe_wrap(Vector{FancoTcoeff}, anco_t_list[], 1000)
        v_vector[]   = unsafe_wrap(Vector{FancoVcoeff}, anco_v_list[], 2000)
        mct_vector[] = unsafe_wrap(Vector{Fmctcoefficient}, mct_list[], NUMBER_OF_MCT_COEFFICIENTS_MAX)
    end


    """
    `AngularCoefficientsRatip2013.number_of_t_coeffs()`  
        ... to read the number of calculated T coefficients from the external library.
    """
    number_of_t_coeffs() = unsafe_load(no_anco_t_list[])


    """
    `AngularCoefficientsRatip2013.number_of_v_coeffs()`  
        ... to read the number of calculated V coefficients from the external library.
    """
    number_of_v_coeffs() = unsafe_load(no_anco_v_list[])


    """
    `AngularCoefficientsRatip2013.number_of_mct_coeffs()`  
        ... to read the number of calculated MCT coefficients from the external library.
    """
    number_of_mct_coeffs() = unsafe_load(number_of_mct_coefficients[])


    """
    `AngularCoefficientsRatip2013.view_t_coeffs()`  
        ... to return the array of T coefficients as a view into the memory of the external library.
    """
    view_t_coeffs() = view(t_vector[], 1:number_of_t_coeffs())


    """
    `AngularCoefficientsRatip2013.view_v_coeffs()`  
        ... to return the array of V coefficients as a view into the memory of the external library.
    """
    view_v_coeffs() = view(v_vector[], 1:number_of_v_coeffs())


    """
    `AngularCoefficientsRatip2013.view_mct_coeffs()`  
        ... to return the array of MCT coefficients as a view into the memory of the external library.
    """
    view_mct_coeffs() = view(mct_vector[], 1:number_of_mct_coeffs())


    """
    `AngularCoefficientsRatip2013.angular_coefficients_pair(r::Integer, s::Integer)`  
        ... to calculate the angular coefficients between the `r`-th and `s`-th configuration state function of the 
            CSF set that was loaded into the external library.
    """
    function angular_coefficients_pair(r::Integer, s::Integer)
        ccall(dlsym(ANCOLIB[], :jac_anco_calculate_csf_pair), Nothing, (Ref{Int32}, Ref{Int32}), r, s)
        return AngularTcoeff.(view_t_coeffs()), AngularVcoeff.(view_v_coeffs())
    end


    """
    `AngularCoefficientsRatip2013.angular_coefficients_pair_1p(ν::Integer, r::Integer, s::Integer)`  
        ... to calculate the angular coefficients between the `r`-th and `s`-th configuration state function of the CSF 
            set that was loaded into the external library, for a non-scalar one-particle operator of rank `ν`.
    """
    function angular_coefficients_pair_1p(ν::Integer, r::Integer, s::Integer)
        ccall(dlsym(ANCOLIB[], :jac_anco_calculate_csf_pair_1p), Nothing, (Ref{Int32}, Ref{Int32}, Ref{Int32}), ν, r, s)
        return AngularTcoeff.(view_t_coeffs())
    end


    """
    `AngularCoefficientsRatip2013.load_csl(nocsf, nwshells, nwcore, number_of_electrons, subshell, 
                                           totalJs, parities, occupations, seniorities, subshellJs, subshellXs)` 
        ... a wrapper function for the call to the subroutine `jac_anco_load_csl` of the external library.
    """
    load_csl(nocsf::Integer, nwshells::Integer, nwcore::Integer, number_of_electrons::Integer,
             subshell::AbstractVector{Fnkappa},
             totalJs::AbstractVector{Int8}, parities::AbstractVector{UInt8},
             occupations::AbstractMatrix{Int8}, seniorities::AbstractMatrix{Int8},
             subshellJs::AbstractMatrix{Int8},  subshellXs::AbstractMatrix{Int8}) = 
        ccall(dlsym(ANCOLIB[], :jac_anco_load_csl), Nothing,
             (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Fnkappa}, Ref{Int8},
              Ref{UInt8}, Ref{Int8}, Ref{Int8}, Ref{Int8}, Ref{Int8}),
              nocsf, nwshells, nwcore, number_of_electrons, subshell, totalJs, parities,
              occupations, seniorities, subshellJs, subshellXs)


    """
    `AngularCoefficientsRatip2013.load_csl(csf::CsfR)`  
        ... to load a single relativistic configuration state function into the external library. Returns the list of
            subshells the CSF is based on.
    """
    function load_csl(csf::CsfR)
        number_of_electrons = sum(csf.occupation)
        subshells = getsubshells(csf)
        nkappas = Fnkappa.(subshells)
        nwshells = length(nkappas)
        occ = reshape(Int8.(csf.occupation), nwshells, 1)
        sen = reshape(Int8.(csf.seniority), nwshells, 1)
        suJ = reshape(Int8.(twice.(csf.subshellJ)), nwshells, 1)
        suX = reshape(Int8.(twice.(csf.subshellX)), nwshells, 1)
        totalJs =  Int8[twice(csf.J)]
        # parities = UInt8[Char(csf.parity)]
        if   csf.parity == Basics.plus    parities = UInt8[Char('+')]   else    parities = UInt8[Char('-')]   end
        nwcore = 0
        for i in eachindex(nkappas)
            if  occ[i]!=maxocc(nkappas[i])
                nwcore = i-1; break
            end
        end
        load_csl(1, nwshells, nwcore, number_of_electrons, nkappas, totalJs, parities, occ, sen, suJ, suX)
        subshells
    end


    """
    `AngularCoefficientsRatip2013.load_csl(csfa::CsfR, csfb::CsfR)`  
        ... to load two relativistic configuration state functions into the external library. Returns a common
            list of subshells for both CSFs.
    """
    function load_csl(csfa::CsfR, csfb::CsfR)
        number_of_electrons = sum(csfa.occupation)
        sum(csfb.occupation) == number_of_electrons || error("CSFs must have the same number of electrons")
        if have_same_subshells(csfa, csfb)
            subshells = getsubshells(csfa)
            nwshells = length(subshells)
            occ = Matrix{Int8}(hcat(csfa.occupation, csfb.occupation))
            sen = Matrix{Int8}(hcat(csfa.seniority, csfb.seniority))
            suJ = Matrix{Int8}(twice.(hcat(csfa.subshellJ, csfb.subshellJ)))
            suX = Matrix{Int8}(twice.(hcat(csfa.subshellX, csfb.subshellX)))
        else
            subshells = sort(unique(vcat(getsubshells(csfa), getsubshells(csfb))))
            nwshells = length(subshells)
            occ = zeros(Int8, nwshells, 2)
            sen = zeros(Int8, nwshells, 2)
            suJ = zeros(Int8, nwshells, 2)
            suX = zeros(Int8, nwshells, 2)
            aindex, bindex = 1, 1
            ashells, bshells = getsubshells.((csfa, csfb))
            for i in eachindex(subshells)
                if ashells[aindex] == subshells[i]
                    occ[i,1] = csfa.occupation[aindex]
                    sen[i,1] = csfa.seniority[aindex]
                    suJ[i,1] = twice(csfa.subshellJ[aindex])
                    suX[i,1] = twice(csfa.subshellX[aindex])
                    endof(ashells)<aindex && (aindex+=1)
                else
                    suX[i,1] = i==1 ? 0 : suX[i-1,1]
                end
                if bshells[bindex] == subshells[i]
                    occ[i,2] = csfb.occupation[bindex]
                    sen[i,2] = csfb.seniority[bindex]
                    suJ[i,2] = twice(csfb.subshellJ[bindex])
                    suX[i,2] = twice(csfb.subshellX[bindex])
                    endof(bshells)<bindex && (bindex+=1)
                else
                    suX[i,2] = i==1 ? 0 : suX[i-1,2]
                end
            end
            aindex==endof(ashells) && bindex==endof(bshells) || error("Currently, only standard subshell ordering is supported")
        end
        totalJs =  Int8[twice(csfa.J), twice(csfb.J)]
        ## parities = UInt8[Char(csfa.parity), Char(csfb.parity)]
        parities = UInt8['a', 'a']
        if csfa.parity == Basics.plus      parities[1] = '+'     else       parities[1] = '-'      end
        if csfb.parity == Basics.plus      parities[2] = '+'     else       parities[2] = '-'      end
        nkappas = Fnkappa.(subshells)
        nwcore = 0
        for i in eachindex(subshells)
            coreocc = maxocc(nkappas[i])
            if  occ[i,1]!=coreocc || occ[i,2]!=coreocc
                nwcore = i-1; break
            end
        end
        load_csl(2, nwshells, nwcore, number_of_electrons, nkappas, totalJs, parities, occ, sen, suJ, suX)
        subshells
    end


    """
    `AngularCoefficientsRatip2013.mct_generate_coefficients(rl::Integer, ru::Integer, sl::Integer, su::Integer, iopar::Integer, rank::Integer)`
        ... to calculate the MCT coefficients between the `rl`..`ru`-th and `sl`..`su`-th configuration state functions of the CSF set
            that was loaded into the external library, for a non-scalar one-particle operator of parity `iopar` and rank `rank`.
    """
    function mct_generate_coefficients(rl::Integer, ru::Integer, sl::Integer, su::Integer, iopar::Integer, rank::Integer)
        ccall(dlsym(ANCOLIB[], :jac_mct_generate_coefficients), Nothing,
              (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}), rl, ru, sl, su, iopar, rank)
        return view_mct_coeffs()
    end

end # module


