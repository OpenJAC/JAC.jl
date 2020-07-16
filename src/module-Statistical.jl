
"""
`module  JAC.Statistical`  
    ... a submodel of JAC that contains all structs and methods to deal with time-dependent statistical tensors.
"""
module Statistical

    using ..ManyElectron


    """
    `struct  Statistical.ResonanceR`  ... defines a type for a resonance state in the continuum with a well-defined bound-ionic core, 
             one or several electrons in the continuum, a widths as well as a loss rate due to additional decay processes.

        + isBound            ::Bool                ... True if it just refers to a bound state with no electron in the continuum.
        + ionLevel           ::Level               ... Level of the bound-state core.
        + widths             ::Float64             ... Widths of the resonance state.
        + lossRate           ::Float64             ... Loss rate due to processes that are not considered explicitly in the time evolution.
    """
    struct ResonanceR
        isBound              ::Bool   
        ionLevel             ::Level
        widths               ::Float64
        lossRate             ::Float64
    end 


    """
    `Statistical.ResonanceR()`  ... constructor for an `empty` instance of Statistical.ResonanceR().
    """
    function ResonanceR()
        ResonanceR( true, Level(), 0., 0. )
    end


    # `Base.show(io::IO, resonance::Statistical.ResonanceR)`  ... prepares a proper printout of the variable resonance::Statistical.ResonanceR.
    function Base.show(io::IO, resonance::Statistical.ResonanceR) 
        println(io, "isBound:              $(resonance.isBound)  ")
        println(io, "ionLevel:             $(resonance.ionLevel)  ")
        println(io, "widths:               $(resonance.widths)  ")
        println(io, "lossRate:             $(resonance.lossRate)  ")
    end


    """
    `struct  Statistical.Tensor`  ... defines a type for a statistical tensor of given rank k, projection q and w.r.t. two resonances.

        + k                  ::AngularJ64          ... Rank of the tensor
        + q                  ::AngularM64          ... Projection of the tensor rank.
        + a                  ::ResonanceR          ... Level/resonance a with total symmetry alpha J.
        + b                  ::ResonanceR          ... Level/resonance a with total symmetry alpha' J'.
        + value              ::Complex{Float64}    ... value of the statistical tensor
    """
    struct Tensor
        k                    ::AngularJ64
        q                    ::AngularM64
        a                    ::ResonanceR
        b                    ::ResonanceR
        value                ::Complex{Float64}
    end 


    """
    `Statistical.Tensor()`  ... constructor for an `empty` instance of Statistical.Tensor().
    """
    function Tensor()
        Tensor( AngularJ64(0), AngularM64(0), ResonanceR(), ResonanceR(), Complex(0.) )
    end


    # `Base.show(io::IO, tensor::Statistical.Tensor)`  ... prepares a proper printout of the variable tensor::Statistical.Tensor.
    function Base.show(io::IO, tensor::Statistical.Tensor) 
        println(io, "k:              $(tensor.k)  ")
        println(io, "q:              $(tensor.q)  ")
        println(io, "a:              $(tensor.a)  ")
        println(io, "b:              $(tensor.b)  ")
        println(io, "value:          $(tensor.value)  ")
    end

end # module
