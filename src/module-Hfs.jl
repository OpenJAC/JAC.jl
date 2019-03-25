
"""
`module  JAC.Hfs`  ... a submodel of JAC that contains all methods for computing HFS A and B coefficients; it is using JAC, JAC.ManyElectron, 
                       JAC.Radial.
"""
module Hfs

    using Printf, JAC, JAC.ManyElectron, JAC.Radial
    global JAC_counter = 0
    

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
        matrixT2   ::Array
    end 


    """
    `JAC.Hfs.InteractionMatrix()`  ... constructor for an `empty` instance of InteractionMatrix.
    """
    function InteractionMatrix()
        InteractionMatrix(false, false, zeros(2,2), zeros(2,2))
    end


    """
    `Base.show(io::IO, im::Hfs.InteractionMatrix)`  ... prepares a proper printout of the variable InteractionMatrix.
    """
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
    `JAC.Hfs.IJF_Vector()`  ... constructor for an `empty` instance of IJF_Vector.
    """
    function IJF_Vector()
        IJF_Vector(AngularJ64(0), AngularJ64(0), Level())
    end


    """
    `Base.show(io::IO, ijfVector::Hfs.IJF_Vector)`  ... prepares a proper printout of the variable ijfVector.
    """
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
   `JAC.Hfs.IJF_Basis()`  ... constructor for an `empty` instance of IJF_Basis.
    """
    function IJF_Basis()
        IJF_Basis(IJF_Vector[], Basis())
    end


    """
    `Base.show(io::IO, ijfBasis::Hfs.IJF_Basis)`  ... prepares a proper printout of the variable ijfBasis.
    """
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
        + hasStateRep    ::Bool             ... Determines whether this IJF-level 'point' to a physical representation of an state (i.e. an basis 
                                                and corresponding mixing coefficient vector), or just to empty instances, otherwise.
        + basis          ::IJF_Basis        ... basis for this level
        + mc             ::Vector{Float64}  ... Vector of mixing coefficients w.r.t basis.
    """
    struct IJF_Level
        I                ::AngularJ64
        F                ::AngularJ64
        M                ::AngularM64
        parity           ::Parity
        energy           ::Float64
        hasStateRep      ::Bool
        basis            ::IJF_Basis
        mc               ::Vector{Float64}
    end 


    """
    `JAC.Hfs.IJF_Level()`  ... constructor for an `empty` instance of IJF_Level.
    """
    function IJF_Level()
        IJF_Level(AngularJ64(0), AngularJ64(0), AngularJ64(0), AngularM64(0), JAC.plus, 0., false, IJF_Basis(), Float64[])
    end


    """
    `Base.show(io::IO, ijfLevel::Hfs.IJF_Level)`  ... prepares a proper printout of the variable ijfLevel::Hfs.IJF_Level.
    """
    function Base.show(io::IO, ijfLevel::Hfs.IJF_Level) 
        println(io, "I:              $(ijfLevel.I)  ")
        println(io, "J:              $(ijfLevel.J)  ")
        println(io, "F:              $(ijfLevel.F)  ")
        println(io, "M:              $(ijfLevel.M)  ")
        println(io, "parity:         $(ijfLevel.parity)  ")
        println(io, "energy:         $(ijfLevel.energy)  ")
        println(io, "hasStateRep:    $(ijfLevel.hasStateRep)  ")
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
    `JAC.Hfs.IJF_Multiplet()`  ... constructor for an `empty` instance of Hfs.IJF_Multiplet.
    """
    function IJF_Multiplet()
        IJF_Multiplet("", IJF_Level[])
    end


    """
    `Base.show(io::IO, ijfMultiplet::Hfs.IJF_Multiplet)`  ... prepares a proper printout of the variable ijfMultiplet::Hfs.IJF_Multiplet.
    """
    function Base.show(io::IO, ijfMultiplet::Hfs.IJF_Multiplet) 
        println(io, "name:           $(ijfMultiplet.name)  ")
        println(io, "levelFs:        $(ijfMultiplet.levelFs)  ")
    end


    """
    `struct  Hfs.Outcome`  ... defines a type to keep the outcome of a HFS computation, such as the HFS A and B coefficients as well other results.

        + Jlevel                    ::Level            ... Atomic level to which the outcome refers to.
        + AIoverMu                  ::Float64          ... HFS A * I / mu value.
        + BoverQ                    ::Float64          ... HFS B / Q value
        + amplitudeT1               ::Complex{Float64} ... T1 amplitude
        + amplitudeT2               ::Complex{Float64} ... T2 amplitude
        + nuclearI                  ::AngularJ64       ... nuclear spin
        + hasFlevels                ::Bool             ... True if the F-values and F-energies are specified explicitly
        + Flevels                   ::Bool             ... Array to keep F-values and F-energies. ... Not yet properly defined.
    """
    struct Outcome 
        Jlevel                      ::Level 
        AIoverMu                    ::Float64
        BoverQ                      ::Float64
        amplitudeT1                 ::Complex{Float64}
        amplitudeT2                 ::Complex{Float64}
        nuclearI                    ::AngularJ64
        hasFlevels                  ::Bool
        Flevels                     ::Bool
    end 


    """
    `JAC.Hfs.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of HFS properties.
    """
    function Outcome()
        Outcome(Level(), 0., 0., 0., 0., AngularJ64(0), false, false)
    end


    """
    `Base.show(io::IO, outcome::Hfs.Outcome)`  ... prepares a proper printout of the variable outcome::Hfs.Outcome.
    """
    function Base.show(io::IO, outcome::Hfs.Outcome) 
        println(io, "Jlevel:                    $(outcome.Jlevel)  ")
        println(io, "AIoverMu:                  $(outcome.AIoverMu)  ")
        println(io, "BoverQ:                    $(outcome.BoverQ)  ")
        println(io, "amplitudeT1:               $(outcome.amplitudeT1)  ")
        println(io, "amplitudeT2:               $(outcome.amplitudeT2)  ")
        println(io, "nuclearI:                  $(outcome.nuclearI)  ")
        println(io, "hasFlevels:                $(outcome.hasFlevels)  ")
        println(io, "Flevels:                   $(outcome.Flevels)  ")
    end


    """
    `struct  Settings`  ... defines a type for the details and parameters of computing HFS A and B coefficients.

        + calcT1                    ::Bool             ... True if T1-amplitudes (HFS A values) need to be calculated, and false otherwise.
        + calcT2                    ::Bool             ... True if T2-amplitudes (HFS B values) need to be calculated, and false otherwise.
        + calcNondiagonal           ::Bool             ... True if also (non-)diagonal hyperfine amplitudes are to be calculated and
                                                           printed, and false otherwise.
        + calcIJFexpansion          ::Bool             ... True if the selected atomic levels are to be represented in a IJF-coupled basis,
                                                           and false otherwise.
        + printBeforeComputation    ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
        + printDeltaEF              ::Bool             ... True if the energy shift of E_F (with regard to E_J) is to be printed, and false otherwise.
        + selectLevels              ::Bool             ... True if individual levels are selected for the computation.
        + selectedLevels            ::Array{Int64,1}   ... List of selected levels.
    """
    struct Settings
        calcT1                      ::Bool
        calcT2                      ::Bool
        calcNondiagonal             ::Bool 
        calcIJFexpansion            ::Bool 
        printBeforeComputation      ::Bool 
        printDeltaEF                ::Bool
        selectLevels                ::Bool
        selectedLevels              ::Array{Int64,1}
    end 


    """
    `JAC.Hfs.Settings(; calcT1::Bool=true,` calcT2::Bool=false, calcNondiagonal::Bool=false, calcIJFexpansion::Bool=false, 
                        printBeforeComputation::Bool=false, printDeltaEF::Bool=false, 
                        selectLevels::Bool=false, selectedLevels::Array{Int64,1}=Int64[]) 
         ... keyword constructor to overwrite selected value of Einstein line computations.
    """
    function Settings(; calcT1::Bool=true, calcT2::Bool=false, calcNondiagonal::Bool=false, calcIJFexpansion::Bool=false, 
                        printBeforeComputation::Bool=false, printDeltaEF::Bool=false, 
                        selectLevels::Bool=false, selectedLevels::Array{Int64,1}=Int64[])
        Settings(calcT1, calcT2, calcNondiagonal, calcIJFexpansion, printBeforeComputation, printDeltaEF, selectLevels, selectedLevels)
    end


    """
    `Base.show(io::IO, settings::Hfs.Settings)`  ... prepares a proper printout of the variable settings::Hfs.Settings.
    """
    function Base.show(io::IO, settings::Hfs.Settings) 
        println(io, "calcT1:                   $(settings.calcT1)  ")
        println(io, "calcT2:                   $(settings.calcT2)  ")
        println(io, "calcNondiagonal:          $(settings.calcNondiagonal)  ")
        println(io, "calcIJFexpansion:         $(settings.calcIJFexpansion)  ")
        println(io, "printBeforeComputation:   $(settings.printBeforeComputation)  ")
        println(io, "printDeltaEF:             $(settings.printDeltaEF)  ")
        println(io, "selectLevels:             $(settings.selectLevels)  ")
        println(io, "selectedLevels:           $(settings.selectedLevels)  ")
    end


    """
    `JAC.Hfs.amplitude(kind::String, rLevel::Level, sLevel::Level, grid::Radial.Grid; printout::Bool=true)` ... to compute either 
         the  T^(1) or T^(2) hyperfine amplitude <alpha_r J_r || T^(n)) || alpha_s J_s>  for a given pair of levels. 
         A value::ComplexF64 is returned.
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
                    wa = compute("angular coefficients: 1-p, Grasp92", 0, 1, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.hfs_t1(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                        #
                        ## println("**  <$(coeff.a) || t1 || $(coeff.b)>  = $(coeff.T * tamp)   = $(coeff.T) * $tamp" )
                        me = me + coeff.T * tamp  
                    end
                #
                elseif  kind == "T^(2) amplitude"
                #--------------------------------
                    wa = compute("angular coefficients: 1-p, Grasp92", 0, 2, rLevel.basis.csfs[r], sLevel.basis.csfs[s])
                    for  coeff in wa
                        ja = JAC.subshell_2j(rLevel.basis.orbitals[coeff.a].subshell)
                        jb = JAC.subshell_2j(sLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.hfs_t2(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                        #
                        ## println("**  <$(coeff.a) || t2 || $(coeff.b)>  = $(coeff.T * tamp)   = $(coeff.T) * $tamp" )
                        me = me + coeff.T * tamp  
                    end
                #
                else    error("stop a")
                end
                #
                matrix[r,s] = me
                println("**                                                                    " *
                        "                      <$r || $kind || $s> = $me" )
            end
        end
        if  printout   printstyled("done.\n", color=:light_green)   end
        amplitude = transpose(rLevel.mc) * matrix * sLevel.mc 
        #
        return( amplitude )
    end


    """
    `JAC.Hfs.computeAmplitudesProperties(outcome::Hfs.Outcome, grid::Radial.Grid, settings::Hfs.Settings, im::Hfs.InteractionMatrix) 
        ... to compute all amplitudes and properties of for a given level; an outcome::Hfs.Outcome is returned for which the 
            amplitudes and properties are now evaluated explicitly.
    """
    function  computeAmplitudesProperties(outcome::Hfs.Outcome, grid::Radial.Grid, settings::Hfs.Settings, im::Hfs.InteractionMatrix)
        AIoverMu = 0.0;   BoverQ = 0.0;   amplitudeT1 = 0.0;   amplitudeT2 = 0.0;    J = JAC.AngularMomentum.oneJ(outcome.Jlevel.J)
        ##x println("++ AngJ = $(outcome.Jlevel.J)  J = $J ")
        if  settings.calcT1  &&  outcome.Jlevel.J != AngularJ64(0)
            if  im.calcT1   amplitudeT1 = transpose(outcome.Jlevel.mc) * im.matrixT1 * outcome.Jlevel.mc
            else            amplitudeT1 = JAC.Hfs.amplitude("T^(1) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
            end
            
            AIoverMu = amplitudeT1 / sqrt(J * (J+1))
        end
        #
        if  settings.calcT2  &&  outcome.Jlevel.J != AngularJ64(0)
            if  im.calcT2   amplitudeT2 = transpose(outcome.Jlevel.mc) * im.matrixT2 * outcome.Jlevel.mc
            else            amplitudeT2 = JAC.Hfs.amplitude("T^(2) amplitude", outcome.Jlevel, outcome.Jlevel, grid)
            end
            
            BoverQ   = 2 * amplitudeT2 * sqrt( (2J-1) / ((J+1)*(2J+3)) )  ## * sqrt(J)
        end
        newOutcome = Hfs.Outcome( outcome.Jlevel, AIoverMu, BoverQ, amplitudeT1, amplitudeT2, 
                                  outcome.nuclearI, outcome.hasFlevels, outcome.Flevels)
        return( newOutcome )
    end


    """
    `JAC.Hfs.computeHyperfineMultiplet(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)`  
         ... to compute a hyperfine multiplet, i.e. a representation of hyperfine levels within a hyperfine-coupled basis as defined by the
         given (electronic) multiplet; a hfsMultiplet::IJF_Multiplet is returned.
    """
    function computeHyperfineMultiplet(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)
        println("")
        printstyled("JAC.Hfs.computeHyperfineMultiplet(): The computation of the hyperfine multiplet starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------- \n", color=:light_green)
        #
        hfsBasis     = JAC.Hfs.defineHyperfineBasis(multiplet, nm, settings)
        hfsMultiplet = JAC.Hfs.computeHyperfineRepresentation(hfsBasis, nm, grid, settings)
        # Print all results to screen
        hfsMultiplet = JAC.sort("multiplet: by energy", hfsMultiplet)
        JAC.tabulate("multiplet: energies",                                      hfsMultiplet) 
        JAC.tabulate("multiplet: energy of each level relative to lowest level", hfsMultiplet) 
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary     
            JAC.tabulate("multiplet: energies",                                      hfsMultiplet, stream=iostream) 
            JAC.tabulate("multiplet: energy of each level relative to lowest level", hfsMultiplet, stream=iostream) 
        end
        #
        if    output    return( hfsMultiplet )
        else            return( nothing )
        end
    end


    """
    `JAC.Hfs.computeHyperfineRepresentation(hfsBasis::IJF_Basis, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings)`  
         ... to set-up and diagonalized the Hamiltonian matrix of H^(DFB) + H^(hfs) within the atomic hyperfine (IJF-coupled) basis;
         a hfsMultiplet::IJF_Multiplet is returned.
    """
    function computeHyperfineRepresentation(hfsBasis::IJF_Basis, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings)
        n = length(hfsBasis.vectors);   matrix = zeros(n,n)
        for  r = 1:n
            for  s = 1:n
                matrix[r,s] = 0.
                if  r == s    matrix[r,s] = matrix[r,s] + hfsBasis.vectors[r].levelJ.energy                  end
                if  hfsBasis.vectors[r].F              !=  hfsBasis.vectors[s].F                continue     end
                if  hfsBasis.vectors[r].levelJ.parity  !=  hfsBasis.vectors[s].levelJ.parity    continue     end
                # Now add the hyperfine matrix element for K = 1 (magnetic) and K = 2 (electric) contributions
                spinI = nm.spinI;   Jr = hfsBasis.vectors[r].levelJ.J;   Js = hfsBasis.vectors[s].levelJ.J;   F = hfsBasis.vectors[r].F
                Ix    = JAC.AngularMomentum.oneJ(spinI)
                #
                wa = JAC.AngularMomentum.phaseFactor([spinI, +1, Jr, +1, F])  * 
                     JAC.AngularMomentum.Wigner_6j(spinI, Jr, F, Js, spinI, AngularJ64(1))
                wb = JAC.Hfs.amplitude("T^(1) amplitude", hfsBasis.vectors[r].levelJ, hfsBasis.vectors[s].levelJ, grid::Radial.Grid)
                wc = nm.mu * sqrt( (Ix+1)/Ix )
                matrix[r,s] = matrix[r,s] + wa * wb * wc     * 1.0e-6 ## fudge-factor to keep HFS interaction small
                #
                if  spinI  in [AngularJ64(0), AngularJ64(1//2)]                                 continue     end
                wa = JAC.AngularMomentum.phaseFactor([spinI, +1, Jr, +1, F]) * 
                     JAC.AngularMomentum.Wigner_6j(spinI, Jr, F, Js, spinI, AngularJ64(2))
                wb = JAC.Hfs.amplitude("T^(2) amplitude", hfsBasis.vectors[r].levelJ, hfsBasis.vectors[s].levelJ, grid::Radial.Grid)
                wc = nm.Q / 2. * sqrt( (Ix+1)*(2Ix+3)/ (Ix*(2Ix-1)) )
                matrix[r,s] = matrix[r,s] + wa * wb * wc     * 1.0e-6 ## fudge-factor to keep HFS interaction small
            end
        end
        #
        # Diagonalize the matrix and set-up the representation
        eigen    = JAC.diagonalize("matrix: Julia, eigfact", matrix)
        levelFs  = JAC.Hfs.IJF_Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evector   = eigen.vectors[ev];    en = eigen.values[ev]
            parity    = JAC.plus;    F = AngularJ64(0);     MF = AngularM64(0)
            ##x println("len = $(length(hfsBasis.vectors))  evector = $evector   ")
            for  r = 1:length(hfsBasis.vectors)
                if  abs(evector[r]) > 1.0e-6    
                    parity = hfsBasis.vectors[r].levelJ.parity
                    F      = hfsBasis.vectors[r].F;     MF = AngularM64(hfsBasis.vectors[r].F);   break    end
            end
            newlevelF = JAC.Hfs.IJF_Level(nm.spinI, F, MF, parity, en, true, hfsBasis, evector) 
            push!( levelFs, newlevelF)
        end
        hfsMultiplet = JAC.Hfs.IJF_Multiplet("hyperfine", levelFs)
        
        return( hfsMultiplet )
    end


    """
    `JAC.Hfs.computeInteractionMatrix(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)` ... to compute the
            T^1 and/or T^2 interaction matrices for the given basis, i.e. (<csf_r || T^(n)) || csf_s>).
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
                    wa = compute("angular coefficients: 1-p, Grasp92", 0, 1, basis.csfs[r], basis.csfs[s])
                    for  coeff in wa
                        ja   = JAC.subshell_2j(basis.orbitals[coeff.a].subshell)
                        jb   = JAC.subshell_2j(basis.orbitals[coeff.b].subshell)
                        tamp = JAC.InteractionStrength.hfs_t1(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
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
                    wa = compute("angular coefficients: 1-p, Grasp92", 0, 2, basis.csfs[r], basis.csfs[s])
                    for  coeff in wa
                        ja   = JAC.subshell_2j(basis.orbitals[coeff.a].subshell)
                        jb   = JAC.subshell_2j(basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.hfs_t2(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                        matrixT2[r,s] = matrixT2[r,s] + coeff.T * tamp  
                    end
                end
            end
        else   
            calcT2 = false;    matrixT2 = zeros(2,2)
        end
        #
        im = Hfs.InteractionMatrix(calcT1, calcT2, matrixT1, matrixT2)
        println("im = $im")
        #
        return( im )
    end
    


    """
    `JAC.Hfs.computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)`  
         ... to compute (as selected) the HFS A and B parameters as well as hyperfine energy splittings for the levels of the given 
         multiplet and as specified by the given settings. The results are printed in neat tables to screen and, if requested, 
         an arrays{Hfs.Outcome,1} with all the results are returned.
    """
    function computeOutcomes(multiplet::Multiplet, nm::JAC.Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)
        println("")
        printstyled("JAC.Hfs.computeOutcomes(): The computation of the Hyperfine amplitudes and parameters starts now ... \n", color=:light_green)
        printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        outcomes = JAC.Hfs.determineOutcomes(multiplet, settings)
        # Display all selected levels before the computations start
        if  settings.printBeforeComputation    JAC.Hfs.displayOutcomes(outcomes)    end
        # Calculate all amplitudes and requested properties
        im = computeInteractionMatrix(multiplet.levels[1].basis, grid, settings)
        newOutcomes = Hfs.Outcome[]
        for  outcome in outcomes
            newOutcome = JAC.Hfs.computeAmplitudesProperties(outcome, grid, settings, im) 
            push!( newOutcomes, newOutcome)
        end
        # Print all results to screen
        JAC.Hfs.displayResults(stdout, newOutcomes, nm, settings)
        # Compute and display the non-diagonal hyperfine amplitudes, if requested
        if  settings.calcNondiagonal    JAC.Hfs.displayNondiagonal(stdout, multiplet, grid, settings)   end
        printSummary, iostream = JAC.give("summary flag/stream")
        if  printSummary    
            JAC.Hfs.displayResults(iostream, newOutcomes, nm, settings) 
            if  settings.calcNondiagonal    JAC.Hfs.displayNondiagonal(iostream, multiplet, grid, settings)   end
        end
        #
        if    output    return( newOutcomes )
        else            return( nothing )
        end
    end



    """
    `JAC.Hfs.computeMatrixT1_obsolete(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)`  ... to compute the matrix T1 = (<csf_r|| T^(1) ||csf_s>) 
         of the hyperfine interaction for the given CSF basis; a (length(basis.csfs) x length(basis.csfs)}-dimensional matrix::Array{Float64,2) 
         is returned. In the HFS theory, the T1 matrix refers to the magnetic-dipole interaction and occurs in the computation hyperfine A 
         parameters.  
    """
    function  computeMatrixT1_obsolete(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings) 
        n = length(basis.csfs)
  
        print("Compute T1 matrix of dimension $n x $n for the given basis ...")
        matrix = zeros(Float64, n, n)
        for  r = 1:n
            for  s = 1:n
                wa = compute("angular coefficients: 1-p, Ratip2013", 1, basis.csfs[r], basis.csfs[s])
                me = 0.
                for  coeff in wa
                    jj = JAC.subshell_2j(basis.orbitals[coeff.a].subshell)
                    me = me + coeff.T * sqrt( jj + 1) * 
                                        JAC.InteractionStrength.hfs_t1(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                end
                matrix[r,s] = me
            end
        end 
        println(" done.")
        return( matrix )
    end


    """
    `JAC.Hfs.computeMatrixT2_obsolete(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)`  ... to compute the matrix T2 = (<csf_r|| T^(2) ||csf_s>) 
         of the hyperfine interaction for the given CSF basis; a (length(basis.csfs) x length(basis.csfs)}-dimensional matrix::Array{Float64,2) 
         is returned. In the HFS theory, the T2 matrix refers to the electric-quadrupole interaction and occurs in the computation hyperfine A 
         parameters. **Not yet implemented !**
    """
    function  computeMatrixT2_obsolete(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)   
        n = length(basis.csfs)
  
        print("Compute T2 matrix of dimension $n x $n for the given basis ...")
        matrix = zeros(Float64, n, n)
        for  r = 1:n
            for  s = 1:n
                wa = compute("angular coefficients: 1-p, Ratip2013", 2, basis.csfs[r], basis.csfs[s])
                me = 0.
                for  coeff in wa
                    jj = JAC.subshell_2j(basis.orbitals[coeff.a].subshell)
                    me = me + coeff.T * sqrt( jj + 1) * 
                                        JAC.InteractionStrength.hfs_t2(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                end
                matrix[r,s] = me
            end
        end 
        println("   ... done.")
        return( matrix )
    end


    """
    `JAC.Hfs.defineHyperfineBasis(multiplet::Multiplet, nm::JAC.Nuclear.Model, settings::Hfs.Settings)`  ... to define/set-up
         an atomic hyperfine (IJF-coupled) basis for the given electronic multipet; a hfsBasis::IJF_Basis is returned.
    """
    function  defineHyperfineBasis(multiplet::Multiplet, nm::JAC.Nuclear.Model, settings::Hfs.Settings)
        vectors = IJF_Vector[]
        #
        for i = 1:length(multiplet.levels)
            Flist = JAC.oplus(nm.spinI, multiplet.levels[i].J)
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
            sa = sa * JAC.TableStrings.ijfVector(k, hfsBasis.vectors[k])
        end
        if   length(sa) > 5    println(sa)   end
        println(" ")
        
        return( hfsBasis )
    end


    """
    `JAC.Hfs.determineOutcomes(multiplet::Multiplet, settings::Hfs.Settings)`  ... to determine a list of Outcomes's for the computation
         of HFS A- and B-parameters for the given multiplet. It takes into account the particular selections and settings. An Array{Hfs.Outcome,1} 
         is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineOutcomes(multiplet::Multiplet, settings::Hfs.Settings) 
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end
    
        outcomes = Hfs.Outcome[]
        for  i = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
            push!( outcomes, Hfs.Outcome(multiplet.levels[i], 0., 0., 0., 0., AngularJ64(0), false, false) )
        end
        return( outcomes )
    end


    """
    `JAC.Hfs.displayNondiagonal(stream::IO, multiplet::Multiplet, grid::Radial.Grid, settings::Hfs.Settings)`  
         ... to compute and display all non-diagonal hyperfine amplitudes for the selected levels.
         A small neat table of all (pairwise) hyperfine amplitudes is printed but nothing is returned otherwise.
    """
    function  displayNondiagonal(stream::IO, multiplet::Multiplet, grid::Radial.Grid, settings::Hfs.Settings)
        if    settings.selectLevels   selectLevels   = true;   selectedLevels = copy(settings.selectedLevels)
        else                          selectLevels   = false
        end
    
        pairs = Tuple{Int64,Int64}[]
        for  f = 1:length(multiplet.levels)
            if  selectLevels  &&  !(haskey(selectedLevels, f))    continue   end
            for  i = f:length(multiplet.levels)
                if  selectLevels  &&  !(haskey(selectedLevels, i))    continue   end
                push!( pairs, (f,i) )
            end
        end
        #
        println(stream, " ")
        println(stream, "  Selected (non-) diagonal hyperfine amplitudes:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(107))
        sa = "  "
        sa = sa * JAC.TableStrings.center(10, "Level_f"; na=2)
        sa = sa * JAC.TableStrings.center(10, "Level_i"; na=2)
        sa = sa * JAC.TableStrings.center(10, "J^P_f";   na=3)
        sa = sa * JAC.TableStrings.center(57, "T1   --   Amplitudes   --   T2"; na=4);              
        sa = sa * JAC.TableStrings.center(10, "J^P_f";   na=4)
        println(stream, sa);    println(stream, "  ", JAC.TableStrings.hLine(107)) 
        #  
        for  (f,i) in pairs
            sa   = "  ";    
            sa   = sa * JAC.TableStrings.center(10, string(f); na=2)
            sa   = sa * JAC.TableStrings.center(10, string(i); na=2)
            symf = LevelSymmetry( multiplet.levels[f].J, multiplet.levels[f].parity)
            symi = LevelSymmetry( multiplet.levels[i].J, multiplet.levels[i].parity)
            sa   = sa * JAC.TableStrings.center(10, string(symf); na=4)
            T1   = JAC.Hfs.amplitude("T^(1) amplitude", multiplet.levels[f], multiplet.levels[i], grid, printout=false)
            T2   = JAC.Hfs.amplitude("T^(2) amplitude", multiplet.levels[f], multiplet.levels[i], grid, printout=false)
            sa   = sa * @sprintf("%.5e %s %.5e", T1.re, "  ", T1.im) * "    "
            sa   = sa * @sprintf("%.5e %s %.5e", T2.re, "  ", T2.im) * "    "
            sa   = sa * JAC.TableStrings.center(10, string(symi); na=4)
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(107))
        #
        return( nothing )
    end


    """
    `JAC.Hfs.displayOutcomes(outcomes::Array{Hfs.Outcome,1})`  ... to display a list of levels that have been selected for the computations
         A small neat table of all selected levels and their energies is printed but nothing is returned otherwise.
    """
    function  displayOutcomes(outcomes::Array{Hfs.Outcome,1})
        println(" ")
        println("  Selected HFS levels:")
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
        println(" ")
        #
        return( nothing )
    end


    """
    `JAC.Hfs.displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1}, nm::JAC.Nuclear.Model, settings::Hfs.Settings)`  
        ... to display the energies, A- and B-values, Delta E_F energy shifts, etc. for the selected levels. All nuclear parameters are taken 
        from the nuclear model. A neat table is printed but nothing is returned otherwise.
    """
    function  displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1}, nm::JAC.Nuclear.Model, settings::Hfs.Settings)
        println(stream, " ")
        println(stream, "  HFS parameters and amplitudes:")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(117))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);              
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(31, "A/mu [mu_N]  --  HFS  --  B/Q [barn]"; na=4);              
        sb = sb * JAC.TableStrings.center(31, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(32, "T1 -- Amplitudes -- T2"    ; na=4);        sb = sb * JAC.TableStrings.hBlank(36)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(117)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            energy = outcome.Jlevel.energy
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", energy))           * "      "
            wa = JAC.convert("energy: from atomic", outcome.AIoverMu) / nm.spinI.num * nm.spinI.den 
            wa = wa / (2 * 1836.15267 ) * JAC.give("alpha")  # take 2m_p and alpha into account
            sa = sa * @sprintf("%.8e", wa)                                                   * "    " 
            wa = JAC.convert("energy: from atomic", outcome.BoverQ) 
            wa = wa / JAC.convert("cross section: from atomic to barn", 1.0)  # take Q [barn] into account
            sa = sa * @sprintf("%.8e", wa)                                                   * "    "
            sa = sa * @sprintf("%.8e %s %.8e", outcome.amplitudeT1.re, "  ", outcome.amplitudeT2.re) * "    "
            println(stream, sa )
        end
        println(stream, "  ", JAC.TableStrings.hLine(117))
        #
        #
        if !settings.printDeltaEF   return( nothing )   end
        #
        # Printout the Delta E_F energy shifts of the hyperfine levels |alpha F> with regard to the (electronic) levels |alpha J>
        println(stream, " ")
        println(stream, "  HFS Delta E_F energy shifts with regard to the (electronic) level energies E_J:")
        println(stream, " ")
        println(stream, "    Nuclear spin I:                        $(nm.spinI) ")
        println(stream, "    Nuclear magnetic-dipole moment    mu = $(nm.mu)    ")
        println(stream, "    Nuclear electric-quadrupole moment Q = $(nm.Q)     ")
        println(stream, " ")
        println(stream, "  ", JAC.TableStrings.hLine(90))
        sa = "  ";   sb = "  "
        sa = sa * JAC.TableStrings.center(10, "Level"; na=2);                             sb = sb * JAC.TableStrings.hBlank(12)
        sa = sa * JAC.TableStrings.center(10, "J^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Energy"; na=4);                            sb = sb * JAC.TableStrings.hBlank(18)
        sa = sa * JAC.TableStrings.center(10, "F^P";   na=4);                             sb = sb * JAC.TableStrings.hBlank(14)
        sa = sa * JAC.TableStrings.center(14, "Delta E_F"; na=4);                         
        sb = sb * JAC.TableStrings.center(14, JAC.TableStrings.inUnits("energy"); na=4)
        sa = sa * JAC.TableStrings.center(14, "C factor"; na=4);                         
        println(stream, sa);    println(stream, sb);    println(stream, "  ", JAC.TableStrings.hLine(90)) 
        #  
        for  outcome in outcomes
            sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * JAC.TableStrings.center(10, JAC.TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * JAC.TableStrings.center(10, string(sym); na=4)
            energy = outcome.Jlevel.energy
            #
            sa = sa * @sprintf("%.8e", JAC.convert("energy: from atomic", energy))           * "    "
            J     = JAC.AngularMomentum.oneJ(outcome.Jlevel.J);   spinI = JAC.AngularMomentum.oneJ(nm.spinI)
            Flist = JAC.oplus(nm.spinI, outcome.Jlevel.J)
            first = true
            for  Fang in Flist
                Fsym    = LevelSymmetry(Fang, outcome.Jlevel.parity)
                F       = JAC.AngularMomentum.oneJ(Fang)
                Cfactor = F*(F+1) - J*(J+1) - spinI*(spinI+1)
                energy  = outcome.AIoverMu * nm.mu / spinI * Cfactor / 2.
                if  abs(outcome.BoverQ) > 1.0e-10
                    energy  = energy +  outcome.BoverQ * nm.Q * 3/4 * (Cfactor*(Cfactor+1) - spinI*(spinI+1)*J*(J+1) ) /
                                        ( 2spinI*(2spinI-1)*J*(2J-1) )
                end
                ##x println("energy-b = $energy  outcome.BoverQ = $(outcome.BoverQ) Q = $(nm.Q)")
                sb = JAC.TableStrings.center(10, string(Fsym); na=2)
                sb = sb * JAC.TableStrings.flushright(16, @sprintf("%.8e", JAC.convert("energy: from atomic", energy))) * "    "
                sb = sb * JAC.TableStrings.flushright(12, @sprintf("%.5e", Cfactor))
                #
                if   first    println(stream,  sa*sb );   first = false
                else          println(stream,  JAC.TableStrings.hBlank( length(sa) ) * sb )
                end
            end
        end
        println(stream, "  ", JAC.TableStrings.hLine(90))
        #
        return( nothing )
    end

end # module


