
"""
`module  JAC.HighHarmonic
    ... a submodel of JAC that contains all methods to set-up and process high-harmonic computations.
"""
module HighHarmonic

    using Printf, ..Basics, ..Defaults, ..Radial, ..ManyElectron, ..Nuclear, ..TableStrings

    """
    `abstract type HighHarmonic.AbstractApproachxx` 
        ... defines an abstract and a number of singleton types for the computational approach/model/scenario that 
            is considered for generating spectra of high harmonics.

      + struct AverageSCA         
        ... to evaluate the level structure and transitions of all involved levels in single-configuration approach but 
            without configuration interaction and from just a single orbital set.
            
      + struct SCA                
        ... to evaluate the level structure and transitions of all involved levels in single-configuration approach but 
            by calculating all fine-structure resolved transitions.
    """
    abstract type  AbstractApproach                   end
    struct         AverageSCA  <:  AbstractApproach   end
    struct         SCA         <:  AbstractApproach   end
    


    """
    `struct  HighHarmonic.Block`  
        ... defines a type for an individual block of configurations that are treatet together within the cascade. Such an block is given 
            by a list of configurations that may occur as initial- or final-state configurations in some step of the canscade and that are 
            treated together in a single multiplet to allow for configuration interaction but to avoid 'double counting' of individual levels.

        + NoElectrons     ::Int64                     ... Number of electrons in this block.
        + confs           ::Array{Configuration,1}    ... List of one or several configurations that define the multiplet.
        + hasMultiplet    ::Bool                      ... true if the (level representation in the) multiplet has already been computed and
                                                          false otherwise.
        + multiplet       ::Multiplet                 ... Multiplet of the this block.
    """
    struct  Block
        NoElectrons       ::Int64
        confs             ::Array{Configuration,1} 
        hasMultiplet      ::Bool
        multiplet         ::Multiplet  
    end 


    """
    `HighHarmonic.Block()`  ... constructor for an 'empty' instance of a HighHarmonic.Block.
    """
    function Block()
        Block( )
    end


    # `Base.show(io::IO, block::HighHarmonic.Block)`  ... prepares a proper printout of the variable block::HighHarmonic.Block.
    function Base.show(io::IO, block::HighHarmonic.Block) 
        println(io, "NoElectrons:        $(block.NoElectrons)  ")
    end


    """
    `struct  HighHarmonic.Step`  
        ... defines a type for an individual step of an 

        + process          ::JBasics.AtomicProcess         ... Atomic process that 'acts' in this step of the cascade.
        + settings         ::Union{PhotoEmission.Settings, AutoIonization.Settings}        
                                                       ... Settings for this step of the cascade.
        + initialConfs     ::Array{Configuration,1}    ... List of one or several configurations that define the initial-state multiplet.
    """
    struct  Step
        process            ::Basics.AtomicProcess
    end 


    """
    `struct  HighHarmonic.Computation`  
        ... defines a type for the computation of high-harmonic spectra, ...

        + name               ::String                         ... A name for the cascade
        + nuclearModel       ::Nuclear.Model                  ... Model, charge and parameters of the nucleus.
        + grid               ::Radial.Grid                    ... The radial grid to be used for the computation.
        + asfSettings        ::AsfSettings                    ... Provides the settings for the SCF process.
        + approach           ::Cascade.AbstractApproach       ... Computational approach/model that is applied to generate and evaluate the 
                                                                  cascade; possible approaches are: {'single-configuration', ...}
        + processes          ::Array{Basics.AtomicProcess,1}   ... List of the atomic processes that are supported and should be included into the 
                                                                  cascade.
        + initialConfs       ::Array{Configuration,1}         ... List of one or several configurations that contain the level(s) from which the 
                                                                  cascade starts.
        + initialLevels      ::Array{Tuple{Int64,Float64},1}  ... List of one or several (tupels of) levels together with their relative population 
                                                                  from which the cascade starts.
        + maxElectronLoss    ::Int64                          ... (Maximum) Number of electrons in which the initial- and final-state 
                                                                  configurations can differ from each other; this also determines the maximal steps 
                                                                  of any particular decay path.
        + NoShakeDisplacements ::Int64                        ... Maximum number of electron displacements due to shake-up or shake-down processes 
                                                                  in any individual step of the cascade.
        + shakeFromShells ::Array{Shell,1}                    ... List of shells from which shake transitions may occur.
        + shakeToShells   ::Array{Shell,1}                    ... List of shells into which shake transitions may occur.
        + steps           ::Array{Cascade.Step,1}             ... List of individual steps between well-defined atomic multiplets that are 
                                                                  included into the cascade.
    """
    struct  Computation
        name                 ::String
        nuclearModel         ::Nuclear.Model
    end 


    """
    `Cascade.Computation()`  ... constructor for an 'empty' instance of a Cascade.Computation.
    """
    function Computation()
        Computation("",  Nuclear.Model(0.), Radial.Grid(), AsfSettings(), averageSCA, [Auger], 
                    Configuration[], [(0, 0.)], 0, 0, Shell[], Shell[], CascadeComputationStep[] )
    end


    # `Base.show(io::IO, computation::Cascade.Computation)`  ... prepares a proper printout of the variable computation::Cascade.Computation.
    function Base.show(io::IO, computation::Cascade.Computation) 
        println(io, "name:                     $(computation.name)  ")
    end

    
    """
    `abstract type  HighHarmonic.Property`  
        ... defines an abstract and various singleton types for the different properties that can be obtained from the simulation of cascade data.

        + struct IonDist               ... simulate the 'ion distribution' as it is found after all cascade processes are completed.
        + struct FinalDist             ... simulate the 'final-level distribution' as it is found after all cascade processes are completed.
    """
    abstract type  Property                      end
    struct   IonDist              <:  Property   end
    struct   FinalDist            <:  Property   end


    """
    abstract type  HighHarmonic.SimulationMethod`  
        ... defines a enumeration for the various methods that can be used to run the simulation of cascade data.

        + struct ProbPropagation     ... to propagate the (occupation) probabilites of the levels until no further changes occur.
        + struct MonteCarlo          ... to simulate the cascade decay by a Monte-Carlo approach of possible pathes (not yet considered).
        + struct RateEquations       ... to solve the cascade by a set of rate equations (not yet considered).
    """
    abstract type  SimulationMethod                  end
    struct   ProbPropagation  <:  SimulationMethod   end
    struct   MonteCarlo       <:  SimulationMethod   end
    struct   RateEquations    <:  SimulationMethod   end


    """
    `struct  HighHarmonic.SimulationSettings`  ... defines settings for performing the simulation of some cascade (data).

        + minElectronEnergy ::Float64     ... Minimum electron energy for the simulation of electron spectra.
        + maxElectronEnergy ::Float64     ... Maximum electron energy for the simulation of electron spectra.
        + minPhotonEnergy   ::Float64     ... Minimum photon energy for the simulation of electron spectra.
        + maxPhotonEnergy   ::Float64     ... Maximum photon energy for the simulation of electron spectra..
    """
    struct  SimulationSettings
        minElectronEnergy   ::Float64
        maxElectronEnergy   ::Float64
        minPhotonEnergy     ::Float64
        maxPhotonEnergy     ::Float64
    end 


    """
    `HighHarmonic.SimulationSettings()`  ... constructor for an 'empty' instance of ...
    """
    function SimulationSettings()
        SimulationSettings(0., 1.0e6,  0., 1.0e6 )
    end


    """
    `Cascade.computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                      linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)` 
        ... to compute the flourescence and Auger yields for a single outcome as specified by its level; an 
            outcome::DecayYield.Outcome is returned in which are physical parameters are now specified.
    """
    function computeDecayYieldOutcome(outcome::DecayYield.Outcome, linesR::Array{PhotoEmission.Line,1}, 
                                      linesA::Array{AutoIonization.Line,1}, settings::DecayYield.Settings)
        # Identify the level key of the given level also in the lists of radiative and Auger lines
        level = outcome.level;    levelKey = LevelKey( LevelSymmetry(level.J, level.parity), level.index, level.energy, 0.)
        similarKey = LevelKey();  rateR = 0.;    rateA = 0.;   NoPhotonLines = 0;   NoAugerLines = 0
        ##x println("** levelKey = $levelKey")
        for  line in linesR
            compareKey = LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
            if   Basics.isSimilar(levelKey, compareKey, 1.0e-3)    println("** compareKey = $compareKey");   similarKey = deepcopy(compareKey)    end
        end
        return( newOutcome )
    end

end # module


