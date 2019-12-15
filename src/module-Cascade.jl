
"""
`module  JAC.Cascade
    ... a submodel of JAC that contains all methods to set-up and process (simple) atomic cascade computations.
"""
module Cascade

    using Printf, ..AutoIonization, ..Basics, ..Defaults, ..DecayYield, ..Radial, ..ManyElectron, ..Nuclear, ..PhotoEmission, 
                  ..PhotoIonization, ..Semiempirical, ..TableStrings

    """
    `abstract type Cascade.AbstractApproach` 
        ... defines an abstract and a number of singleton types for the computational approach/model that is applied to 
            generate and evaluate a cascade.

      + struct AverageSCA         
        ... to evaluate the level structure and transitions of all involved levels in single-configuration approach but 
            without configuration interaction and from just a single orbital set.
            
      + struct SCA                
        ... to evaluate the level structure and transitions of all involved levels in single-configuration approach but 
            by calculating all fine-structure resolved transitions.
    """
    abstract type  AbstractCascadeApproach                   end
    struct         AverageSCA  <:  AbstractCascadeApproach   end
    struct         SCA         <:  AbstractCascadeApproach   end
    struct         UserMCA     <:  AbstractCascadeApproach   end
    


    """
    `struct  Cascade.Block`  
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
    `Cascade.Block()`  ... constructor for an 'empty' instance of a Cascade.Block.
    """
    function Block()
        Block( )
    end


    # `Base.show(io::IO, block::Cascade.Block)`  ... prepares a proper printout of the variable block::Cascade.Block.
    function Base.show(io::IO, block::Cascade.Block) 
        println(io, "NoElectrons:        $(block.NoElectrons)  ")
        println(io, "confs :             $(block.confs )  ")
        println(io, "hasMultiplet:       $(block.hasMultiplet)  ")
        println(io, "multiplet:          $(block.multiplet)  ")
    end


    """
    `struct  Cascade.Step`  
        ... defines a type for an individual step of an excitation and/or decay cascade. Such an individual step is given by a well-defined 
            process, such as Auger, PhotoEmission, or others and two lists of initial- and final-state configuration that are (each) treated 
            together in a multiplet to allow for configuration interaction but to avoid 'double counting' of individual levels.

        + process          ::JBasics.AtomicProcess         ... Atomic process that 'acts' in this step of the cascade.
        + settings         ::Union{PhotoEmission.Settings, AutoIonization.Settings}        
                                                       ... Settings for this step of the cascade.
        + initialConfigs   ::Array{Configuration,1}    ... List of one or several configurations that define the initial-state multiplet.
        + finalConfigs     ::Array{Configuration,1}    ... List of one or several configurations that define the final-state multiplet.
        + initialMultiplet ::Multiplet                 ... Multiplet of the initial-state levels of this step of the cascade.
        + finalMultiplet   ::Multiplet                 ... Multiplet of the final-state levels of this step of the cascade.
    """
    struct  Step
        process            ::Basics.AtomicProcess
        settings           ::Union{PhotoEmission.Settings, AutoIonization.Settings}
        initialConfigs     ::Array{Configuration,1}
        finalConfigs       ::Array{Configuration,1}
        initialMultiplet   ::Multiplet
        finalMultiplet     ::Multiplet
    end 


    """
    `Cascade.Step()`  ... constructor for an 'empty' instance of a Cascade.Step.
    """
    function Step()
        Step( Basics.NoProcess, PhotoEmission.Settings, Configuration[], Configuration[], Multiplet(), Multiplet())
    end


    # `Base.show(io::IO, step::Cascade.Step)`  ... prepares a proper printout of the variable step::Cascade.Step.
    function Base.show(io::IO, step::Cascade.Step) 
        println(io, "process:                $(step.process)  ")
        println(io, "settings:               $(step.settings)  ")
        println(io, "initialConfigs:         $(step.initialConfigs)  ")
        println(io, "finalConfigs:           $(step.finalConfigs)  ")
        println(io, "initialMultiplet :      $(step.initialMultiplet )  ")
        println(io, "finalMultiplet:         $(step.finalMultiplet)  ")
    end


    """
    `struct  Cascade.Computation`  
        ... defines a type for a cascade computation, i.e. for the computation of a whole excitation and/or decay cascade. The data 
            from this computation can be modified, adapted and refined to the practical needs before the actual computations are 
            carried out. Initially, this struct contains the physical metadata about the cascade to be calculated but gets enlarged 
            in course of the computation to keep also wave functions, level multiplets, etc.

        + name               ::String                          ... A name for the cascade
        + nuclearModel       ::Nuclear.Model                   ... Model, charge and parameters of the nucleus.
        + grid               ::Radial.Grid                     ... The radial grid to be used for the computation.
        + asfSettings        ::AsfSettings                     ... Provides the settings for the SCF process.
        + approach           ::Cascade.AbstractCascadeApproach ... Computational approach/model that is applied to generate and evaluate the 
                                                                   cascade; possible approaches are: {'single-configuration', ...}
        + processes          ::Array{Basics.AtomicProcess,1}   ... List of the atomic processes that are supported and should be included into the 
                                                                   cascade.
        + initialConfs       ::Array{Configuration,1}          ... List of one or several configurations that contain the level(s) from which the 
                                                                   cascade starts.
        + initialLevels      ::Array{Tuple{Int64,Float64},1}   ... List of one or several (tupels of) levels together with their relative population 
                                                                  from which the cascade starts.
        + maxElectronLoss    ::Int64                           ... (Maximum) Number of electrons in which the initial- and final-state 
                                                                   configurations can differ from each other; this also determines the maximal steps 
                                                                   of any particular decay path.
        + NoShakeDisplacements ::Int64                         ... Maximum number of electron displacements due to shake-up or shake-down processes 
                                                                   in any individual step of the cascade.
        + shakeFromShells ::Array{Shell,1}                     ... List of shells from which shake transitions may occur.
        + shakeToShells   ::Array{Shell,1}                     ... List of shells into which shake transitions may occur.
        ## + steps           ::Array{Cascade.Step,1}              ... List of individual steps between well-defined atomic multiplets that are 
        ##                                                           included into the cascade.
    """
    struct  Computation
        name                 ::String
        nuclearModel         ::Nuclear.Model
        grid                 ::Radial.Grid
        asfSettings          ::AsfSettings
        approach             ::Cascade.AbstractCascadeApproach
        processes            ::Array{Basics.AtomicProcess,1}
        initialConfigs       ::Array{Configuration,1}
        initialLevels        ::Array{Tuple{Int64,Float64},1}
        maxElectronLoss      ::Int64
        NoShakeDisplacements ::Int64
        shakeFromShells      ::Array{Shell,1}
        shakeToShells        ::Array{Shell,1}
        ## steps                ::Array{Cascade.Step,1}
    end 


    """
    `Cascade.Computation()`  ... constructor for an 'default' instance of a Cascade.Computation.
    """
    function Computation()
        Computation("Default cascade computation",  Nuclear.Model(10.), Radial.Grid(), AsfSettings(), AverageSCA(), [Radiative], 
                    Configuration[], [(0, 0.)], 0, 0, Shell[], Shell[] )
    end


    """
    `Cascade.Computation(comp::Cascade.Computation;`
        
                name=..,               nuclearModel=..,             grid=..,              asfSettings=..,     
                approach=..,           processes=..,                initialConfigs=..,    initialLevels=..,
                maxElectronLoss=..,    NoShakeDisplacements=..,     shakeFromShells=..,   shakeToShells=.. )
                
        ... constructor for re-defining the computation::Cascade.Computation.
    """
    function Computation(comp::Cascade.Computation;                              
        name::Union{Nothing,String}=nothing,                                  nuclearModel::Union{Nothing,Nuclear.Model}=nothing,
        grid::Union{Nothing,Radial.Grid}=nothing,                             asfSettings::Union{Nothing,AsfSettings}=nothing,    
        approach::Union{Nothing,Cascade.AbstractCascadeApproach}=nothing,     processes::Union{Nothing,Array{Basics.AtomicProcess,1}}=nothing,  
        initialConfigs::Union{Nothing,Array{Configuration,1}}=nothing,        initialLevels::Union{Nothing,Array{Tuple{Int64,Float64},1}}=nothing,  
        maxElectronLoss::Union{Nothing,Int64}=nothing,                        NoShakeDisplacements::Union{Nothing,Int64}=nothing,  
        shakeFromShells::Union{Nothing,Array{Shell,1}}=nothing,               shakeToShells::Union{Nothing,Array{Shell,1}}=nothing ) 
        ##x steps::Union{Nothing,Array{Cascade.Step,1}}=nothing )

        if  name                 == nothing   namex                 = comp.name                    else  namex = name                                   end 
        if  nuclearModel         == nothing   nuclearModelx         = comp.nuclearModel            else  nuclearModelx = nuclearModel                   end 
        if  grid                 == nothing   gridx                 = comp.grid                    else  gridx = grid                                   end 
        if  asfSettings          == nothing   asfSettingsx          = comp.asfSettings             else  asfSettingsx = asfSettings                     end 
        if  approach             == nothing   approachx             = comp.approach                else  approachx = approach                           end 
        if  processes            == nothing   processesx            = comp.processes               else  processesx = processes                         end 
        if  initialConfigs       == nothing   initialConfigsx       = comp.initialConfigs          else  initialConfigsx = initialConfigs               end 
        if  initialLevels        == nothing   initialLevelsx        = comp.initialLevels           else  initialLevelsx = initialLevels                 end 
        if  maxElectronLoss      == nothing   maxElectronLossx      = comp.maxElectronLoss         else  maxElectronLossx = maxElectronLoss             end 
        if  NoShakeDisplacements == nothing   NoShakeDisplacementsx = comp.NoShakeDisplacements    else  NoShakeDisplacementsx = NoShakeDisplacements   end 
        if  shakeFromShells      == nothing   shakeFromShellsx      = comp.shakeFromShells         else  shakeFromShellsx = shakeFromShells             end 
        if  shakeToShells        == nothing   shakeToShellsx        = comp.shakeToShells           else  shakeToShellsx = shakeToShells                 end 
        ##x if  steps                == nothing   stepsx                = comp.steps                   else  stepsx = steps                                 end 
    	
    	Computation(namex, nuclearModelx, gridx, asfSettingsx, approachx, processesx, initialConfigsx, initialLevelsx, 
    	            maxElectronLossx, NoShakeDisplacementsx, shakeFromShellsx, shakeToShellsx)
    end


    # `Base.string(comp::Cascade.Computation)`  ... provides a String notation for the variable comp::Cascade.Computation.
    function Base.string(comp::Cascade.Computation)
        sa = "Cascade computation   $(comp.name)  in $(comp.approach) approach  for Z = $(comp.nuclearModel.Z) and initial configurations: \n "
        for  config  in  comp.initialConfigs    sa = sa * string(config) * ",  "     end
        sa = sa * "\n Initial (levels, weights) are: "
        for  ilevel  in  comp.initialLevels     sa = sa * string(ilevel) * ",  "     end
        return( sa )
    end


    # `Base.show(io::IO, comp::Cascade.Computation)`  ... prepares a proper printout comp::Cascade.Computation.
    function Base.show(io::IO, comp::Cascade.Computation)
        sa = Base.string(comp)
        sa = sa * "\n ... in addition, the following parameters/settings are defined: ";       print(io, sa, "\n")
        println(io, "processes:                $(comp.processes)  ")
        println(io, "maxElectronLoss:          $(comp.maxElectronLoss)  ")
        println(io, "NoShakeDisplacements:     $(comp.NoShakeDisplacements)  ")
        println(io, "shakeFromShells:          $(comp.shakeFromShells)  ")
        println(io, "shakeToShells:            $(comp.shakeToShells)  ")
        ##x println(io, "steps:                    $(comp.steps)  \n")
        println(io, "nuclearModel:             $(comp.nuclearModel)  ")
        println(io, "grid:                     $(comp.grid)  ")
        println(io, "asfSettings:\n------------\n$(comp.asfSettings) ")
    end


    """
    `struct  Cascade.Data`  ... defines a type for an atomic cascade, i.e. lists of radiative, Auger and photoionization lines.

        + name           ::String                               ... A name for the cascade.
        + linesR         ::Array{PhotoEmission.Line,1}          ... List of radiative lines.
        + linesA         ::Array{AutoIonization.Line,1}         ... List of Auger lines.
        + linesP         ::Array{PhotoIonization.Line,1}        ... List of photoionization lines.
    """  
    struct  Data
        name             ::String
        linesR           ::Array{PhotoEmission.Line,1}
        linesA           ::Array{AutoIonization.Line,1}
        linesP           ::Array{PhotoIonization.Line,1}
    end 


    """
    `Cascade.Data()`  ... (simple) constructor for cascade data.
    """
    function Data()
        Data("", Array{PhotoEmission.Line,1}[], Array{AutoIonization.Line,1}[], Array{PhotoIonization.Line,1}[] )
    end


    # `Base.show(io::IO, data::Cascade.Data)`  ... prepares a proper printout of the variable data::Cascade.Data.
    function Base.show(io::IO, data::Cascade.Data) 
        println(io, "name:                    $(data.name)  ")
        println(io, "linesR:                  $(data.linesR)  ")
        println(io, "linesA:                  $(data.linesA)  ")
        println(io, "linesP:                  $(data.linesP)  ")
    end

    
    """
    `abstract type  Cascade.AbstractCascadeProperty`  
        ... defines an abstract and various singleton types for the different properties that can be obtained from the simulation of cascade data.

        + struct IonDist               ... simulate the 'ion distribution' as it is found after all cascade processes are completed.
        + struct FinalDist             ... simulate the 'final-level distribution' as it is found after all cascade processes are completed.
        + struct DecayPathes           ... determine the major 'decay pathes' of the cascade.
        + struct ElectronIntensity     ... simulate the electron-line intensities as function of electron energy.
        + struct PhotonIntensity       ... simulate  the photon-line intensities as function of electron energy. 
        + struct ElectronCoincidence   ... simulate electron-coincidence spectra.
    """
    abstract type  AbstractCascadeProperty                      end
    struct   IonDist              <:  AbstractCascadeProperty   end
    struct   FinalDist            <:  AbstractCascadeProperty   end
    struct   DecayPathes          <:  AbstractCascadeProperty   end
    struct   ElectronIntensity    <:  AbstractCascadeProperty   end
    struct   PhotonIntensity      <:  AbstractCascadeProperty   end
    struct   ElectronCoincidence  <:  AbstractCascadeProperty   end


    """
    abstract type  Cascade.AbstractSimulationMethod`  
        ... defines a enumeration for the various methods that can be used to run the simulation of cascade data.

        + struct ProbPropagation     ... to propagate the (occupation) probabilites of the levels until no further changes occur.
        + struct MonteCarlo          ... to simulate the cascade decay by a Monte-Carlo approach of possible pathes (not yet considered).
        + struct RateEquations       ... to solve the cascade by a set of rate equations (not yet considered).
    """
    abstract type  AbstractSimulationMethod                  end
    struct   ProbPropagation  <:  AbstractSimulationMethod   end
    struct   MonteCarlo       <:  AbstractSimulationMethod   end
    struct   RateEquations    <:  AbstractSimulationMethod   end


    """
    `struct  Cascade.SimulationSettings`  ... defines settings for performing the simulation of some cascade (data).

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
    `Cascade.SimulationSettings()`  ... constructor for an 'empty' instance of a Cascade.Block.
    """
    function SimulationSettings()
        SimulationSettings(0., 1.0e6,  0., 1.0e6 )
    end


    # `Base.show(io::IO, settings::SimulationSettings)`  ... prepares a proper printout of the variable settings::SimulationSettings.
    function Base.show(io::IO, settings::SimulationSettings) 
        println(io, "minElectronEnergy:        $(settings.minElectronEnergy)  ")
        println(io, "maxElectronEnergy:        $(settings.maxElectronEnergy)  ")
        println(io, "minPhotonEnergy:          $(settings.minPhotonEnergy)  ")
        println(io, "maxPhotonEnergy:          $(settings.maxPhotonEnergy)  ")
    end


    """
    `struct  Cascade.Simulation`  ... defines a simulation on some given cascade (data).

        + name            ::String                                   ... Name of the simulation
        + properties      ::Array{Cascade.AbstractCascadeProperty,1} ... Properties that are considered in this simulation of the cascade (data).
        + method          ::Cascade.AbstractSimulationMethod         ... Method that is used in the cascade simulation; cf. Cascade.SimulationMethod.
        + settings        ::Cascade.SimulationSettings               ... Settings for performing these simulations.
        + cascadeData     ::Cascade.Data                             ... Date on which the simulations are performed
    """
    struct  Simulation
        name              ::String
        properties        ::Array{Cascade.AbstractCascadeProperty,1}
        method            ::Cascade.AbstractSimulationMethod
        settings          ::Cascade.SimulationSettings 
        cascadeData       ::Cascade.Data 
    end 


    """
    `Cascade.Simulation()`  ... constructor for an 'default' instance of a Cascade.Simulation.
    """
    function Simulation()
        Simulation("Default cascade simulation", Cascade.AbstractCascadeProperty[], Cascade.ProbPropagation(), 
                   Cascade.SimulationSettings(), Cascade.Data() )
    end


    """
    `Cascade.Simulation(sim::Cascade.Simulation;`
        
                name=..,               properties=..,             method=..,              settings=..,     cascadeData=.. )
                
        ... constructor for re-defining the computation::Cascade.Computation.
    """
    function Simulation(sim::Cascade.Simulation;                              
        name::Union{Nothing,String}=nothing,                                  properties::Union{Nothing,Array{Cascade.AbstractCascadeProperty,1}}=nothing,
        method::Union{Nothing,Cascade.AbstractSimulationMethod}=nothing,      settings::Union{Nothing,Cascade.SimulationSettings}=nothing,    
        cascadeData::Union{Nothing,Cascade.Data}=nothing )
 
        if  name         == nothing   namex         = sim.name            else  namex = name                  end 
        if  properties   == nothing   propertiesx   = sim.properties      else  propertiesx = properties      end 
        if  method       == nothing   methodx       = sim.method          else  methodx = method              end 
        if  settings     == nothing   settingsx     = sim.settings        else  settingsx = settings          end 
        if  cascadeData  == nothing   cascadeDatax  = sim.cascadeData     else  cascadeDatax = cascadeData    end 
    	
    	Simulation(namex, propertiesx, methodx, settingsx, cascadeDatax)
    end


    # `Base.show(io::IO, simulation::Cascade.Simulation)`  ... prepares a proper printout of the variable simulation::Cascade.Simulation.
    function Base.show(io::IO, simulation::Cascade.Simulation) 
        println(io, "properties:        $(simulation.properties)  ")
        println(io, "method:            $(simulation.method)  ")
        println(io, "settings:          $(simulation.settings)  ")
    end


    """
    `struct  Cascade.LineIndex`  ... defines a line index with regard to the various lineLists of data::Cascade.Data.

        + process      ::Basics.AtomicProcess    ... refers to the particular lineList of cascade (data).
        + index        ::Int64                       ... index of the corresponding line.
    """
    struct  LineIndex
        process        ::Basics.AtomicProcess
        index          ::Int64 
    end 


    # `Base.show(io::IO, index::Cascade.LineIndex)`  ... prepares a proper printout of the variable index::Cascade.LineIndex.
    function Base.show(io::IO, index::Cascade.LineIndex) 
        println(io, "process:        $(index.process)  ")
        println(io, "index:          $(index.index)  ")
    end


    """
    `mutable struct  Cascade.Level`  ... defines a level specification for dealing with cascade transitions.

        + energy       ::Float64                     ... energy of the level.
        + J            ::AngularJ64                  ... total angular momentum of the level
        + parity       ::Basics.Parity               ... total parity of the level
        + NoElectrons  ::Int64                       ... total number of electrons of the ion to which this level belongs.
        + relativeOcc  ::Float64                     ... relative occupation  
        + parents      ::Array{Cascade.LineIndex,1}  ... list of parent lines that (may) populate the level.     
        + daugthers    ::Array{Cascade.LineIndex,1}  ... list of daugther lines that (may) de-populate the level.     
    """
    mutable struct  Level
        energy         ::Float64 
        J              ::AngularJ64 
        parity         ::Basics.Parity 
        NoElectrons    ::Int64 
        relativeOcc    ::Float64 
        parents        ::Array{Cascade.LineIndex,1} 
        daugthers      ::Array{Cascade.LineIndex,1} 
    end 



    # `Base.show(io::IO, level::Cascade.Level)`  ... prepares a proper printout of the variable level::Cascade.Level.
    function Base.show(io::IO, level::Cascade.Level) 
        println(io, "energy:        $(level.energy)  ")
        println(io, "J:             $(level.J)  ")
        println(io, "parity:        $(level.parity)  ")
        println(io, "NoElectrons:   $(level.NoElectrons)  ")
        println(io, "relativeOcc:   $(level.relativeOcc)  ")
        println(io, "parents:       $(level.parents)  ")
        println(io, "daugthers:     $(level.daugthers)  ")
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
        if   similarKey == LevelKey()    error("No similar level found !")   end
        
        for  line in linesR
            ##x println("** Radiative: similarKey =$similarKey     $(line.initialLevel)")
            if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
                rateR = rateR + line.photonRate;   NoPhotonLines = NoPhotonLines + 1    
                ##x println("NoPhotonLines = $NoPhotonLines")
            end
        end
        for  line in linesA
            ##x println("** Auger: similarKey =$similarKey     $(line.initialLevel)")
            if  similarKey == LevelKey( LevelSymmetry(line.initialLevel.J, line.initialLevel.parity), line.initialLevel.index, line.initialLevel.energy, 0.)
                rateA = rateA + line.totalRate;   NoAugerLines = NoAugerLines + 1    
            end
        end
        
        omegaR = rateR / (rateR + rateA);   omegaA = rateA / (rateR + rateA)
        newOutcome = DecayYield.Outcome(level, NoPhotonLines, NoAugerLines, rateR, rateA, omegaR, omegaA)
        return( newOutcome )
    end


    """
    `Cascade.computeSteps(comp::Cascade.Computation, stepList::Array{Cascade.Step,1})` 
        ... computes in turn all the requested transition amplitudes and PhotoEmission.Line's, AutoIonization.Line's, 
            etc. for all pre-specified decay steps of the cascade. When compared with the standard atomic process computations, 
            however, the amount of output is largely reduced. A set of  data::Cascade.Data  is returned.
    """
    function computeSteps(comp::Cascade.Computation, stepList::Array{Cascade.Step,1})
        linesA = AutoIonization.Line[];    linesR = PhotoEmission.Line[];    linesP = PhotoIonization.Line[]    
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        nt = 0;   st = 0
        for  step  in  stepList
            st = st + 1
            nc = length(step.initialMultiplet.levels) * length(step.finalMultiplet.levels) 
            println("\n  $st) Perform $(string(step.process)) amplitude computations for up to $nc lines (without selection rules): ")
            if  printSummary   println(iostream, "\n* $st) Perform $(string(step.process)) amplitude computations for " *
                                                 "up to $nc lines (without selection rules): ")   end      
            if      step.process == Basics.Auger 
                newLines = AutoIonization.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.nuclearModel, comp.grid, 
                                                              step.settings, output=true, printout=false) 
                append!(linesA, newLines);    nt = length(linesA)
            elseif  step.process == Basics.Radiative
                newLines = PhotoEmission.computeLinesCascade(step.finalMultiplet, step.initialMultiplet, comp.grid, 
                                                             step.settings, output=true, printout=false) 
                append!(linesR, newLines);    nt = length(linesR)
            else   error("Unsupported atomic process for cascade computations.")
            end
            println("     Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, giving now rise " *
                    "to a total of $nt $(string(step.process)) lines." )
            if  printSummary   println(iostream, "\n*    Step $st:: A total of $(length(newLines)) $(string(step.process)) lines are calculated, " *
                                                 "giving now rise to a total of $nt $(string(step.process)) lines." )   end      
        end
        #
        data = Cascade.Data(comp.name, linesR, linesA, linesP)
    end


    """
    `Cascade.determineSteps(comp::Cascade.Computation, blockList::Array{Cascade.Block,1})`  
        ... determines all steps::Cascade.Step that need to be computed for this cascade. It cycles through the given processes and 
            distinguished between the different cascade approaches. It also checks that the averaged energies of the configuration 
            allows such a step energetically. A stepList::Array{Cascade.Step,1} is returned
    """
    function determineSteps(comp::Cascade.Computation, blockList::Array{Cascade.Block,1})
        stepList = Cascade.Step[]
        ##x println("comp.approach = $(comp.approach)")
        if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
            for  a = 1:length(blockList)
                for  b = 1:length(blockList)
                    minEn = 100000.;   maxEn = -100000.
                    for  p = 1:length(blockList[a].multiplet.levels),  q = 1:length(blockList[b].multiplet.levels)
                        minEn = min(minEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                        maxEn = max(maxEn, blockList[a].multiplet.levels[p].energy - blockList[b].multiplet.levels[q].energy)
                    end
                    for  process  in  comp.processes
                        if      process == Basics.Radiative   
                            if  a == b   ||   minEn < 0.    continue   end
                            if  blockList[a].NoElectrons == blockList[b].NoElectrons
                                settings = PhotoEmission.Settings([E1], [UseBabushkin], false, false, false, Tuple{Int64,Int64}[], 0., 0., 1.0e6)
                                push!( stepList, Cascade.Step(process, settings, blockList[a].confs, blockList[b].confs, 
                                                              blockList[a].multiplet, blockList[b].multiplet) )
                            end
                        elseif  process == Basics.Auger       
                            if  a == b   ||   minEn < 0.    continue   end
                            if  blockList[a].NoElectrons == blockList[b].NoElectrons + 1
                                settings = AutoIonization.Settings(false, false, false, Tuple{Int64,Int64}[], 0., 1.0e6, 2, "Coulomb")
                                push!( stepList, Cascade.Step(process, settings, blockList[a].confs, blockList[b].confs, 
                                                              blockList[a].multiplet, blockList[b].multiplet) )
                            end
                        else    error("stop a")
                        end
                    end
                end
            end
            #
        else  error("Unsupported cascade approach.")
        end
        return( stepList )
    end



    """
    `Cascade.displayBlocks(stream::IO, blockList::Array{Cascade.Block,1})` 
        ... group & display the blocks of the cascade with same No. of electrons; this blocks are displayed with the
            minimum and maximum energy of each multiplet.
    """
    function displayBlocks(stream::IO, blockList::Array{Cascade.Block,1})
        #
        println(stream, "\n* Configuration 'blocks' (multiplets) in the given cascade model: \n")
        println(stream, "  ", TableStrings.hLine(134))
        println(stream, "      No.   Configurations                                                                       " *
                        "      Range of total energies " * TableStrings.inUnits("energy") ) 
        println(stream, "  ", TableStrings.hLine(134))
        i = 0
        for  block  in  blockList
            i = i + 1;    
            sa = "   " * TableStrings.flushright( 6, string(i); na=2)
            sb = " ";         for conf  in blockList[i].confs   sb = sb * string(conf) * ", "    end
            en = Float64[];   for level in  block.multiplet.levels    push!(en, level.energy)    end
            minEn = minimum(en);   minEn = Defaults.convertUnits("energy: from atomic", minEn)
            maxEn = maximum(en);   maxEn = Defaults.convertUnits("energy: from atomic", maxEn)
            sa = sa * TableStrings.flushleft(90, sb[1:end-2]; na=2) 
            sa = sa * TableStrings.flushleft(30, string( round(minEn)) * " ... " * string( round(maxEn)); na=2)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(134))

        return( nothing )
    end
   

    """
    `Cascade.displayInitialLevels(stream::IO, multiplet::Multiplet, initialLevels::Array{Tuple{Int64,Float64},1})`  
        ... display the calculated initial levels to screen together with their given relative occupation.
    """
    function displayInitialLevels(stream::IO, multiplet::Multiplet, initialLevels::Array{Tuple{Int64,Float64},1})
        println(stream, " ")
        println(stream, "* Initial levels of the given cascade, relative to the lowest, and their given occupation:")
        ##x println(stream, " ")
        println(stream, "  ", TableStrings.hLine(64))
        println(stream, "    Level  J Parity          Energy " * TableStrings.inUnits("energy") * "         rel. occupation ") 
        println(stream, "  ", TableStrings.hLine(64))
        for  i = 1:length(multiplet.levels)
            lev = multiplet.levels[i]
            en  = lev.energy - multiplet.levels[1].energy;    en_requested = Defaults.convertUnits("energy: from atomic", en)
            wx = 0.
            for  ilevel in initialLevels
                if  i in ilevel   wx = ilevel[2];    break    end
            end
            sb  = "          "   * string(wx)
            sc  = "   "  * TableStrings.level(i) * "     " * string(LevelSymmetry(lev.J, lev.parity)) * "     "
            @printf(stream, "%s %.15e %s %s", sc, en_requested, sb, "\n")
        end
        println(stream, "  ", TableStrings.hLine(64))
        return( nothing )
    end
   

    """
    `Cascade.displaySteps(stream::IO, steps::Array{Cascade.Step,1})` 
        ... displays all predefined steps in a neat table and supports to delete individual steps from the list.
    """
    function displaySteps(stream::IO, steps::Array{Cascade.Step,1})
        println(stream, " ")
        println(stream, "* Steps that are defined for the current cascade due to the given approach:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(170))
        sa = "  "
        sa = sa * TableStrings.center( 9, "Step-No"; na=2)
        sa = sa * TableStrings.flushleft(11, "Process"; na=1)
        sa = sa * TableStrings.flushleft(55, "Initial:  No CSF, configuration(s)"; na=4)
        sa = sa * TableStrings.flushleft(55, "Final:  No CSF, configuration(s)"; na=4)
        sa = sa * TableStrings.flushleft(40, "Energies from ... to in " * TableStrings.inUnits("energy"); na=4)
        println(stream, sa)
        println(stream, "  ", TableStrings.hLine(170))
        #
        for  i = 1:length(steps)
            sa = " " * TableStrings.flushright( 7, string(i); na=5)
            sa = sa  * TableStrings.flushleft( 11, string(steps[i].process); na=1)
            sb = "";   for conf in steps[i].initialConfigs   sb = sb * string(conf) * ", "    end
            sa = sa  * TableStrings.flushright( 5, string( length(steps[i].initialMultiplet.levels[1].basis.csfs) )*", "; na=0) 
            sa = sa  * TableStrings.flushleft( 50, sb[1:end-2]; na=4)
            sb = "";   for conf in steps[i].finalConfigs     sb = sb * string(conf) * ", "    end
            sa = sa  * TableStrings.flushright( 5, string( length(steps[i].finalMultiplet.levels[1].basis.csfs) )*", "; na=0) 
            sa = sa  * TableStrings.flushleft( 50, sb[1:end-2]; na=4)
            minEn = 1000.;   maxEn = -1000.;
            for  p = 1:length(steps[i].initialMultiplet.levels),  q = 1:length(steps[i].finalMultiplet.levels)
                minEn = min(minEn, steps[i].initialMultiplet.levels[p].energy - steps[i].finalMultiplet.levels[q].energy)
                maxEn = max(maxEn, steps[i].initialMultiplet.levels[p].energy - steps[i].finalMultiplet.levels[q].energy)
            end
            minEn = Defaults.convertUnits("energy: from atomic", minEn);   maxEn = Defaults.convertUnits("energy: from atomic", maxEn)
            sa = sa * string( round(minEn)) * " ... " * string( round(maxEn))
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(170))
    end
    
    
    
    """
    `Cascade.generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1})`  
        ... generate all block::Cascade.Block that need to be computed for this cascade and compute the corresponding multiplets.
            The different cascade approches follow different strategies in defining and computing these blocks. 
            A blockList::Array{Cascade.Block,1} is returned.
    """
    function generateBlocks(comp::Cascade.Computation, confs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital})
        blockList = Cascade.Block[]
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        if    comp.approach == AverageSCA()
            println("\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println("    + orbitals from the initial multiplet are applied throughout; ")
            println("    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println("    + only E1 dipole transitions are applied in all radiative decay stets; ")
            println("    + for each decay step, a (single) set of continuum orbitals with an configuration averaged energy is applied " *
                           "for all transitions of the same step. \n")
            if  printSummary   
            println(iostream, "\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println(iostream, "    + orbitals from the initial multiplet are applied throughout; ")
            println(iostream, "    + all blocks (multiplets) are generated from single-CSF levels and without any configuration mixing even in the SC; ")
            println(iostream, "    + only E1 dipole transitions are applied in all radiative decay stets; ")
            println(iostream, "    + for each decay step, a (single) set of continuum orbitals with an configuration averaged energy is applied " *
                              "for all transitions of the same step. \n")
            end
            #
            for  confa  in confs
                print("  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")
                if  printSummary   println(iostream, "\n*  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")   end
                multiplet = perform("computation: mutiplet from orbitals, no CI, CSF diagonal", [confa],  initalOrbitals, 
                                    comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
                println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
                ##x if  printSummary   println(iostream, "* ... and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")   end
            end
        elseif    comp.approach == SCA()
            println("\n  In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println("    + orbitals are generated independently for each multiplet (block); ")
            println("    + configuration interaction is included for each block; ")
            println("    + only E1 dipole transitions are applied in all radiative decay stets; \n")
            if  printSummary   
            println(iostream, "\n* In the cascade approach $(comp.approach), the following assumptions/simplifications are made: ")
            println(iostream, "    + orbitals are generated independently for each multiplet (block); ")
            println(iostream, "    + configuration interaction is included for each block; ")
            println(iostream, "    + only E1 dipole transitions are applied in all radiative decay stets; \n")
            end
            #
            i = 0
            for  confa  in confs
                ## i = i + 1;    if   i in [1,2, 4,5,6,7,8,9,10,11,12,13,14]  ||  i > 15   println("  Block $i omitted.");    continue    end
                ## i = i + 1;    if   i < 11  ||  i > 11   println("  Block $i omitted.");    continue    end
                print("  Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")
                if  printSummary   print(iostream, "* Multiplet computations for $(string(confa)[1:end]) with $(confa.NoElectrons) electrons ... ")   end
                basis     = perform("computation: SCF", [confa], comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                ##x multiplet = perform("computation: CI",  basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                multiplet = Basics.performCI(basis, comp.nuclearModel, comp.grid, comp.asfSettings; printout=false)
                push!( blockList, Cascade.Block(confa.NoElectrons, [confa], true, multiplet) )
                println("and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")
                ##x if  printSummary   println(iostream, "and $(length(multiplet.levels[1].basis.csfs)) CSF done. ")   end
            end
        else  error("Unsupported cascade approach.")
        end

        return( blockList )
    end
    
    

    """
    `Cascade.generateConfigurationList(initialConfs::Array{Configuration,1}, further::Int64, NoShake::Int64)`  
        ... generates all possible (decay) configurations with up to further holes and with NoShake displacements. First, all 
            configuratons are generated for which the hole is either moved 'outwards' or is moved and a second 'outer' hole is 
            created; this step is repated further + 2 times. From the generated list, only those configurations are kept with 
            up to further holes, when compared with the initial configuration. A confList::Array{Configuration,1} is returned.
    """
    function generateConfigurationList(initialConfs::Array{Configuration,1}, further::Int64, NoShake::Int64)
        confList = copy(initialConfs);    cList = copy(initialConfs);   initialNoElectrons = initialConfs[1].NoElectrons
        # First, move and generate new 'outer' hole without displacements
        for  fur = 1:further+1
            newConfList = Configuration[]
            for conf  in cList
                holeList = Basics.determineHoleShells(conf)
                for  holeShell in holeList
                    wa = generateConfigurationsWith1OuterHole(conf,  holeShell);   append!(newConfList, wa)
                    wa = generateConfigurationsWith2OuterHoles(conf, holeShell);   append!(newConfList, wa)
                end
            end
            if  length(newConfList) > 0    newConfList = Basics.excludeDoubles(newConfList)    end
            cList = newConfList
            append!(confList, newConfList)
        end
        # Make sure that only configurations with up to further holes are returned
        newConfList = Configuration[]
        for   conf in confList   
            if  conf.NoElectrons + further >= initialNoElectrons   push!(newConfList, conf)    end
        end
        # Add further shake-displacements if appropriate
        newConfList = Basics.excludeDoubles(newConfList)
        return( newConfList )
    end


    """
    `Cascade.generateConfigurationsWith1OuterHole(conf,  holeShell)`  
        ... generates all possible (decay) configurations where the hole in holeShell is moved 'outwards'. 
            A confList::Array{Configuration,1} is returned.
    """
    function generateConfigurationsWith1OuterHole(conf::Configuration,  holeShell::Shell)
         shList = Basics.generate("shells: ordered list for NR configurations", [conf]);   i0 = 0
         for  i = 1:length(shList)
             if   holeShell == shList[i]    i0 = i;    break    end
         end
         if  i0 == 0   error("stop a")   end
         #
         # Now move the hole 'outwards'
         confList = Configuration[]
         for  i = i0+1:length(shList)
             if  haskey(conf.shells, shList[i])  &&  conf.shells[ shList[i] ] >= 1  
                 newshells = copy( conf.shells )
                 newshells[ shList[i] ] = newshells[ shList[i] ] - 1
                 newshells[ holeShell ] = newshells[ holeShell ] + 1
                 push!(confList, Configuration( newshells, conf.NoElectrons ) )
             end
         end
         return( confList )
    end


    """
    `Cascade.generateConfigurationsWith2OuterHoles(conf,  holeShell)`  
        ... generates all possible (decay) configurations where the hole in holeShell is moved 'outwards'. 
            A confList::Array{Configuration,1} is returned.
    """
    function generateConfigurationsWith2OuterHoles(conf::Configuration,  holeShell::Shell)
         shList = Basics.generate("shells: ordered list for NR configurations", [conf]);   i0 = 0
         for  i = 1:length(shList)
             if   holeShell == shList[i]    i0 = i;    break    end
         end
         if  i0 == 0   error("stop a")   end
         #
         # Now move the hole 'outwards'
         confList = Configuration[]
         for  i = i0+1:length(shList)
             if  haskey(conf.shells, shList[i])  &&  conf.shells[ shList[i] ] >= 2  
                 newshells = copy( conf.shells )
                 newshells[ shList[i] ] = newshells[ shList[i] ] - 2
                 newshells[ holeShell ] = newshells[ holeShell ] + 1
                 push!(confList, Configuration( newshells, conf.NoElectrons - 1 ) )
             end
             #
             for  j = i0+1:length(shList)
                 if  i != j   &&   haskey(conf.shells, shList[i])  &&  conf.shells[ shList[i] ] >= 1   &&
                                   haskey(conf.shells, shList[j])  &&  conf.shells[ shList[j] ] >= 1 
                     newshells = copy( conf.shells )
                     newshells[ shList[i] ] = newshells[ shList[i] ] - 1
                     newshells[ shList[j] ] = newshells[ shList[j] ] - 1
                     newshells[ holeShell ] = newshells[ holeShell ] + 1
                     push!(confList, Configuration( newshells, conf.NoElectrons - 1 ) )
                 end
             end
         end
         return( confList )
    end



    """
    `Cascade.groupDisplayConfigurationList(Z::Float64, confs::Array{Configuration,1})` 
        ... group & display the configuration list into sublists with the same No. of electrons; this lists are displayed together 
            with an estimated total energy.
    """
    function groupDisplayConfigurationList(Z::Float64, confs::Array{Configuration,1})
        minNoElectrons = 1000;   maxNoElectrons = 0  
        for  conf in confs
            minNoElectrons = min(minNoElectrons, conf.NoElectrons)
            maxNoElectrons = max(maxNoElectrons, conf.NoElectrons)
        end
        #
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        #
        println("\n* Electron configuration used in the cascade:")
        @warn "*** Limit to just two configurations for each No. of electrons. ***"                        ## delete nxx
        if  printSummary   println(iostream, "\n* Electron configuration used in the cascade:")    end
        confList = Configuration[];   nc = 0
        for  n = maxNoElectrons:-1:minNoElectrons
            nxx = 0                                                                                        ## delete nxx
            println("\n  Configuration(s) with $n electrons:")
            if  printSummary   println(iostream, "\n    Configuration(s) with $n electrons:")      end
            for  conf in confs
                if n == conf.NoElectrons  
                    nxx = nxx + 1;    if nxx > 2   break    end                                            ## delete nxx
                    nc = nc + 1
                    push!(confList, conf ) 
                    wa = Semiempirical.estimate("binding energy", round(Int64, Z), conf);    wa = Defaults.convertUnits("energy: from atomic", wa)
                    sa = "   av. BE = "  * string( round(-wa) ) * "  " * TableStrings.inUnits("energy")
                    println("      " * string(conf) * sa * "      ($nc)" )
                    if  printSummary   println(iostream, "      " * string(conf) * sa * "      ($nc)")      end
                end  
            end
        end
        
        println("\n  A total of $nc configuration have been defined for this cascade, and selected configurations could be " *
                "removed here:  [currently not supported]")
        if  printSummary   println(iostream, "\n* A total of $nc configuration have been defined for this cascade, and selected " *
                                             "configurations could be removed here:  [currently not supported]")      end
        return( confList )
    end


    """
    `Cascade.modifySteps(stepList::Array{Cascade.Step,1})` 
        ... allows the user to modify the steps, for instance, by deleting selected steps of the cascade or by modifying the settings of
            one or several steps. A newStepList::Array{Cascade.Step,1} for which the transition data are eventually computed.
    """
    function modifySteps(stepList::Array{Cascade.Step,1})
        #
        newStepList = Cascade.Step[]
        #
        println("\n* Here, modify the individual steps explicitly in the code, if needed, ...... and just do it !!")
        # 
        #  Delete individual steps from stepList
        #  if  i in [1,2,5, ...] modify the particular settings, etc.
        for  i = 1:length(stepList)
            step = stepList[i]
            #
            if  i in []
                println("  Modify step $i :")
                newStep = Cascade.Step(step.process, step.settings, step.initialConfs, step.finalConfs, step.initialMultiplet, step.initialMultiplet)
                push!(newStepList, newStep)
            else
                push!(newStepList, step)
            end
        end
        #
        # wa = [1,2,3]
        # delete from list
        #
        println("\n  A total of $(length(newStepList)) steps are still defined in the cascade.")
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary   println(iostream, "\n* A total of $(length(newStepList)) steps are still defined in the cascade.")    end      
        
        return( newStepList )
    end


    """
    `Cascade.simulateLevelDistribution(simulation::Cascade.Simulation, data::Cascade.Data)` 
        ... sorts all levels as given by data and propagates their (occupation) probability until no further changes occur. For this 
            propagation, it runs through all levels and propagates the probabilty until no level probability changes anymore. The final 
            level distribution is then used to derive the ion distribution or the level distribution, if appropriate. Nothing is returned.
    """
    function simulateLevelDistribution(simulation::Cascade.Simulation, data::Cascade.Data)
        levels = Cascade.extractLevels(data)
        Cascade.displayLevelTree(data.name, levels, data)
        Cascade.propagateProbability!(levels, data)
        println("aaa")
        if   Cascade.IonDist()   in simulation.properties    Cascade.displayIonDistribution(data.name, levels)     end
        if   Cascade.FinalDist() in simulation.properties    Cascade.displayLevelDistribution(data.name, levels)   end

        return( nothing )
    end


    """
    `Cascade.displayIonDistribution(sc::String, levels::Array{Cascade.Level,1})` 
        ... displays the (current or final) ion distribution in a neat table. Nothing is returned.
    """
    function displayIonDistribution(sc::String, levels::Array{Cascade.Level,1})
        minElectrons = 1000;   maxElectrons = 0
        for  level in levels   minElectrons = min(minElectrons, level.NoElectrons);   maxElectrons = max(maxElectrons, level.NoElectrons)   end
        println(" ")
        println("  (Final) Ion distribution for the cascade:  $sc   mine = $minElectrons maxe = $maxElectrons ")
        println(" ")
        println("  ", TableStrings.hLine(31))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=4)        
        sa = sa * TableStrings.center(10,"Rel. occ.";      na=2)
        println(sa)
        println("  ", TableStrings.hLine(31))
        for n = maxElectrons:-1:minElectrons
            sa = "             " * string(n);   sa = sa[end-10:end];   prob = 0.
            for  level in levels    if  n == level.NoElectrons   prob = prob + level.relativeOcc    end    end
            sa = sa * "         " * @sprintf("%.5e", prob)
            println(sa)
        end
        println("  ", TableStrings.hLine(31))

        return( nothing )
    end


    """
    `Cascade.displayLevelDistribution(sc::String, levels::Array{Cascade.Level,1})` 
        ... displays the (current or final) level distribution in a neat table. Only those levels with a non-zero 
            occupation are displayed here. Nothing is returned.
    """
    function displayLevelDistribution(sc::String, levels::Array{Cascade.Level,1})
        minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels))
        for  i = 1:length(levels)
            minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
            energies[i]  = levels[i].energy   
        end
        enIndices = sortperm(energies, rev=true)
        # Now printout the results
        println(" ")
        println("  (Final) Level distribution for the cascade:  $sc")
        println(" ")
        println("  ", TableStrings.hLine(57))
        sa = "  "
        sa = sa * TableStrings.center(14, "No. electrons"; na=2)        
        sa = sa * TableStrings.center( 6, "J^P"          ; na=2);               
        sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=2)
        println(sa)
        println("  ", TableStrings.hLine(57))
        for n = maxElectrons:-1:minElectrons
            sa = "            " * string(n);   sa = sa[end-10:end]
            for  en in enIndices
                if  n == levels[en].NoElectrons  ##    &&  levels[en].relativeOcc > 0
                    sb = sa * "       " * string( LevelSymmetry(levels[en].J, levels[en].parity) )     * "    "
                    sb = sb * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", levels[en].energy))  * "      "
                    sb = sb * @sprintf("%.5e", levels[en].relativeOcc) 
                    sa = "           "
                    println(sb)
                end
            end
        end
        println("  ", TableStrings.hLine(57))

        return( nothing )
    end


    """
    `Cascade.displayLevelTree(sc::String, levels::Array{Cascade.Level,1}, data::Cascade.Data)` 
        ... displays all defined levels  in a neat table, together with their No. of electrons, symmetry, level energy, 
            current (relative) population as well as analogue information about their parents and daugther levels. This 
            enables one to recognize (and perhaps later add) missing parent and daughter levels. Nothing is returned.
    """
    function displayLevelTree(sc::String, levels::Array{Cascade.Level,1}, data::Cascade.Data)
        minElectrons = 1000;   maxElectrons = 0;   energies = zeros(length(levels))
        for  i = 1:length(levels)
            minElectrons = min(minElectrons, levels[i].NoElectrons);   maxElectrons = max(maxElectrons, levels[i].NoElectrons)
            energies[i]  = levels[i].energy   
        end
        enIndices = sortperm(energies, rev=true)
        # Now printout the results
        println(" ")
        println("  Level tree of this cascade:  $sc")
        println(" ")
        println("  ", TableStrings.hLine(175))
        sa = "  "
        sa = sa * TableStrings.center( 6, "No. e-"; na=2)        
        sa = sa * TableStrings.center( 6, "J^P"          ; na=2);               
        sa = sa * TableStrings.center(16, "Energy " * TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(10, "Rel. occ.";                                    na=5)
        sb = "Parents P(A: No_e, sym, energy) and Daughters D(R: No_e, sym, energy);  all energies in " * TableStrings.inUnits("energy")
        sa = sa * TableStrings.flushleft(100, sb; na=2)
        println(sa)
        println("  ", TableStrings.hLine(175))
        for n = maxElectrons:-1:minElectrons
            sa = "            " * string(n);   sa = sa[end-5:end]
            for  en in enIndices
                if  n == levels[en].NoElectrons
                    sb = sa * "      " * string( LevelSymmetry(levels[en].J, levels[en].parity) )      * "   "
                    sb = sb * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", levels[en].energy))  * "      "
                    sb = sb * @sprintf("%.5e", levels[en].relativeOcc)                                 * "    "
                    pProcessSymmetryEnergyList = Tuple{AtomicProcess,Int64,LevelSymmetry,Float64}[]
                    dProcessSymmetryEnergyList = Tuple{Basics.AtomicProcess,Int64,LevelSymmetry,Float64}[]
                    for  p in levels[en].parents
                        idx = p.index
                        if      p.process == Basics.Auger         lev = data.linesA[idx].initialLevel
                        elseif  p.process == Basics.Radiative    lev = data.linesR[idx].initialLevel
                        elseif  p.process == Basics.Photo         lev = data.linesP[idx].initialLevel
                        else    error("stop a")    end
                        push!( pProcessSymmetryEnergyList, (p.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                    end
                    for  d in levels[en].daugthers
                        idx = d.index
                        if      d.process == Basics.Auger         lev = data.linesA[idx].finalLevel
                        elseif  d.process == Basics.Radiative    lev = data.linesR[idx].finalLevel
                        elseif  d.process == Basics.Photo         lev = data.linesP[idx].finalLevel
                        else    error("stop b")    end
                        push!( dProcessSymmetryEnergyList, (d.process, lev.basis.NoElectrons, LevelSymmetry(lev.J, lev.parity), lev.energy) )
                    end
                    wa = TableStrings.processSymmetryEnergyTupels(120, pProcessSymmetryEnergyList, "P")
                    if  length(wa) > 0    sc = sb * wa[1];    println( sc )    else    println( sb )   end  
                    for  i = 2:length(wa)
                        sc = TableStrings.hBlank( length(sb) ) * wa[i];    println( sc )
                    end
                    wa = TableStrings.processSymmetryEnergyTupels(120, dProcessSymmetryEnergyList, "D")
                    for  i = 1:length(wa)
                        sc = TableStrings.hBlank( length(sb) ) * wa[i];    println( sc )
                    end
                    sa = "      "
                end
            end
        end
        println("  ", TableStrings.hLine(175))

        return( nothing )
    end


    """
    `Cascade.propagateProbability!(levels::Array{Cascade.Level,1}, data::Cascade.Data)` 
        ... propagates the relative level occupation through the levels of the cascade until no further change occur in the 
            relative level occupation. The argument levels is modified during the propagation, but nothing is returned otherwise.
    """
    function propagateProbability!(levels::Array{Cascade.Level,1}, data::Cascade.Data)
        n = 0
        println("\n  Probability propagation through $(length(levels)) levels of the cascade:")
        while true
            n = n + 1;    totalProbability = 0.
            print("    $n-th round ... ")
            for  level in levels
                if   level.relativeOcc > 0.   && length(level.daugthers) > 0
                    # A level with relative occupation > 0 has still 'daugther' levels; collect all decay rates for this level
                    prob  = level.relativeOcc;   totalProbability = totalProbability + prob;   rates = zeros(length(level.daugthers))
                    level.relativeOcc = 0.
                    for  i = 1:length(level.daugthers)
                        idx = level.daugthers[i].index
                        if      level.daugthers[i].process == Basics.Radiative    rate[i] = data.lineR[idx].photonRate.Babushkin
                        elseif  level.daugthers[i].process == Basics.Auger         rate[i] = data.lineA[idx].totalRate
                        else    error("stop a; process = $(level.daugthers[i].process) ")
                        end
                    end
                    totalRate = sum(rates)
                    # Shift the relative occupation to the 'daugther' levels due to the different decay pathes
                    for  i = 1:length(level.daugthers)
                        idx = level.daugthers[i].index
                        if      level.daugthers[i].process == Basics.Radiative    line = data.lineR[idx]
                        elseif  level.daugthers[i].process == Basics.Auger         line = data.lineA[idx]
                        else    error("stop b; process = $(level.daugthers[i].process) ")
                        end
                        level = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, 
                                               line.finalLevel.basis.NoElectrons, 0., Cascade.LineIndex[], Cascade.LineIndex[] )
                        kk    = Cascade.findLevelIndex(level, levels)
                        levels[kk].relativeOcc = levels[kk].relativeOcc + prob * rates[i] / totalRate
                    end
                end
            end
            println("has propagated a total of $totalProbability level occupation.")
            # Cycle once more if the relative occupation has still changed
            if  totalProbability == 0.    break    end
        end

        return( nothing )
    end


    """
    `Cascade.findLevelIndex(level::Cascade.Level, levels::Array{Cascade.Level,1})` 
        ... find the index of the given level within the given list of levels; an idx::Int64 is returned and an error message is 
            issued if the level is not found in the list.
    """
    function findLevelIndex(level::Cascade.Level, levels::Array{Cascade.Level,1})
        for  k = 1:length(levels)
            if  level.energy == levels[k].energy  &&   level.J == levels[k].J   &&   level.parity == levels[k].parity   &&
                level.NoElectrons == levels[k].NoElectrons
                kk = k;   return( kk )
            end
        end
        error("findLevelIndex():  No index was found;\n   level = $(level) ")
    end


    """
    `Cascade.extractLevels(data::Cascade.Data)` 
        ... extracts and sorts all levels from the given cascade data into a new levelList::Array{Cascade.Level,1} to simplify the 
            propagation of the probabilities. In this list, every level of the overall cascade just occurs just once, together 
            with its parent lines (which may populate the level) and the daugther lines (to which the pobability may decay). 
            A levelList::Array{Cascade.Level,1} is returned.
    """
    function extractLevels(data::Cascade.Data)
        levels = Cascade.Level[]
        print("\n  Extract, sort and unify the list of levels of the cascade ... ")
        for  i = 1:length(data.linesR)
            line = data.linesR[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(Basics.Radiative, i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(Basics.Radiative, i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end

        for  i = 1:length(data.linesA)
            line = data.linesA[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(Basics.Auger, i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(Basics.Auger, i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end

        for  i = 1:length(data.linesP)
            line = data.linesP[i]
            iLevel = Cascade.Level( line.initialLevel.energy, line.initialLevel.J, line.initialLevel.parity, line.initialLevel.basis.NoElectrons,
                                    line.initialLevel.relativeOcc, Cascade.LineIndex[], [ Cascade.LineIndex(Basics.Photo, i)] ) 
            Cascade.pushLevels!(levels, iLevel)  
            fLevel = Cascade.Level( line.finalLevel.energy, line.finalLevel.J, line.finalLevel.parity, line.finalLevel.basis.NoElectrons,
                                    line.finalLevel.relativeOcc, [ Cascade.LineIndex(Basics.Photo, i)], Cascade.LineIndex[] ) 
            Cascade.pushLevels!(levels, fLevel)  
        end

        println("a total of $(length(levels)) levels were found.")
        return( levels )
    end


    """
    `Cascade.pushLevels!(levels::Array{Cascade.Level,1}, newLevel::Cascade.Level)` 
        ... push's the information of newLevel of levels. This is the standard 'push!(levels, newLevel)' if newLevel is not yet 
            including in levels, and the proper modification of the parent and daugther lines of this level otherwise. The argument 
            levels::Array{Cascade.Level,1} is modified and nothing is returned otherwise.
    """
    function pushLevels!(levels::Array{Cascade.Level,1}, newLevel::Cascade.Level)
        for  i = 1:length(levels)
            if  newLevel.energy == levels[i].energy  &&  newLevel.J == levels[i].J  &&  newLevel.parity == levels[i].parity
                append!(levels[i].parents,   newLevel.parents)
                append!(levels[i].daugthers, newLevel.daugthers)
                return( nothing )
            end
        end
        push!( levels, newLevel)
        ##x info("... one level added, n = $(length(levels)) ")
        return( nothing )
    end

end # module


