
"""
`module  JAC.Atomic`  ... a submodel of JAC that contains all methods to set-up and process (simple) atomic and SCF computations as well as 
                          atomic cascade computations; it is using JAC, JAC.Radial, JAC.ManyElectron, JAC.Nuclear, JAC.Einstein, 
                          JAC.Hfs, JAC.IsotopeShift, JAC.PlasmaShift, JAC.LandeZeeman, JAC.PhotoEmission, JAC.PhotoIonization, JAC.AutoIonization.
"""
module Atomic

    using Interact
    using JAC, JAC.Radial, JAC.ManyElectron, JAC.Nuclear
    using JAC.Einstein, JAC.Hfs, JAC.IsotopeShift, JAC.PlasmaShift, JAC.LandeZeeman, JAC.AlphaVariation, JAC.FormFactor, JAC.DecayYield,
          JAC.GreenFunction,  JAC.MultipolePolarizibility,
          JAC.PhotoEmission, JAC.PhotoIonization, JAC.PhotoExcitation, JAC.AutoIonization ##x , JAC.ElectricDipoleMoment


    """
    `struct  Atomic.CasStep`  ... a type for defining an individual step of a (relativistic) complete active space computation for a specified 
                                  set of levels. An instance of this struct provides all information to generate the atomic basis and to perform 
                                  the associated SCF and multiplet computations for a selected No. of levels.

        + name             ::String                    ... to assign a name to the given SCF step.
        + NoElectrons      ::Int64                     ... Number of electrons.
        + refConfigs       ::Array{Configuration,1}    ... List of references configurations, at least 1.
        + Sfrom            ::Array{Shell,1}            ... Single-excitations from shells   [sh_1, sh_2, ...]
        + Sto              ::Array{Shell,1}            ... Single-excitations to shells  [sh_1, sh_2, ...]
        + Dfrom            ::Array{Shell,1}            ... Double-excitations from shells   [sh_1, sh_2, ...]
        + Dto              ::Array{Shell,1}            ... Double-excitations to shells  [sh_1, sh_2, ...]
        + Tfrom            ::Array{Shell,1}            ... Triple-excitations from shells   [sh_1, sh_2, ...]
        + Tto              ::Array{Shell,1}            ... Triple-excitations to shells  [sh_1, sh_2, ...]
        + Qfrom            ::Array{Shell,1}            ... Quadrupole-excitations from shells   [sh_1, sh_2, ...]
        + Qto              ::Array{Shell,1}            ... Quadrupole-excitations to shells  [sh_1, sh_2, ...]
        + restrictions     ::Array{String,1}           ... List of Strings to define 'restrictions' to applied to the CSF basis.
        + NoIterations     ::Int64                     ... Number of SCf iterations to be applied in this step of computations.
        + frozenSubshells  ::Array{Subshell,1}         ... List of subshells that are kept 'frozen' in this step.
        + failureHandling  ::Array{String,1}           ... List of Strings to define subsequent steps in case of failure.
    """
    struct  CasStep
        name               ::String                    
        NoElectrons        ::Int64  
        refConfigs         ::Array{Configuration,1}
        Sfrom              ::Array{Shell,1}
        Sto                ::Array{Shell,1}
        Dfrom              ::Array{Shell,1}
        Dto                ::Array{Shell,1} 
        Tfrom              ::Array{Shell,1}
        Tto                ::Array{Shell,1}
        Qfrom              ::Array{Shell,1}
        Qto                ::Array{Shell,1}
        restrictions       ::Array{String,1} 
        NoIterations       ::Int64 
        frozenShells       ::Array{Subshell,1}
        failureHandling    ::Array{String,1}
    end


    """
    `JAC.Atomic.CasStep(name::String, NoElectrons::Int64)`  ... constructor for an 'initial' instance of a variable::Atomic.CasStep with just 
                                                                a given name and No. of Electrons.
    """
    function CasStep(name::String, NoElectrons::Int64)
        CasStep(name, NoElectrons, Configuration[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], 
                String[], 0, Subshell[], String[])
    end


    """
    `JAC.Atomic.CasStep("interactive"; refine::Atomic.CasStep=..)`  ... constructor to generate a new instance of Atomic.CasStep
                                       interactively by replying to some detailed dialog. This constructor enables one to add, modify or
                                       refine an individual step in the complete active space computations. **Not yet implemented !**
    """
    function CasStep(sa::String; refine::Atomic.CasStep = Atomic.CasStep())
        sa != "interactive"   &&    error("Unsupported keystring = $sa.")

        #= print("Enter the nuclear charge Z:  ");    Zx = Base.parse(Float64, readline(STDIN) )

        yes = yesno("The ...; modify ? [N|y]") 
        while  yes     
            print("Enter a charge distribution model {point, Fermi, uniform}:  ");    modelx = strip(readline(STDIN))
            !(modelx in [""])   &&   println("Unsupported charge distribution model $modelx ... redo:")
        end =#
    end


    # `Base.show(io::IO, step::Atomic.CasStep)`  ... prepares a proper printout of the (individual step of computations) step::Atomic.CasStep.
    function Base.show(io::IO, step::Atomic.CasStep)
        sa = Base.string(step);                print(io, sa, "\n")
        sa = "Reference configurations: ";     print(io, sa, "\n")   
        if  length(step.refConfigs) > 0
           for  conf  in step.refConfigs          print(io, conf )  end;   print(io, "\n")
        end
        if  length(step.Sfrom) > 0
           sa = "Singles from:              { ";      for  sh in step.Sfrom    sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] *  " }   ... to { ";      for  sh in step.Sto      sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] *  " }";            print(io, sa, "\n")
        end   
        if  length(step.Dfrom) > 0
           sa = "Doubles from:              { ";      for  sh in step.Dfrom    sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] *  " }   ... to { ";      for  sh in step.Dto      sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] *  " }";            print(io, sa, "\n")
        end   
        if  length(step.Tfrom) > 0
           sa = "Triples from:              { ";      for  sh in step.Tfrom    sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] *  " }   ... to { ";      for  sh in step.Tto      sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] *  " }";            print(io, sa, "\n")
        end   
        if  length(step.Qfrom) > 0
           sa = "Quadruples from:           { ";      for  sh in step.Qfrom    sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] *  " }   ... to { ";      for  sh in step.Qto      sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] *  " }";            print(io, sa, "\n")
        end   
    end


    # `Base.string(step::Atomic.CasStep)`  ... provides a String notation for the variable step::Atomic.CasStep.
    function Base.string(step::Atomic.CasStep)
        sa = "Step: Computational model for $(step.name) with $(step.NoElectrons) electrons and with (step.NoCoreShells) "
        sa = sa * "core shells:"
        return( sa )
    end



    """
    `struct  Atomic.CasComputation`  ... a struct for defining a series of complete active space computations for a specified set of levels. 
                                         This type provides all information to generate the atomic basis and to perform the corresponding
                                         SCF and multiplet computations for a selected No. of levels.

        + name             ::String                     ... to assign a name to the given model.
        + previousStep     ::Int64                      ... 0, if no previous steps was yet done in this series.
        + steps            ::Array{Atomic.CasStep,1}    ... List of SCF steps that are to be done in this model computation.
    """
    struct  CasComputation
        name               ::String                   
        previousStep       ::Int64 
        steps              ::Array{Atomic.CasStep,1}
    end


    """
    `JAC.Atomic.CasComputation()`  ... constructor for an 'empty' instance of the a variable::Atomic.CasComputation
    """
    function CasComputation()
        CasComputation("", 0, Array{CasStep,1}[])
    end


    """
    `JAC.Atomic.CasComputation("interactive"; refine::Atomic.CasComputation=..)`  ... constructor to generate a new instance of
                               Atomic.CasComputation interactively by replying to some detailed dialog. This constructor enables one to add,
                               modify or refine the model for a particular complete active space computations. **Not yet implemented !**
    """
    function CasComputation(sa::String; refine::Atomic.CasComputation = Atomic.CasComputation())
        sa != "interactive"   &&    error("Unsupported keystring = $sa.")

        #= print("Enter the nuclear charge Z:  ");    Zx = Base.parse(Float64, readline(STDIN) )

        yes = yesno("The ...; modify ? [N|y]") 
        while  yes     
            print("Enter a charge distribution model {point, Fermi, uniform}:  ");    modelx = strip(readline(STDIN))
            !(modelx in [""])   &&   println("Unsupported charge distribution model $modelx ... redo:")
        end =#
    end


    # `Base.show(io::IO, model::CasComputation)`  ... prepares a proper printout of the (individual step of computations) model::CasComputation.
    function Base.show(io::IO, model::CasComputation)
        sa = Base.string(model)
        print(io, sa)
    end


    # `Base.string(model::CasComputation)`  ... provides a String notation for the variable model::CasComputation.
    function Base.string(model::CasComputation)
        sa = "Model: $(model.name) includes $(length(model.steps)) steps and has been calculated up to the step $(model.previousStep)."
        return( sa )
    end



    """
    `struct  Atomic.CasSettings`  ... a struct for defining the settings for complete active space computations

        + generateScf       ::Bool               ... True, if a SCF need to be generated, and false if just the start orbitals should 
                                                     be applied.
        + startOrbitals     ::String             ... Specify how the start orbitals are obtained ["fromNRorbitals", "fromGrasp", "hydrogenic"].
        + rwfFilename       ::String             ... Filename of orbitals, if taken from Grasp.
        + levels            ::Array{Int64,1}     ... Levels on which the optimization need to be carried out.
        + includeBreit      ::Bool               ... True, if the Breit interaction is included into the SCF generation.
        + maxIterations     ::Int64              ... maximum number of SCF iterations
        + scfAccuracy       ::Float64            ... convergence criterion for the SCF field.
        + iterationSequence ::Array{Subshell,1}  ... Sequence of subshells to be optimized.
    """
    struct  CasSettings
        generateScf         ::Bool
        startOrbitals       ::String
        rwfFilename         ::String
        levels              ::Array{Int64,1}
        includeBreit        ::Bool
        maxIterations       ::Int64
        scfAccuracy         ::Float64
        iterationSequence   ::Array{Subshell,1}
    end


    """
    `JAC.Atomic.CasSettings()`  ... constructor for setting the default values.
    """
    function CasSettings()
        CasSettings(false, "fromNRorbitals", "", Int64[1], false, 24, 1.0e-8, Subshell[] )
    end


    # `Base.show(io::IO, settings::Atomic.CasSettings)`  ... prepares a proper printout of the settings::Atomic.CasSettings.
    function Base.show(io::IO, settings::Atomic.CasSettings)
        println(io, "generateScf:        $(settings.generateScf)  ")
        println(io, "startOrbitals:      $(settings.startOrbitals)  ")
        println(io, "rwfFilename:        $(settings.rwfFilename)  ")
        println(io, "levels:             $(settings.levels)  ")
        println(io, "includeBreit:       $(settings.includeBreit)  ")
        println(io, "maxIterations:      $(settings.maxIterations)  ")
        println(io, "scfAccuracy:        $(settings.scfAccuracy)  ")
        println(io, "iterationSequence:  $(settings.iterationSequence)  ")
    end


    # `Base.string(settings::Atomic.CasSettings)`  ... provides a String notation for the variable settings::Atomic.CasSettings.
    function Base.string(settings::Atomic.CasSettings)
        error("Not yet implemented.")
        sa = "Cas settings: maximum No. of iterations = $(settings.maxIterations), accuracy = (settings.scfAccuracy)"
        return( sa )
    end


    """
    `struct  Computation`  ... defines a type for defining  (the model of simple) atomic computation of a single multiplet, 
                               including the SCF and CI as well as level properties and transition property calculations.

        + name                           ::String                        ... A name associated to the computation.
        + nuclearModel                   ::Nuclear.Model          ... Model, charge and parameters of the nucleus.
        + grid                           ::Radial.Grid                   ... The radial grid to be used for the computation.
        ##x + calcLevelProperties            ::Bool                          ... True, if level structures and properties are to be calculated
        + properties                     ::Array{AtomicLevelProperty,1}  ... List of atomic properties to be calculated.
        + configs                        ::Array{Configuration,1}        ... A list of non-relativistic configurations.
        + asfSettings                    ::AsfSettings                   ... Provides the settings for the SCF process and for the CI and QED calculations.
        + initialConfigs                 ::Array{Configuration,1}        ... A list of initial-state configurations for some transition 
                                                                             property calculation, such as radiative transition, Auger, etc. 
        + initialAsfSettings             ::AsfSettings                   ... Provides the SCF and CI settings for the initial-state multiplet.
        + intermediateConfigs            ::Array{Configuration,1}        ... A list of initial-state configurations.
        + intermediateAsfSettings        ::AsfSettings                   ... Provides the SCF settings for the intermediate-state multiplet.
        + finalConfigs                   ::Array{Configuration,1}        ... A list of final-state configurations.
        + finalAsfSettings               ::AsfSettings                   ... Provides the SCF and CI settings for the final-state multiplet.
        + alphaSettings                  ::AlphaVariation.Settings       ... Settings for alpha-variation parameter calculations.
        + einsteinSettings               ::Einstein.Settings             ... Settings for Einstein coefficient calculations.
        + formSettings                   ::FormFactor.Settings           ... Settings for atomic form factor calculations.
        + greenSettings                  ::GreenFunction.Settings        ... Settings for approximate Green function calculations.
        + hfsSettings                    ::Hfs.Settings                  ... Settings for hyperfine parameter calculations.
        + isotopeSettings                ::IsotopeShift.Settings         ... Settings for isotope shift parameter calculations.
        + plasmaSettings                 ::PlasmaShift.Settings          ... Settings for plasma-shift calculations.
        + polaritySettings               ::MultipolePolarizibility.Settings   ... Settings for plasma-shift calculations.
        + yieldSettings                  ::DecayYield.Settings           ... Settings for fluoresence and Auger yield calculations.
        + zeemanSettings                 ::LandeZeeman.Settings          ... Settings for Lande-Zeeman coefficient calculations.
        + process                        ::JAC.AtomicProcess             ... An (additional) process for which the properties are to be evaluated 
                                                                             for the given initial- and final-state configurations.
        + processSettings                ::Union{JAC.PhotoEmission.Settings, JAC.AutoIonization.Settings, JAC.PlasmaShift.AugerSettings, 
                                                 JAC.PhotoIonization.Settings, JAC.PlasmaShift.PhotoSettings, 
                                                 JAC.PhotoExcitation.Settings, JAC.PhotoExcitationAutoion.Settings, JAC.PhotoRecombination.Settings, 
                                                 JAC.ImpactExcitation.Settings, JAC.Dielectronic.Settings, RadiativeAuger.Settings,
                                                 JAC.PairAnnihilation1Photon.Settings, JAC.ImpactExcitationAutoion.Settings, 
                                                 JAC.MultiPhotonDeExcitation.Settings, JAC.CoulombExcitation.Settings, 
                                                 JAC.CoulombIonization.Settings} ... Provides the settings for the selected process.
    """
    struct  Computation
        name                           ::String
        nuclearModel                   ::Nuclear.Model
        grid                           ::Radial.Grid
        ##x calcLevelProperties            ::Bool
        properties                     ::Array{AtomicLevelProperty,1}
        configs                        ::Array{Configuration,1}
        asfSettings                    ::AsfSettings
        initialConfigs                 ::Array{Configuration,1} 
        initialAsfSettings             ::AsfSettings
        intermediateConfigs            ::Array{Configuration,1} 
        intermediateAsfSettings        ::AsfSettings
        finalConfigs                   ::Array{Configuration,1}
        finalAsfSettings               ::AsfSettings
        alphaSettings                  ::AlphaVariation.Settings
        einsteinSettings               ::Einstein.Settings
        formSettings                   ::FormFactor.Settings
        greenSettings                  ::GreenFunction.Settings
        hfsSettings                    ::Hfs.Settings
        isotopeSettings                ::IsotopeShift.Settings
        plasmaSettings                 ::PlasmaShift.Settings
        polaritySettings               ::MultipolePolarizibility.Settings
        yieldSettings                  ::DecayYield.Settings
        zeemanSettings                 ::LandeZeeman.Settings
        process                        ::JAC.AtomicProcess
        processSettings                ::Union{JAC.PhotoEmission.Settings, JAC.AutoIonization.Settings, JAC.PlasmaShift.AugerSettings, 
                                               JAC.PhotoIonization.Settings, JAC.PlasmaShift.PhotoSettings,
                                               JAC.PhotoRecombination.Settings, JAC.Dielectronic.Settings, JAC.ImpactExcitation.Settings,
                                               JAC.CoulombExcitation.Settings, JAC.CoulombIonization.Settings, 
                                               JAC.PhotoExcitation.Settings, JAC.PhotoExcitationFluores.Settings, 
                                               JAC.PhotoExcitationAutoion.Settings, JAC.RayleighCompton.Settings, 
                                               JAC.MultiPhotonDeExcitation.Settings, JAC.PhotoIonizationFluores.Settings, 
                                               JAC.PhotoIonizationAutoion.Settings, JAC.ImpactExcitationAutoion.Settings,
                                               RadiativeAuger.Settings, JAC.MultiPhotonIonization.Settings,
                                               JAC.MultiPhotonDoubleIon.Settings, JAC.InternalConversion.Settings,} #= , 
                                               #
                                               JAC.PairAnnihilation1Photon.Settings } =#
    end 


    """
    `JAC.Atomic.Computation()`  ... constructor for an 'empty' instance::Atomic.Computation.
    """
    function Computation()
        Computation("", Nuclear.Model(1.), Radial.Grid(), AtomicLevelProperty[], 
                    Configuration[], AsfSettings(),
                    Configuration[], AsfSettings(),
                    Configuration[], AsfSettings(),
                    Configuration[], AsfSettings(),
                    AlphaVariation.Settings(), Einstein.Settings(), FormFactor.Settings(), GreenFunction.Settings(), 
                    Hfs.Settings(), 
                    IsotopeShift.Settings(), PlasmaShift.Settings(), MultipolePolarizibility.Settings(), 
                    DecayYield.Settings(),  LandeZeeman.Settings(), 
                    ##x NoAmplitude, ElectricDipoleMoment.Settings(), 
                    NoProcess, PhotoEmission.Settings() )
    end


    """
    `JAC.Atomic.Computation(sa::String, nm::Nuclear.Model; grid::Radial.Grid = JAC.Radial.Grid('grid: exponential'), ##x calc::Bool = false,
         properties::Array{AtomicLevelProperty,1} = AtomicLevelProperty[], configs::Array{Configuration,1} = Configuration[],
         asfSettings::AsfSettings = AsfSettings(),  
         initialConfigs::Array{Configuration,1} = Configuration[], initialAsfSettings::AsfSettings = AsfSettings(), 
         intermediateConfigs::Array{Configuration,1} = Configuration[], intermediateAsfSettings::AsfSettings = AsfSettings(), 
         finalConfigs::Array{Configuration,1} = Configuration[], finalAsfSettings::AsfSettings = AsfSettings(), 
         alphaSettings::AlphaVariation.Settings = AlphaVariation.Settings(), einsteinSettings::Einstein.Settings = Einstein.Settings(), 
         formSettings::FormFactor.Settings = FormFactor.Settings(), greenSettings::GreenFunction.Settings = GreenFunction.Settings(), 
         hfsSettings::Hfs.Settings = Hfs.Settings(),
         isotopeSettings::IsotopeShift.Settings = IsotopeShift.Settings(), plasmaSettings::PlasmaShift.Settings = PlasmaShift.Settings(), 
         polaritySettings::MultipolePolarizibility.Settings = MultipolePolarizibility.Settings(), 
         yieldSettings::DecayYield.Settings = DecayYield.Settings(), zeemanSettings::LandeZeeman.Settings = LandeZeeman.Settings(), 
         ##x amplitude::JAC.AtomicAmplitude = JAC.NoAmplitude, amplitudeSettings::Any = true, 
         process::JAC.AtomicProcess = JAC.NoProcess, processSettings::Any = true)`   ... constructor for an instance::Atomic.Computation for 
         which all requested details are given by proper keyword specification. A few internal checks are  made. 
    """
    function Computation(sa::String, nm::Nuclear.Model; grid::Radial.Grid = JAC.Radial.Grid("grid: exponential"), ##x calc::Bool = false,
                         properties::Array{AtomicLevelProperty,1} = [JAC.NoProperty], configs::Array{Configuration,1} = Configuration[],
                         asfSettings::AsfSettings = AsfSettings(),  
                         initialConfigs::Array{Configuration,1} = Configuration[], initialAsfSettings::AsfSettings = AsfSettings(), 
                         intermediateConfigs::Array{Configuration,1} = Configuration[], intermediateAsfSettings::AsfSettings = AsfSettings(), 
                         finalConfigs::Array{Configuration,1} = Configuration[], finalAsfSettings::AsfSettings = AsfSettings(), 
                         alphaSettings::AlphaVariation.Settings = AlphaVariation.Settings(), 
                         einsteinSettings::Einstein.Settings = Einstein.Settings(), 
                         formSettings::FormFactor.Settings = FormFactor.Settings(), 
                         greenSettings::GreenFunction.Settings = GreenFunction.Settings(), 
                         hfsSettings::Hfs.Settings = Hfs.Settings(),
                         isotopeSettings::IsotopeShift.Settings = IsotopeShift.Settings(), 
                         plasmaSettings::PlasmaShift.Settings = PlasmaShift.Settings(), 
                         polaritySettings::MultipolePolarizibility.Settings = MultipolePolarizibility.Settings(), 
                         yieldSettings::DecayYield.Settings = DecayYield.Settings(), 
                         zeemanSettings::LandeZeeman.Settings = LandeZeeman.Settings(), 
                         ##x amplitude::JAC.AtomicAmplitude = JAC.NoAmplitude, amplitudeSettings::Any = true, 
                         process::JAC.AtomicProcess = JAC.NoProcess, processSettings::Any = true)
        if      process == JAC.NoProcess         && typeof(processSettings) == Bool && processSettings             procSettings = PhotoEmission.Settings()
        elseif  process == JAC.Auger             && typeof(processSettings) == Bool && processSettings             procSettings = AutoIonization.Settings()
        elseif  process == JAC.Auger             && typeof(processSettings) == AutoIonization.Settings             procSettings = processSettings 
        elseif  process == JAC.AugerInPlasma     && typeof(processSettings) == Bool && processSettings             procSettings = PlasmaShift.AugerSettings()
        elseif  process == JAC.AugerInPlasma     && typeof(processSettings) == PlasmaShift.AugerSettings           procSettings = processSettings 
        elseif  process == JAC.Radiative        && typeof(processSettings) == Bool && processSettings             procSettings = PhotoEmission.Settings()
        elseif  process == JAC.Radiative        && typeof(processSettings) == PhotoEmission.Settings                  procSettings = processSettings 
        elseif  process == JAC.PhotoExc          && typeof(processSettings) == Bool && processSettings             procSettings = PhotoExcitation.Settings()
        elseif  process == JAC.PhotoExc          && typeof(processSettings) == PhotoExcitation.Settings            procSettings = processSettings 
        elseif  process == JAC.Photo             && typeof(processSettings) == Bool && processSettings             procSettings = PhotoIonization.Settings()
        elseif  process == JAC.Photo             && typeof(processSettings) == PhotoIonization.Settings            procSettings = processSettings 
        elseif  process == JAC.PhotoInPlasma     && typeof(processSettings) == Bool && processSettings             procSettings = PlasmaShift.PhotoSettings()
        elseif  process == JAC.PhotoInPlasma     && typeof(processSettings) == PlasmaShift.PhotoSettings           procSettings = processSettings 
        elseif  process == JAC.Rec               && typeof(processSettings) == Bool && processSettings             procSettings = PhotoRecombination.Settings()
        elseif  process == JAC.Rec               && typeof(processSettings) == PhotoRecombination.Settings         procSettings = processSettings 
        elseif  process == JAC.Dierec            && typeof(processSettings) == Bool && processSettings             procSettings = Dielectronic.Settings()
        elseif  process == JAC.Dierec            && typeof(processSettings) == Dielectronic.Settings               procSettings = processSettings 
        elseif  process == JAC.PhotoExcFluor     && typeof(processSettings) == Bool && processSettings             procSettings = PhotoExcitationFluores.Settings()
        elseif  process == JAC.PhotoExcFluor     && typeof(processSettings) == PhotoExcitationFluores.Settings     procSettings = processSettings 
        elseif  process == JAC.PhotoExcAuto      && typeof(processSettings) == Bool && processSettings             procSettings = PhotoExcitationAutoion.Settings()
        elseif  process == JAC.PhotoExcAuto      && typeof(processSettings) == PhotoExcitationAutoion.Settings     procSettings = processSettings 
        elseif  process == JAC.PhotoIonFluor     && typeof(processSettings) == Bool && processSettings             procSettings = PhotoIonizationFluores.Settings()
        elseif  process == JAC.PhotoIonFluor     && typeof(processSettings) == PhotoIonizationFluores.Settings     procSettings = processSettings 
        elseif  process == JAC.PhotoIonAuto      && typeof(processSettings) == Bool && processSettings             procSettings = PhotoIonizationAutoion.Settings()
        elseif  process == JAC.PhotoIonAuto      && typeof(processSettings) == PhotoIonizationAutoion.Settings     procSettings = processSettings 
        elseif  process == JAC.MultiPhotonDE     && typeof(processSettings) == Bool && processSettings             procSettings = MultiPhotonDeExcitation.Settings()
        elseif  process == JAC.MultiPhotonDE     && typeof(processSettings) == MultiPhotonDeExcitation.Settings    procSettings = processSettings 
        elseif  process == JAC.Compton           && typeof(processSettings) == Bool && processSettings             procSettings = RayleighCompton.Settings()
        elseif  process == JAC.Compton           && typeof(processSettings) == RayleighCompton.Settings            procSettings = processSettings 
        elseif  process == JAC.Eimex             && typeof(processSettings) == Bool && processSettings             procSettings = ImpactExcitation.Settings()
        elseif  process == JAC.Eimex             && typeof(processSettings) == ImpactExcitation.Settings           procSettings = processSettings 
        elseif  process == JAC.ImpactExcAuto     && typeof(processSettings) == Bool && processSettings             procSettings = ImpactExcitationAutoion.Settings()
        elseif  process == JAC.ImpactExcAuto     && typeof(processSettings) == ImpactExcitationAutoion.Settings    procSettings = processSettings 
        elseif  process == JAC.RAuger            && typeof(processSettings) == Bool && processSettings             procSettings = RadiativeAuger.Settings()
        elseif  process == JAC.RAuger            && typeof(processSettings) == RadiativeAuger.Settings             procSettings = processSettings 
        elseif  process == JAC.MultiPI           && typeof(processSettings) == Bool && processSettings             procSettings = MultiPhotonIonization.Settings()
        elseif  process == JAC.MultiPI           && typeof(processSettings) == MultiPhotonIonization.Settings      procSettings = processSettings 
        elseif  process == JAC.MultiPDI          && typeof(processSettings) == Bool && processSettings             procSettings = MultiPhotonDoubleIon.Settings()
        elseif  process == JAC.MultiPDI          && typeof(processSettings) == MultiPhotonDoubleIon.Settings       procSettings = processSettings 
        elseif  process == JAC.InternalConv      && typeof(processSettings) == Bool && processSettings             procSettings = InternalConversion.Settings()
        elseif  process == JAC.InternalConv      && typeof(processSettings) == InternalConversion.Settings         procSettings = processSettings 
        elseif  process == JAC.Coulex            && typeof(processSettings) == Bool && processSettings             procSettings = CoulombExcitation.Settings()
        elseif  process == JAC.Coulex            && typeof(processSettings) == CoulombExcitation.Settings          procSettings = processSettings 
        elseif  process == JAC.Coulion           && typeof(processSettings) == Bool && processSettings             procSettings = CoulombIonization.Settings()
        elseif  process == JAC.Coulion           && typeof(processSettings) == CoulombIonization.Settings          procSettings = processSettings 
        else    error("The processSettings must fit to the given process or should not occur if NoProcess is to be calculated.")
        end
        
        Computation(sa, nm, grid, properties, configs, asfSettings, 
                    initialConfigs, initialAsfSettings, 
                    intermediateConfigs, intermediateAsfSettings,
                    finalConfigs, finalAsfSettings,
                    alphaSettings, einsteinSettings, formSettings, greenSettings, hfsSettings, isotopeSettings, plasmaSettings, 
                    polaritySettings, yieldSettings, zeemanSettings, 
                    ##x amplitude, ampSettings, 
                    process, procSettings) 
    end


    # `Base.show(io::IO, computation::Atomic.Computation)`  ... prepares a proper printout of the variable computation::Atomic.Computation.
    function Base.show(io::IO, computation::Atomic.Computation) 
        println(io, "name:                           $(computation.name)  ")
        println(io, "nuclearModel:                   $(computation.nuclearModel)  ")
        println(io, "grid:                           $(computation.grid)  ")
        ##x println(io, "calcLevelProperties:            $(computation.calcLevelProperties)  ")
        println(io, "properties:                     $(computation.properties)  ")
        println(io, "configs:                        $(computation.configs)  ")
        println(io, "asfSettings:                      computation.asfSettings  ")
        println(io, "initialConfigs:                 $(computation.initialConfigs)  ")
        println(io, "initialAsfSettings:               computation.initialAsfSettings  ")
        println(io, "intermediateConfigs:            $(computation.intermediateConfigs)  ")
        println(io, "intermediateAsfSettings:          computation.intermediateAsfSettings  ")
        println(io, "finalConfigs:                   $(computation.finalConfigs)  ")
        println(io, "finalAsfSettings:                 computation.finalAsfSettings  ")
        println(io, "alphaSettings:                    computation.alphaSettings  ")
        println(io, "einsteinSettings:                 computation.einsteinSettings  ")
        println(io, "formSettings:                     computation.formSettings  ")
        println(io, "greenSettings:                    computation.greenSettings  ")
        println(io, "hfsSettings:                      computation.hfsSettings  ")
        println(io, "isotopeSettings:                  computation.isotopeSettings  ")
        println(io, "plasmaSettings:                   computation.plasmaSettings  ")
        println(io, "polaritySettings:                 computation.polaritySettings  ")
        println(io, "yieldSettings:                    computation.yieldSettings  ")
        println(io, "zeemanSettings:                   computation.zeemanSettings  ")
        println(io, "process:                        $(computation.process)  ")
        println(io, "processSettings:                  computation.processSettings  ")
    end


    """
    `JAC.Atomic.Computation(gui::Guint; comp::Atomic.Computation=Atomic.Computation())`  ... constructor that is defined by a graphical user interface.
    """
    function Computation(gui::Guint; comp::Atomic.Computation=Atomic.Computation())
        
        if      gui == Gui
            println("aa")
            output = Atomic.ComputationGui(comp)
        elseif  gui == GuiSettings
            println("bb")
            output = Atomic.ComputationGuiSettings(comp)
        else  error("Unsupported Guint = $gui.")
        end
        
        return( output )
    end


    """
    `JAC.Atomic.ComputationGui(comp::Atomic.Computation)`  ... constructor that defines the standard part of a computation by a graphical 
         user interface apart from the setting and the configurations.
    """
    function ComputationGui(comp::Atomic.Computation)
        nmd = comp.nuclearModel
        
        title = comp.name;    nmd = comp.nuclearModel;    grid = comp.grid;    calcLevelProperties = comp.calcLevelProperties
        properties          = comp.properties
        configs             = comp.configs;                   asfSettings             = comp.asfSettings
        initialConfigs      = comp.initialConfigs;            initialAsfSettings      = comp.initialAsfSettings
        intermediateConfigs = comp.intermediateConfigs;       intermediateAsfSettings = comp.intermediateAsfSettings
        finalConfigs        = comp.finalConfigs;              finalAsfSettings        = comp.finalAsfSettings
        alphaSettings       = comp.alphaSettings;             einsteinSettings        = comp.einsteinSettings
        formSettings        = comp.formSettings;              greenSettings           = comp.greenSettings
        hfsSettings         = comp.hfsSettings;               isotopeSettings         = comp.isotopeSettings
        plasmaSettings      = comp.plasmaSettings;            polaritySettings        = comp.polaritySettings
        yieldSettings       = comp.yieldSettings;             zeemanSettings          = comp.zeemanSettings 
        process             = comp.process;                   processSettings         = comp.processSettings
        
        # Update the nuclear model
        tn1 = "Update nuclear model: "
        bn1 = slider(1:110, label = "charge", value = nmd.Z)
        bn2 = dropdown([nmd.model, "Fermi", "point", "uniform"])
        bn3 = spinbox(label="mass  ";   value=nmd.mass)
        bn4 = spinbox(label="radius";   value=nmd.radius)
        bn5 = spinbox(label="2*spin";   value=nmd.spinI.num)
        bn6 = spinbox(label="mu    ";   value=nmd.mu)
        bn7 = spinbox(label="Q     ";   value=nmd.Q)
        #
        #
        update = button("Update all")
        ui = vbox( hbox( pad(0em, tn1) ),
                   hbox( pad(1em, bn1) ), 
                   hbox( pad(1em, bn2), pad(1em, bn3), pad(1em, bn4) ), 
                   hbox( pad(1em, bn5), pad(1em, bn6), pad(1em, bn7) ), 
                   hbox( pad(2em, update) )
                  )
        #
        Interact.display(ui) 
        output = Interact.@map  (&update;     Computation(title, 
                    Nuclear.Model( observe(bn1)[], observe(bn2)[], observe(bn3)[], observe(bn4)[], 
                                   AngularJ64( Int64(observe(bn5)[])//2 ), observe(bn6)[], observe(bn7)[] ),
                    grid, calcLevelProperties, properties, configs, asfSettings, 
                    initialConfigs, initialAsfSettings, 
                    intermediateConfigs, intermediateAsfSettings, 
                    finalConfigs, finalAsfSettings, 
                    alphaSettings, einsteinSettings, formSettings, greenSettings, hfsSettings, isotopeSettings, plasmaSettings, 
                    polaritySettings, yieldSettings, zeemanSettings, 
                    process, processSettings) )
        #
        return( output )
    end


    """
    `JAC.Atomic.ComputationGuiSettings(comp::Atomic.Computation)`  ... constructor that defines the settings and configurations of a computation 
         by a graphical user interface.
    """
    function ComputationGuiSettings(comp::Atomic.Computation)
        ac = comp;   nmd = ac.nuclearModel
        
        t1 = "A slider: "
        b1 = slider(1:110, label = "charge", value = nmd.Z)
        update = button("Update")
        ui = vbox( hbox( pad(0em, t1) ),
                   hbox( pad(0em, b1) ),
                   hbox( pad(1em, update) )
                  )
        #
        Interact.display(ui) 
        output = Interact.@map  (&update; observe(b1)[] )
        return( output )
    end

end # module
