
"""
`module  JAC.Atomic`  ... a submodel of JAC that contains all methods to set-up and process (simple) atomic and SCF computations as well as 
                          atomic cascade computations; it is using JAC, JAC.Radial, JAC.ManyElectron, JAC.Nuclear, JAC.Einstein, 
                          JAC.Hfs, JAC.IsotopeShift, JAC.PlasmaShift, JAC.LandeZeeman, JAC.PhotoEmission, JAC.PhotoIonization, JAC.AutoIonization.
"""
module Atomic

    using Interact
    using JAC.Basics, JAC.Radial, JAC.ManyElectron, JAC.Nuclear
    using JAC.Einstein, JAC.Hfs, JAC.IsotopeShift, JAC.PlasmaShift, JAC.LandeZeeman, JAC.AlphaVariation, JAC.FormFactor, JAC.DecayYield,
          JAC.GreenFunction,  JAC.MultipolePolarizibility,
          JAC.PhotoEmission, JAC.PhotoIonization, JAC.PhotoExcitation, JAC.PhotoRecombination, 
          JAC.AutoIonization, JAC.Dielectronic, JAC.ImpactExcitation, JAC.CoulombExcitation, JAC.CoulombIonization,
          JAC.PhotoExcitationFluores, JAC.PhotoExcitationAutoion, JAC.RayleighCompton, JAC.MultiPhotonDeExcitation,
          JAC.PhotoIonizationFluores, JAC.PhotoIonizationAutoion, JAC.ImpactExcitationAutoion, JAC.RadiativeAuger, 
          JAC.MultiPhotonIonization, JAC.MultiPhotonDoubleIon, JAC.InternalConversion
          ##x , JAC.ElectricDipoleMoment

export  RasSettings, RasStep, RasComputation

    """
    `struct  Atomic.RasSettings`  
        ... a struct for defining the settings for a restricted active-space computations.

        + levelsScf            ::Array{Int64,1}     ... Levels on which the optimization need to be carried out.
        + maxIterationsScf     ::Int64              ... maximum number of SCF iterations in each RAS step.
        + accuracyScf          ::Float64            ... convergence criterion for the SCF field.
        
    	+ breitCI              ::Bool               ... logical flag to include Breit interactions.
    	+ selectLevelsCI       ::Bool               ... true, if specific level (number)s have been selected.
    	+ selectedLevelsCI     ::Array{Int64,1}     ... Level number that have been selected.
    """
    struct  RasSettings
        levelsScf              ::Array{Int64,1}
        maxIterationsScf       ::Int64  
        accuracyScf            ::Float64 
    	breitCI                ::Bool 
    	selectLevelsCI         ::Bool 
    	selectedLevelsCI       ::Array{Int64,1} 
     end


    """
    `Atomic.RasSettings()`  ... constructor for setting the default values.
    """
    function RasSettings()
    	RasSettings(Int64[1], 24, 1.0e-6, false, false, Int64[])
    end
    
    
    # `Base.show(io::IO, settings::RasSettings)`  ... prepares a proper printout of the settings::RasSettings.
    function Base.show(io::IO, settings::RasSettings)
    	  println(io, "levelsScf:            $(settings.levelsScf)  ")
    	  println(io, "maxIterationsScf:     $(settings.maxIterationsScf)  ")
    	  println(io, "accuracyScf:          $(settings.accuracyScf)  ")
    	  println(io, "breitCI:              $(settings.breitCI)  ")
    	  println(io, "selectLevelsCI:       $(settings.selectLevelsCI)  ")
    	  println(io, "selectedLevelsCI:     $(settings.selectedLevelsCI)  ")
    end


    """
    `struct  Atomic.RasStep`  
        ... specifies an individual step of a (relativistic) restricted active space computation for a set of levels. This struct 
            comprises all information to generate the orbital basis and to perform the associated SCF and multiplet computations for a 
            selected number of levels.

        + seFrom            ::Array{Shell,1}        ... Single-excitations from shells   [sh_1, sh_2, ...]
        + seTo              ::Array{Shell,1}        ... Single-excitations to shells  [sh_1, sh_2, ...]
        + deFrom            ::Array{Shell,1}        ... Double-excitations from shells   [sh_1, sh_2, ...]
        + deTo              ::Array{Shell,1}        ... Double-excitations to shells  [sh_1, sh_2, ...]
        + teFrom            ::Array{Shell,1}        ... Triple-excitations from shells   [sh_1, sh_2, ...]
        + teTo              ::Array{Shell,1}        ... Triple-excitations to shells  [sh_1, sh_2, ...]
        + qeFrom            ::Array{Shell,1}        ... Quadrupole-excitations from shells   [sh_1, sh_2, ...]
        + QeTo              ::Array{Shell,1}        ... Quadrupole-excitations to shells  [sh_1, sh_2, ...]
        + frozenShells      ::Array{Shell,1}        ... List of shells that are kept 'frozen' in this step.
        + constraints       ::Array{String,1}       ... List of Strings to define 'constraints/restrictions' to the generated CSF basis.
    """
    struct  RasStep
        seFrom              ::Array{Shell,1}
        seTo                ::Array{Shell,1}
        deFrom              ::Array{Shell,1}
        deTo                ::Array{Shell,1}
        teFrom              ::Array{Shell,1}
        teTo                ::Array{Shell,1}
        qeFrom              ::Array{Shell,1}
        qeTo                ::Array{Shell,1}
        frozenShells        ::Array{Shell,1}
        constraints         ::Array{String,1}
    end

    """
    `JAC.Atomic.RasStep(b::Bool)`  ... constructor for an 'empty' instance of a variable::Atomic.RasStep for b = true/false
    """
    function RasStep(b::Bool)
        RasStep(Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[],    Shell[], String[])
    end


    """
    `JAC.Atomic.RasStep(;seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[], 
                         deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[], 
                         teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[], 
                         qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[], 
                         frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[], refine::RasStep=RasStep(true)`  
        ... constructor for specifying all excitations, frozen shells and constraints optionally.
    """
    function RasStep(;seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[], 
                      deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[], 
                      teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[], 
                      qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[], 
                      frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[], refine::RasStep=RasStep(true))
        if  seFrom == Shell[]   sxFrom = refine.seFrom   else      sxFrom = seFrom      end
        if  seTo   == Shell[]   sxTo   = refine.seTo     else      sxTo   = seTo        end
        if  deFrom == Shell[]   dxFrom = refine.deFrom   else      dxFrom = deFrom      end
        if  deTo   == Shell[]   dxTo   = refine.deTo     else      dxTo   = deTo        end
        if  teFrom == Shell[]   txFrom = refine.teFrom   else      txFrom = teFrom      end
        if  teTo   == Shell[]   txTo   = refine.teTo     else      txTo   = teTo        end
        if  qeFrom == Shell[]   qxFrom = refine.qeFrom   else      qxFrom = qeFrom      end
        if  qeTo   == Shell[]   qxTo   = refine.qeTo     else      qxTo   = qeTo        end
        if  frozen == Shell[]        frozx  = refine.frozenShells   else      frozx  = frozen             end
        if  constraints ==String[]   consx  = refine.constraints    else      consx  = constraints        end
        
        RasStep( sxFrom, sxTo, dxFrom, dxTo, txFrom, txTo, qxFrom, qxTo, frozx, consx)
    end


    # `Base.string(step::Atomic.RasStep)`  ... provides a String notation for the variable step::Atomic.RasStep.
    function Base.string(step::Atomic.RasStep)
        sa = "\nRAS step with $(length(step.frozenShells)) frozen shell(s): $(step.frozenShells)  ... and virtual excitations"
        return( sa )
    end


    # `Base.show(io::IO, step::Atomic.RasStep)`  ... prepares a proper printout of the (individual step of computations) step::Atomic.RasStep.
    function Base.show(io::IO, step::Atomic.RasStep)
        sa = Base.string(step);                 print(io, sa, "\n")
        if  length(step.seFrom) > 0
           sa = "   Singles from:          { ";      for  sh in step.seFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.seTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
        if  length(step.deFrom) > 0
           sa = "   Doubles from:          { ";      for  sh in step.deFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.deTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
        if  length(step.teFrom) > 0
           sa = "   Triples from:          { ";      for  sh in step.teFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.teTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
        if  length(step.qeFrom) > 0
           sa = "   Quadruples from:       { ";      for  sh in step.qeFrom   sa = sa * string(sh) * ", "  end
           sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.qeTo     sa = sa * string(sh) * ", "  end;   
           sa = sa[1:end-2] * " }";            print(io, sa, "\n")
        end   
    end



    """
    `struct  Atomic.RasComputation`  
        ... a struct for defining a series of restricted active space computations for a specified set of levels. 
            This type provides all information to generate the atomic basis and to perform the corresponding SCF and 
            multiplet computations for a specified number of levels.

        + name             ::String                    ... to assign a name to the given model.
        + nuclearModel     ::Nuclear.Model             ... Model, charge and parameters of the nucleus.
        + grid             ::Radial.Grid               ... The radial grid to be used for the computation.
        + refConfigs       ::Array{Configuration,1}    ... List of references configurations, at least 1.
        + symmetry         ::LevelSymmetry             ... Symmetry of the levels/CSF in the many-electron basis.
        + NoElectrons      ::Int64                     ... Number of electrons.
        + NoIterations     ::Int64                     ... Number of SCf iterations to be applied in this step of computations.
        + steps            ::Array{Atomic.RasStep,1}   ... List of SCF steps that are to be done in this model computation.
        + settings         ::Atomic.RasSettings        ... Settings for the given RAS computation
    """
    struct  RasComputation
        name               ::String  
        nuclearModel       ::Nuclear.Model
        grid               ::Radial.Grid
        refConfigs         ::Array{Configuration,1} 
        symmetry           ::LevelSymmetry
        NoElectrons        ::Int64
        NoIterations       ::Int64 
        steps              ::Array{Atomic.RasStep,1}
        settings           ::Atomic.RasSettings 
    end


    """
    `JAC.Atomic.RasComputation()`  ... constructor for an 'empty' instance of the a variable::Atomic.RasComputation
    """
    function RasComputation()
        RasComputation("", Nuclear.Model(1.0), Radial.Grid(), ManyElectron.Configuration[], Basics.LevelSymmetry(0, Basics.plus),
                       0, 0, Atomic.RasStep[], Atomic.RasSettings())
    end


    # `Base.string(comp::RasComputation)`  ... provides a String notation for the variable comp::RasComputation.
    function Base.string(comp::RasComputation)
        sa = "RAS computation: '$(comp.name)' for symmétry $(comp.symmetry) and with $(length(comp.steps)) steps as well as with reference configurations: \n   "
        for  refConfig  in  comp.refConfigs     sa = sa * string(refConfig) * ",  "     end
        sa = sa * "\n" * "and the current settings:"
        return( sa )
    end


    # `Base.show(io::IO, comp::RasComputation)`  ... prepares a proper printout of the (individual step of computations) comp::RasComputation.
    function Base.show(io::IO, comp::RasComputation)
        sa = Base.string(comp);            print(io, sa, "\n")
        println(io, "NoElectrons:          $(comp.NoElectrons)  ")
        println(io, "steps:                $(comp.steps)  ")
        println(io, "settings:             $(comp.settings)  ")
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
        process                        ::AtomicProcess
        processSettings                ::Union{PhotoEmission.Settings, AutoIonization.Settings, PlasmaShift.AugerSettings, 
                                               PhotoIonization.Settings, PlasmaShift.PhotoSettings,
                                               PhotoRecombination.Settings, Dielectronic.Settings, ImpactExcitation.Settings,
                                               CoulombExcitation.Settings, CoulombIonization.Settings, 
                                               PhotoExcitation.Settings, PhotoExcitationFluores.Settings, 
                                               PhotoExcitationAutoion.Settings, RayleighCompton.Settings, 
                                               MultiPhotonDeExcitation.Settings, PhotoIonizationFluores.Settings, 
                                               PhotoIonizationAutoion.Settings, ImpactExcitationAutoion.Settings,
                                               RadiativeAuger.Settings, MultiPhotonIonization.Settings,
                                               MultiPhotonDoubleIon.Settings, InternalConversion.Settings} #= , 
                                               #
                                               PairAnnihilation1Photon.Settings } =#
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
    function Computation(sa::String, nm::Nuclear.Model; grid::Radial.Grid = Radial.Grid("grid: exponential"), ##x calc::Bool = false,
                         properties::Array{AtomicLevelProperty,1} = [Basics.NoProperty], configs::Array{Configuration,1} = Configuration[],
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
                         process::Basics.AtomicProcess = Basics.NoProcess, processSettings::Any = true)
        if      process == Basics.NoProcess         && typeof(processSettings) == Bool && processSettings             procSettings = PhotoEmission.Settings()
        elseif  process == Basics.Auger             && typeof(processSettings) == Bool && processSettings             procSettings = AutoIonization.Settings()
        elseif  process == Basics.Auger             && typeof(processSettings) == AutoIonization.Settings             procSettings = processSettings 
        elseif  process == Basics.AugerInPlasma     && typeof(processSettings) == Bool && processSettings             procSettings = PlasmaShift.AugerSettings()
        elseif  process == Basics.AugerInPlasma     && typeof(processSettings) == PlasmaShift.AugerSettings           procSettings = processSettings 
        elseif  process == Basics.Radiative         && typeof(processSettings) == Bool && processSettings             procSettings = PhotoEmission.Settings()
        elseif  process == Basics.Radiative         && typeof(processSettings) == PhotoEmission.Settings              procSettings = processSettings 
        elseif  process == Basics.PhotoExc          && typeof(processSettings) == Bool && processSettings             procSettings = PhotoExcitation.Settings()
        elseif  process == Basics.PhotoExc          && typeof(processSettings) == PhotoExcitation.Settings            procSettings = processSettings 
        elseif  process == Basics.Photo             && typeof(processSettings) == Bool && processSettings             procSettings = PhotoIonization.Settings()
        elseif  process == Basics.Photo             && typeof(processSettings) == PhotoIonization.Settings            procSettings = processSettings 
        elseif  process == Basics.PhotoInPlasma     && typeof(processSettings) == Bool && processSettings             procSettings = PlasmaShift.PhotoSettings()
        elseif  process == Basics.PhotoInPlasma     && typeof(processSettings) == PlasmaShift.PhotoSettings           procSettings = processSettings 
        elseif  process == Basics.Rec               && typeof(processSettings) == Bool && processSettings             procSettings = PhotoRecombination.Settings()
        elseif  process == Basics.Rec               && typeof(processSettings) == PhotoRecombination.Settings         procSettings = processSettings 
        elseif  process == Basics.Dierec            && typeof(processSettings) == Bool && processSettings             procSettings = Dielectronic.Settings()
        elseif  process == Basics.Dierec            && typeof(processSettings) == Dielectronic.Settings               procSettings = processSettings 
        elseif  process == Basics.PhotoExcFluor     && typeof(processSettings) == Bool && processSettings             procSettings = PhotoExcitationFluores.Settings()
        elseif  process == Basics.PhotoExcFluor     && typeof(processSettings) == PhotoExcitationFluores.Settings     procSettings = processSettings 
        elseif  process == Basics.PhotoExcAuto      && typeof(processSettings) == Bool && processSettings             procSettings = PhotoExcitationAutoion.Settings()
        elseif  process == Basics.PhotoExcAuto      && typeof(processSettings) == PhotoExcitationAutoion.Settings     procSettings = processSettings 
        elseif  process == Basics.PhotoIonFluor     && typeof(processSettings) == Bool && processSettings             procSettings = PhotoIonizationFluores.Settings()
        elseif  process == Basics.PhotoIonFluor     && typeof(processSettings) == PhotoIonizationFluores.Settings     procSettings = processSettings 
        elseif  process == Basics.PhotoIonAuto      && typeof(processSettings) == Bool && processSettings             procSettings = PhotoIonizationAutoion.Settings()
        elseif  process == Basics.PhotoIonAuto      && typeof(processSettings) == PhotoIonizationAutoion.Settings     procSettings = processSettings 
        elseif  process == Basics.MultiPhotonDE     && typeof(processSettings) == Bool && processSettings             procSettings = MultiPhotonDeExcitation.Settings()
        elseif  process == Basics.MultiPhotonDE     && typeof(processSettings) == MultiPhotonDeExcitation.Settings    procSettings = processSettings 
        elseif  process == Basics.Compton           && typeof(processSettings) == Bool && processSettings             procSettings = RayleighCompton.Settings()
        elseif  process == Basics.Compton           && typeof(processSettings) == RayleighCompton.Settings            procSettings = processSettings 
        elseif  process == Basics.Eimex             && typeof(processSettings) == Bool && processSettings             procSettings = ImpactExcitation.Settings()
        elseif  process == Basics.Eimex             && typeof(processSettings) == ImpactExcitation.Settings           procSettings = processSettings 
        elseif  process == Basics.ImpactExcAuto     && typeof(processSettings) == Bool && processSettings             procSettings = ImpactExcitationAutoion.Settings()
        elseif  process == Basics.ImpactExcAuto     && typeof(processSettings) == ImpactExcitationAutoion.Settings    procSettings = processSettings 
        elseif  process == Basics.RAuger            && typeof(processSettings) == Bool && processSettings             procSettings = RadiativeAuger.Settings()
        elseif  process == Basics.RAuger            && typeof(processSettings) == RadiativeAuger.Settings             procSettings = processSettings 
        elseif  process == Basics.MultiPI           && typeof(processSettings) == Bool && processSettings             procSettings = MultiPhotonIonization.Settings()
        elseif  process == Basics.MultiPI           && typeof(processSettings) == MultiPhotonIonization.Settings      procSettings = processSettings 
        elseif  process == Basics.MultiPDI          && typeof(processSettings) == Bool && processSettings             procSettings = MultiPhotonDoubleIon.Settings()
        elseif  process == Basics.MultiPDI          && typeof(processSettings) == MultiPhotonDoubleIon.Settings       procSettings = processSettings 
        elseif  process == Basics.InternalConv      && typeof(processSettings) == Bool && processSettings             procSettings = InternalConversion.Settings()
        elseif  process == Basics.InternalConv      && typeof(processSettings) == InternalConversion.Settings         procSettings = processSettings 
        elseif  process == Basics.Coulex            && typeof(processSettings) == Bool && processSettings             procSettings = CoulombExcitation.Settings()
        elseif  process == Basics.Coulex            && typeof(processSettings) == CoulombExcitation.Settings          procSettings = processSettings 
        elseif  process == Basics.Coulion           && typeof(processSettings) == Bool && processSettings             procSettings = CoulombIonization.Settings()
        elseif  process == Basics.Coulion           && typeof(processSettings) == CoulombIonization.Settings          procSettings = processSettings 
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
