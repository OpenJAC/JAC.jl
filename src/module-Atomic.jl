
"""
`module  JAC.Atomic`  ... a submodel of JAC that contains all methods to set-up and process (simple) atomic and SCF computations as well as 
                          atomic representations.
"""
module Atomic

    using Interact
    using ..Basics, ..Radial, ..ManyElectron, ..Nuclear
    using ..Einstein, ..Hfs, ..IsotopeShift, ..PlasmaShift, ..LandeZeeman, ..AlphaVariation, ..FormFactor, ..DecayYield,
          ..MultipolePolarizibility, ..PhotoEmission, ..PhotoIonization, ..PhotoExcitation, ..PhotoRecombination, 
          ..AutoIonization, ..Dielectronic, ..ImpactExcitation, ..CoulombExcitation, ..CoulombIonization,
          ..PhotoExcitationFluores, ..PhotoExcitationAutoion, ..RayleighCompton, ..MultiPhotonDeExcitation,
          ..PhotoIonizationFluores, ..PhotoIonizationAutoion, ..ImpactExcitationAutoion, ..RadiativeAuger, 
          ..MultiPhotonIonization, ..MultiPhotonDoubleIon, ..InternalConversion
          ##x , ..ElectricDipoleMoment

export  MeanFieldSettings, MeanFieldBasis, CiSettings, CiExpansion, RasSettings, RasStep, RasExpansion, 
        GreenSettings, GreenChannel, GreenExpansion, Representation

    include("../src/module-Atomic-inc-Representation.jl")
    

    
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
        ##x + greenSettings                  ::GreenFunction.Settings        ... Settings for approximate Green function calculations.
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
        ##x greenSettings                  ::GreenFunction.Settings
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
                    AlphaVariation.Settings(), Einstein.Settings(), FormFactor.Settings(), ##x GreenFunction.Settings(), 
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
         formSettings::FormFactor.Settings = FormFactor.Settings(), ##x greenSettings::GreenFunction.Settings = GreenFunction.Settings(), 
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
                         ##x greenSettings::GreenFunction.Settings = GreenFunction.Settings(), 
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
                    alphaSettings, einsteinSettings, formSettings, ##x greenSettings, 
                    hfsSettings, isotopeSettings, plasmaSettings, 
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
        ##x println(io, "greenSettings:                    computation.greenSettings  ")
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
    `Atomic.Computation(gui::Guint; comp::Atomic.Computation=Atomic.Computation())`  ... constructor that is defined by a graphical user interface.
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
    `Atomic.ComputationGui(comp::Atomic.Computation)`  ... constructor that defines the standard part of a computation by a graphical 
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
        formSettings        = comp.formSettings;              ##x greenSettings           = comp.greenSettings
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
                    alphaSettings, einsteinSettings, formSettings, ##x greenSettings, 
                    hfsSettings, isotopeSettings, plasmaSettings, 
                    polaritySettings, yieldSettings, zeemanSettings, 
                    process, processSettings) )
        #
        return( output )
    end


    """
    `Atomic.ComputationGuiSettings(comp::Atomic.Computation)`  ... constructor that defines the settings and configurations of a computation 
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
