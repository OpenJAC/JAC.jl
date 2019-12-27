
"""
`module  JAC.Atomic`  ... a submodel of JAC that contains all methods to set-up and process (simple) atomic and SCF computations as well as 
                          atomic representations.
"""
module Atomic

    ## using Interact
    using ..Basics, ..Radial, ..ManyElectron, ..Nuclear
    using ..Einstein, ..Hfs, ..IsotopeShift, ..PlasmaShift, ..LandeZeeman, ..AlphaVariation, ..FormFactor, ..DecayYield,
          ..MultipolePolarizibility, ..PhotoEmission, ..PhotoIonization, ..PhotoExcitation, ..PhotoRecombination, 
          ..AutoIonization, ..Dielectronic, ..ImpactExcitation, ..CoulombExcitation, ..CoulombIonization,
          ..PhotoExcitationFluores, ..PhotoExcitationAutoion, ..RayleighCompton, ..MultiPhotonDeExcitation,
          ..PhotoIonizationFluores, ..PhotoIonizationAutoion, ..ImpactExcitationAutoion, ..RadiativeAuger, 
          ..MultiPhotonIonization, ..MultiPhotonDoubleIon, ..InternalConversion

export  MeanFieldSettings, MeanFieldBasis, CiSettings, CiExpansion, RasSettings, RasStep, RasExpansion, 
        GreenSettings, GreenChannel, GreenExpansion, Representation

    include("../src/module-Atomic-inc-Representation.jl")
    

    
    """
    `struct  Computation`  ... defines a type for defining  (the model of simple) atomic computation of a single multiplet, 
                               including the SCF and CI as well as level properties and transition property calculations.

        + name                           ::String                        ... A name associated to the computation.
        + nuclearModel                   ::Nuclear.Model          ... Model, charge and parameters of the nucleus.
        + grid                           ::Radial.Grid                   ... The radial grid to be used for the computation.
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
                    AlphaVariation.Settings(), Einstein.Settings(), FormFactor.Settings(), Hfs.Settings(), IsotopeShift.Settings(), 
                    PlasmaShift.Settings(), MultipolePolarizibility.Settings(), DecayYield.Settings(),  LandeZeeman.Settings(), 
                    NoProcess, PhotoEmission.Settings() )
    end

    
    """
    `Atomic.Computation(comp::Atomic.Computation;`
    
        name=..,                nuclearModel=..,            grid=..,                    configs=..,             asfSettings=..,     
        initialConfigs=..,      initialAsfSettings=..,      intermediateConfigs=..,     intermediateAsfSettings=.., 
        finalConfigs=..,        finalAsfSettings=..,        alphaSettings=..,           einsteinSettings=.., 
        formSettings=..,        hfsSettings=..,             isotopeSettings=..,         plasmaSettings=..,
        polaritySettings=..,    yieldSettings::=..,         zeemanSettings=..,
        process=..,             processSettings=..,         printout::Bool=false)
                        
        ... constructor for modifying the given Atomic.Computation by 'overwriting' the previously selected parameters.
    """
    function Computation(comp::Atomic.Computation;
        name::Union{Nothing,String}=nothing,                                        nuclearModel::Union{Nothing,Nuclear.Model}=nothing,
        grid::Union{Nothing,Radial.Grid}=nothing,                                   properties::Union{Nothing,Array{AtomicLevelProperty,1}}=nothing,   
        configs::Union{Nothing,Array{Configuration,1}}=nothing,                     asfSettings::Union{Nothing,AsfSettings}=nothing, 
        initialConfigs::Union{Nothing,Array{Configuration,1}}=nothing,              initialAsfSettings::Union{Nothing,AsfSettings}=nothing, 
        intermediateConfigs::Union{Nothing,Array{Configuration,1}}=nothing,         intermediateAsfSettings::Union{Nothing,AsfSettings}=nothing, 
        finalConfigs::Union{Nothing,Array{Configuration,1}}=nothing,                finalAsfSettings::Union{Nothing,AsfSettings}=nothing, 
        alphaSettings::Union{Nothing,AlphaVariation.Settings}=nothing,              einsteinSettings::Union{Nothing,Einstein.Settings}=nothing, 
        formSettings::Union{Nothing,FormFactor.Settings}=nothing,                   hfsSettings::Union{Nothing,Hfs.Settings}=nothing,
        isotopeSettings::Union{Nothing,IsotopeShift.Settings}=nothing,              plasmaSettings::Union{Nothing,PlasmaShift.Settings}=nothing, 
        polaritySettings::Union{Nothing,MultipolePolarizibility.Settings}=nothing,  yieldSettings::Union{Nothing,DecayYield.Settings}=nothing, 
        zeemanSettings::Union{Nothing,LandeZeeman.Settings}=nothing, 
        process::Union{Nothing,Basics.AtomicProcess}=nothing,                       processSettings::Union{Nothing,Any}=nothing,            
        printout::Bool=false)
        
        if  name                    == nothing  namex                    = comp.name                    else  namex                    = name                     end 
        if  nuclearModel            == nothing  nuclearModelx            = comp.nuclearModel            else  nuclearModelx            = nuclearModel             end 
        if  grid                    == nothing  gridx                    = comp.grid                    else  gridx                    = grid                     end 
        if  properties              == nothing  propertiesx              = comp.properties              else  propertiesx              = properties               end 
        if  configs                 == nothing  configsx                 = comp.configs                 else  configsx                 = configs                  end 
        if  asfSettings             == nothing  asfSettingsx             = comp.asfSettings             else  asfSettingsx             = asfSettings              end 
        if  initialConfigs          == nothing  initialConfigsx          = comp.initialConfigs          else  initialConfigsx          = initialConfigs           end 
        if  initialAsfSettings      == nothing  initialAsfSettingsx      = comp.initialAsfSettings      else  initialAsfSettingsx      = initialAsfSettings       end 
        if  intermediateConfigs     == nothing  intermediateConfigsx     = comp.intermediateConfigs     else  intermediateConfigsx     = intermediateConfigs      end 
        if  intermediateAsfSettings == nothing  intermediateAsfSettingsx = comp.intermediateAsfSettings else  intermediateAsfSettingsx = intermediateAsfSettings  end 
        if  finalConfigs            == nothing  finalConfigsx            = comp.finalConfigs            else  finalConfigsx            = finalConfigs             end 
        if  finalAsfSettings        == nothing  finalAsfSettingsx        = comp.finalAsfSettings        else  finalAsfSettingsx        = finalAsfSettings         end 
        if  alphaSettings           == nothing  alphaSettingsx           = comp.alphaSettings           else  alphaSettingsx           = alphaSettings            end 
        if  einsteinSettings        == nothing  einsteinSettingsx        = comp.einsteinSettings        else  einsteinSettingsx        = einsteinSettings         end 
        if  formSettings            == nothing  formSettingsx            = comp.formSettings            else  formSettingsx            = formSettings             end 
        if  hfsSettings             == nothing  hfsSettingsx             = comp.hfsSettings             else  hfsSettingsx             = hfsSettings              end 
        if  isotopeSettings         == nothing  isotopeSettingsx         = comp.isotopeSettings         else  isotopeSettingsx         = isotopeSettings          end 
        if  plasmaSettings          == nothing  plasmaSettingsx          = comp.plasmaSettings          else  plasmaSettingsx          = plasmaSettings           end 
        if  polaritySettings        == nothing  polaritySettingsx        = comp.polaritySettings        else  polaritySettingsx        = polaritySettings         end 
        if  yieldSettings           == nothing  yieldSettingsx           = comp.yieldSettings           else  yieldSettingsx           = yieldSettings            end 
        if  zeemanSettings          == nothing  zeemanSettingsx          = comp.zeemanSettings          else  zeemanSettingsx          = zeemanSettings           end 
        if  process                 == nothing  processx                 = comp.process                 else  processx                 = process                  end 
        if  processSettings         == nothing  prsx                     = comp.processSettings         else  prsx                     = processSettings          end 
        
        if      processx == Basics.NoProcess                                                              prsx = PhotoEmission.Settings()
        elseif  processx == Basics.Auger            && typeof(prsx) != AutoIonization.Settings            prsx = AutoIonization.Settings()
        elseif  processx == Basics.AugerInPlasma    && typeof(prsx) != PlasmaShift.AugerSettings          prsx = PlasmaShift.AugerSettings()
        elseif  processx == Basics.Radiative        && typeof(prsx) != PhotoEmission.Settings             prsx = PhotoEmission.Settings()
        elseif  processx == Basics.PhotoExc         && typeof(prsx) != PhotoExcitation.Settings           prsx = PhotoExcitation.Settings()
        elseif  process == Basics.Photo             && typeof(prsx) != PhotoIonization.Settings           prsx = PhotoIonization.Settings()
        elseif  process == Basics.PhotoInPlasma     && typeof(prsx) != PlasmaShift.PhotoSettings          prsx = PlasmaShift.PhotoSettings()
        elseif  process == Basics.Rec               && typeof(prsx) != PhotoRecombination.Settings        prsx = PhotoRecombination.Settings()
        elseif  process == Basics.Dierec            && typeof(prsx) != Dielectronic.Settings              prsx = Dielectronic.Settings()
        elseif  process == Basics.PhotoExcFluor     && typeof(prsx) != PhotoExcitationFluores.Settings    prsx = PhotoExcitationFluores.Settings()
        elseif  process == Basics.PhotoExcAuto      && typeof(prsx) != PhotoExcitationAutoion.Settings    prsx = PhotoExcitationAutoion.Settings()
        elseif  process == Basics.PhotoIonFluor     && typeof(prsx) != PhotoIonizationFluores.Settings    prsx = PhotoIonizationFluores.Settings()
        elseif  process == Basics.PhotoIonAuto      && typeof(prsx) != PhotoIonizationAutoion.Settings    prsx = PhotoIonizationAutoion.Settings()
        elseif  process == Basics.MultiPhotonDE     && typeof(prsx) != MultiPhotonDeExcitation.Settings   prsx = MultiPhotonDeExcitation.Settings()
        elseif  process == Basics.Compton           && typeof(prsx) != RayleighCompton.Settings           prsx = RayleighCompton.Settings()
        elseif  process == Basics.Eimex             && typeof(prsx) != ImpactExcitation.Settings          prsx = ImpactExcitation.Settings()
        elseif  process == Basics.ImpactExcAuto     && typeof(prsx) != ImpactExcitationAutoion.Settings   prsx = ImpactExcitationAutoion.Settings()
        elseif  process == Basics.RAuger            && typeof(prsx) != RadiativeAuger.Settings            prsx = RadiativeAuger.Settings()
        elseif  process == Basics.MultiPI           && typeof(prsx) != MultiPhotonIonization.Settings     prsx = MultiPhotonIonization.Settings()
        elseif  process == Basics.MultiPDI          && typeof(prsx) != MultiPhotonDoubleIon.Settings      prsx = MultiPhotonDoubleIon.Settings()
        elseif  process == Basics.InternalConv      && typeof(prsx) != InternalConversion.Settings        prsx = InternalConversion.Settings()
        elseif  process == Basics.Coulex            && typeof(prsx) != CoulombExcitation.Settings         prsx = CoulombExcitation.Settings()
        elseif  process == Basics.Coulion           && typeof(prsx) != CoulombIonization.Settings         prsx = CoulombIonization.Settings()
        end
        
        cp = Computation(namex, nuclearModelx, gridx, propertiesx, configsx, asfSettingsx, initialConfigsx, initialAsfSettingsx,     
                         intermediateConfigsx,  intermediateAsfSettingsx, finalConfigsx, finalAsfSettingsx,        
                         alphaSettingsx, einsteinSettingsx , formSettingsx, hfsSettingsx, isotopeSettingsx, plasmaSettingsx,
                         polaritySettingsx, yieldSettingsx, zeemanSettingsx,  processx, prsx) 
                         
        if printout  Base.show(cp)      end
        return( cp )
    end


    # `Base.string(comp::Atomic.Computation)`  ... provides a String notation for the variable comp::Atomic.Computation.
    function Base.string(comp::Atomic.Computation)
        sa = "Atomic computation:    $(comp.name) for Z = $(comp.nuclearModel.Z), "
        if  comp.process != NoProcess   sa = sa * "for the process $(comp.process) and with the \ninitial configurations:    "
            for  config  in  comp.initialConfigs   sa = sa * string(config) * ",  "         end
            if  length(comp.intermediateConfigs) > 0
                sa = sa * "\nintermediate configurations:"
                for  config  in  comp.finalConfigs     sa = sa * string(config) * ",  "     end
            end
            sa = sa * "\nfinal configurations:      "
            for  config  in  comp.finalConfigs     sa = sa * string(config) * ",  "         end
        else                                sa = sa * "for the properties $(comp.properties) and with the \nconfigurations:        "
            for  config  in  comp.configs   sa = sa * string(config) * ",  "                end
        end
        return( sa )
    end


    # `Base.show(io::IO, comp::Atomic.Computation)`  ... prepares a printout of comp::Atomic.Computation.
    function Base.show(io::IO, comp::Atomic.Computation)
        sa = Base.string(comp);             print(io, sa, "\n")
        println(io, "nuclearModel:          $(comp.nuclearModel)  ")
        println(io, "grid:                  $(comp.grid)  ")
        #
        # For the computation of some given atomic process
        if  comp.process != NoProcess  
        println(io, "processSettings:              \n$(comp.processSettings)  ")
        if  comp.initialAsfSettings != AsfSettings()     
        println(io, "initialAsfSettings:           \n$(comp.initialAsfSettings)  ")         end
        if  comp.intermediateAsfSettings != AsfSettings()     
        println(io, "intermediateAsfSettings:      \n$(comp.intermediateAsfSettings)  ")    end
        if  comp.finalAsfSettings != AsfSettings()     
        println(io, "finalAsfSettings:             \n$(comp.finalAsfSettings)  ")           end
        #
        else
        # For the computation of no or several atomic properties
        if  comp.asfSettings != AsfSettings()     
        println(io, "asfSettings:                  \n$(comp.asfSettings)  ")                end
        if  AlphaX in comp.properties     &&  comp.alphaSettings    != AlphaVariation.Settings()
        println(io, "alphaSettings:                \n$(comp.alphaSettings)  ")              end
        if  EinsteinX in comp.properties  &&  comp.einsteinSettings != Einstein.Settings()
        println(io, "einsteinSettings:             \n$(comp.einsteinSettings)  ")           end
        if  FormF in comp.properties      &&  comp.formSettings     != FormFactor.Settings()
        println(io, "formSettings:                 \n$(comp.formSettings)  ")               end
        if  HFS in comp.properties        &&  comp.hfsSettings      != Hfs.Settings()
        println(io, "hfsSettings:                  \n$(comp.hfsSettings)  ")                end
        if  Isotope in comp.properties    &&  comp.isotopeSettings  != IsotopeShift.Settings()
        println(io, "isotopeSettings:              \n$(comp.isotopeSettings)  ")            end
        if  Plasma in comp.properties     &&  comp.plasmaSettings   != PlasmaShift.Settings()
        println(io, "plasmaSettings:               \n$(comp.plasmaSettings)  ")             end
        if  Polarity in comp.properties   &&  comp.polaritySettings != MultipolePolarizibility.Settings()
        println(io, "polaritySettings:             \n$(comp.polaritySettings)  ")           end
        if  Yields in comp.properties     &&  comp.yieldSettings    != DecayYield.Settings()
        println(io, "yieldSettings:                \n$(comp.yieldSettings)  ")              end
        if  Zeeman in comp.properties     &&  comp.zeemanSettings   != LandeZeeman.Settings()
        println(io, "zeemanSettings:               \n$(comp.zeemanSettings)  ")             end
        end
    end

end # module
