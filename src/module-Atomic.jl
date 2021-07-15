
"""
`module  JAC.Atomic`  
    ... a submodel of JAC that contains all methods to set-up and process (simple) atomic and SCF computations as well as 
        atomic representations.
"""
module Atomic

    ## using Interact
    using ..Basics, ..Radial, ..ManyElectron, ..Nuclear
    using ..Einstein, ..Hfs, ..IsotopeShift, ..PlasmaShift, ..LandeZeeman, ..AlphaVariation, ..FormFactor, ..DecayYield,
          ..MultipolePolarizibility, ..PhotoEmission, ..PhotoIonization, ..PhotoExcitation, ..PhotoRecombination, 
          ..AutoIonization, ..DoubleAutoIonization, ..Dielectronic, ..ImpactExcitation, ..CoulombExcitation, ..CoulombIonization,
          ..PhotoDoubleIonization, ..PhotoExcitationFluores, ..PhotoExcitationAutoion, ..RayleighCompton, ..MultiPhotonDeExcitation,
          ..PhotoIonizationFluores, ..PhotoIonizationAutoion, ..ImpactExcitationAutoion, ..RadiativeAuger, 
          ..MultiPhotonIonization, ..MultiPhotonDoubleIon, ..InternalConversion

    
    """
    `struct  Computation`  
        ... defines a type for defining  (the model of simple) atomic computation of a single multiplet, 
            including the SCF and CI as well as level properties and transition property calculations.

        + name                           ::String                          ... A name associated to the computation.
        + nuclearModel                   ::Nuclear.Model                   ... Model, charge and parameters of the nucleus.
        + grid                           ::Radial.Grid                     ... The radial grid to be used for the computation.
        + propertySettings               ::Array{Basics.AbstractPropertySettings,1}  ... List of atomic properties to be calculated.
        + configs                        ::Array{Configuration,1}          ... A list of non-relativistic configurations.
        + asfSettings                    ::AsfSettings                     
            ... Provides the settings for the SCF process and for the CI and QED calculations.
        + initialConfigs                 ::Array{Configuration,1}        
            ... A list of initial-state configurations for some transition property calculation, such as radiative transition, Auger, etc. 
        + initialAsfSettings             ::AsfSettings                     ... Provides the SCF and CI settings for the initial-state multiplet.
        + intermediateConfigs            ::Array{Configuration,1}          ... A list of initial-state configurations.
        + intermediateAsfSettings        ::AsfSettings                     ... Provides the SCF settings for the intermediate-state multiplet.
        + finalConfigs                   ::Array{Configuration,1}          ... A list of final-state configurations.
        + finalAsfSettings               ::AsfSettings                     ... Provides the SCF and CI settings for the final-state multiplet.
        ##x + alphaSettings                  ::AlphaVariation.Settings         ... Settings for alpha-variation parameter calculations.
        ##x + einsteinSettings               ::Einstein.Settings               ... Settings for Einstein coefficient calculations.
        ##x + formSettings                   ::FormFactor.Settings             ... Settings for atomic form factor calculations.
        ##x + hfsSettings                    ::Hfs.Settings                    ... Settings for hyperfine parameter calculations.
        ##x + isotopeSettings                ::IsotopeShift.Settings           ... Settings for isotope shift parameter calculations.
        ##x + plasmaSettings                 ::PlasmaShift.Settings            ... Settings for plasma-shift calculations.
        ##x + polaritySettings               ::MultipolePolarizibility.Settings .. Settings for polarizibility calculations.
        ##x + yieldSettings                  ::DecayYield.Settings             ... Settings for fluoresence and Auger yield calculations.
        ##x + zeemanSettings                 ::LandeZeeman.Settings            ... Settings for Lande-Zeeman coefficient calculations.
        ##x + process                        ::Basic.AbstractProcess           
        ##x     ... An (additional) process for which the properties are to be evaluated for the given initial- and final-state configurations.
        + processSettings                ::Basics.AbstracProcessSettings   ... Provides the settings for the selected process.
    """
    struct  Computation
        name                           ::String
        nuclearModel                   ::Nuclear.Model
        grid                           ::Radial.Grid
        propertySettings               ::Array{Basics.AbstractPropertySettings,1}
        configs                        ::Array{Configuration,1}
        asfSettings                    ::AsfSettings
        initialConfigs                 ::Array{Configuration,1} 
        initialAsfSettings             ::AsfSettings
        intermediateConfigs            ::Array{Configuration,1} 
        intermediateAsfSettings        ::AsfSettings
        finalConfigs                   ::Array{Configuration,1}
        finalAsfSettings               ::AsfSettings
        processSettings                ::Basics.AbstractProcessSettings
    end 


    """
    `Atomic.Computation()`  ... constructor for an 'empty' instance::Atomic.Computation.
    """
    function Computation()
        Computation("", Nuclear.Model(1.), Radial.Grid(), Basics.AbstractPropertySettings[], 
                    Configuration[], AsfSettings(),
                    Configuration[], AsfSettings(),
                    Configuration[], AsfSettings(),
                    Configuration[], AsfSettings(), PhotoEmission.Settings() )
    end

    
    """
    `Atomic.Computation(comp::Atomic.Computation;`
    
        name=..,                nuclearModel=..,            grid=..,                    configs=..,                   asfSettings=..,     
        initialConfigs=..,      initialAsfSettings=..,      intermediateConfigs=..,     intermediateAsfSettings=.., 
        finalConfigs=..,        finalAsfSettings=..,        alphaSettings=..,           einsteinSettings=.., 
        formSettings=..,        hfsSettings=..,             isotopeSettings=..,         plasmaSettings=..,
        polaritySettings=..,    yieldSettings::=..,         zeemanSettings=..,
        process=..,             processSettings=..,         printout::Bool=false)
                        
        ... constructor for modifying the given Atomic.Computation by 'overwriting' the previously selected parameters.
    """
    function Computation(comp::Atomic.Computation;
        name::Union{Nothing,String}=nothing,                                        nuclearModel::Union{Nothing,Nuclear.Model}=nothing,
        grid::Union{Nothing,Radial.Grid}=nothing,                                   propertySettings::Union{Nothing,Array{Basics.AbstractPropertySettings,1},Any}=nothing,   
        configs::Union{Nothing,Array{Configuration,1}}=nothing,                     asfSettings::Union{Nothing,AsfSettings}=nothing, 
        initialConfigs::Union{Nothing,Array{Configuration,1}}=nothing,              initialAsfSettings::Union{Nothing,AsfSettings}=nothing, 
        intermediateConfigs::Union{Nothing,Array{Configuration,1}}=nothing,         intermediateAsfSettings::Union{Nothing,AsfSettings}=nothing, 
        finalConfigs::Union{Nothing,Array{Configuration,1}}=nothing,                finalAsfSettings::Union{Nothing,AsfSettings}=nothing, 
        ##x alphaSettings::Union{Nothing,AlphaVariation.Settings}=nothing,              einsteinSettings::Union{Nothing,Einstein.Settings}=nothing, 
        ##x formSettings::Union{Nothing,FormFactor.Settings}=nothing,                   hfsSettings::Union{Nothing,Hfs.Settings}=nothing,
        ##x isotopeSettings::Union{Nothing,IsotopeShift.Settings}=nothing,              plasmaSettings::Union{Nothing,PlasmaShift.Settings}=nothing, 
        ##x polaritySettings::Union{Nothing,MultipolePolarizibility.Settings}=nothing,  yieldSettings::Union{Nothing,DecayYield.Settings}=nothing, 
        ##x zeemanSettings::Union{Nothing,LandeZeeman.Settings}=nothing, 
        ##x process::Union{Nothing,Basics.AbstractProcess}=nothing,                     
        processSettings::Union{Nothing,Any}=nothing,            
        printout::Bool=false)
        
        if  name                    == nothing  namex                    = comp.name                    else  namex                    = name                     end 
        if  nuclearModel            == nothing  nuclearModelx            = comp.nuclearModel            else  nuclearModelx            = nuclearModel             end 
        if  grid                    == nothing  gridx                    = comp.grid                    else  gridx                    = grid                     end 
        if  propertySettings        == nothing  propertySettingsx        = comp.propertySettings        else  propertySettingsx        = propertySettings         end 
        if  configs                 == nothing  configsx                 = comp.configs                 else  configsx                 = configs                  end 
        if  asfSettings             == nothing  asfSettingsx             = comp.asfSettings             else  asfSettingsx             = asfSettings              end 
        if  initialConfigs          == nothing  initialConfigsx          = comp.initialConfigs          else  initialConfigsx          = initialConfigs           end 
        if  initialAsfSettings      == nothing  initialAsfSettingsx      = comp.initialAsfSettings      else  initialAsfSettingsx      = initialAsfSettings       end 
        if  intermediateConfigs     == nothing  intermediateConfigsx     = comp.intermediateConfigs     else  intermediateConfigsx     = intermediateConfigs      end 
        if  intermediateAsfSettings == nothing  intermediateAsfSettingsx = comp.intermediateAsfSettings else  intermediateAsfSettingsx = intermediateAsfSettings  end 
        if  finalConfigs            == nothing  finalConfigsx            = comp.finalConfigs            else  finalConfigsx            = finalConfigs             end 
        if  finalAsfSettings        == nothing  finalAsfSettingsx        = comp.finalAsfSettings        else  finalAsfSettingsx        = finalAsfSettings         end 
        ##x if  alphaSettings           == nothing  alphaSettingsx           = comp.alphaSettings           else  alphaSettingsx           = alphaSettings            end 
        ##x if  einsteinSettings        == nothing  einsteinSettingsx        = comp.einsteinSettings        else  einsteinSettingsx        = einsteinSettings         end 
        ##x if  formSettings            == nothing  formSettingsx            = comp.formSettings            else  formSettingsx            = formSettings             end 
        ##x if  hfsSettings             == nothing  hfsSettingsx             = comp.hfsSettings             else  hfsSettingsx             = hfsSettings              end 
        ##x if  isotopeSettings         == nothing  isotopeSettingsx         = comp.isotopeSettings         else  isotopeSettingsx         = isotopeSettings          end 
        ##x if  plasmaSettings          == nothing  plasmaSettingsx          = comp.plasmaSettings          else  plasmaSettingsx          = plasmaSettings           end 
        ##x if  polaritySettings        == nothing  polaritySettingsx        = comp.polaritySettings        else  polaritySettingsx        = polaritySettings         end 
        ##x if  yieldSettings           == nothing  yieldSettingsx           = comp.yieldSettings           else  yieldSettingsx           = yieldSettings            end 
        ##x if  zeemanSettings          == nothing  zeemanSettingsx          = comp.zeemanSettings          else  zeemanSettingsx          = zeemanSettings           end 
        ##x if  process                 == nothing  processx                 = comp.process                 else  processx                 = process                  end 
        if  processSettings         == nothing  prsx                     = comp.processSettings         else  prsx                     = processSettings          end 
        
        #==
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
        end  ==#
        
        cp = Computation(namex, nuclearModelx, gridx, propertySettingsx, configsx, asfSettingsx, initialConfigsx, initialAsfSettingsx,     
                         intermediateConfigsx,  intermediateAsfSettingsx, finalConfigsx, finalAsfSettingsx,        
                         ##x alphaSettingsx, einsteinSettingsx , formSettingsx, hfsSettingsx, isotopeSettingsx, plasmaSettingsx,
                         ##x polaritySettingsx, yieldSettingsx, zeemanSettingsx,  processx, 
                         prsx) 
                         
        if printout  Base.show(cp)      end
        return( cp )
    end

    
    """
    `Atomic.Computation( ... example for SCF computations)`  
    
            grid     = Radial.Grid(true)
            nuclearM = Nuclear.Model(18., "Fermi")
            settings = AsfSettings(AsfSettings(), selectLevelsCI = true, selectedLevelsCI = [1,2, 4,5, 7,8], jjLS = LSjjSettings(false) )
            configs  = [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")]
            Atomic.Computation(Atomic.Computation(), name="Example", grid=grid, nuclearModel=nuclearM, configs=configs, asfSettings=settings )
        
    `Atomic.Computation( ... example for the computation of atomic properties)`  
    
            grid        = Radial.Grid(true)
            nuclearM    = Nuclear.Model(26., "Fermi", 58., 3.81, AngularJ64(5//2), 1.0, 1.0)
            hfsSettings = Hfs.Settings(true, true, false, false, false, false, false, Int64[] )
            configs     = [Configuration("[Ne] 3s"), Configuration("[Ne] 3p"), Configuration("[Ne] 3d")]
            Atomic.Computation(Atomic.Computation(), name="Example", grid=grid, nuclearModel=nuclearM, configs=configs, properties=[HFS()],
                               hfsSettings=hfsSettings )
        
    `Atomic.Computation( ... example for the computation of one atomic process)`  
    
            grid           = Radial.Grid(true)
            initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")]
            finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")] 
            photoSettings  = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1, M1], gauges=[UseCoulomb], printBefore=true)
            Atomic.Computation(Atomic.Computation(), name="Example", grid=grid, nuclearModel=nuclearM;
                              initialConfigs=initialConfigs, finalConfigs=finalConfigs, 
                              process = Radiative(), processSettings=photoSettings ); 
        
        ... These simple examples can be further improved by overwriting the corresponding parameters.
    """
    function Computation(wa::Bool)    
        Atomic.Computation()    
    end


    # `Base.string(comp::Atomic.Computation)`  ... provides a String notation for the variable comp::Atomic.Computation.
    function Base.string(comp::Atomic.Computation)
        sa = "Atomic computation:    $(comp.name) for Z = $(comp.nuclearModel.Z), "
        if  comp.processSettings != Nothing   sa = sa * "for the process (comp.process) and with the \ninitial configurations:    "
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
        if  comp.processSettings != Nothing
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
