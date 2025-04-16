
"""
`module  JAC.Atomic`  
	... a submodel of JAC that contains all methods to set-up and process (simple) atomic and SCF computations as well as 
	    atomic representations.
"""
module Atomic


## using Interact
using   ..Basics, ..Defaults, ..LSjj, ..Radial, ..ManyElectron, ..Nuclear, ..PhotoEmission, ..SelfConsistent
export  perform

#== using  ..Einstein, ..Hfs, ..IsotopeShift, ..LandeZeeman, ..AlphaVariation, ..FormFactor, ..DecayYield,
       ..MultipolePolarizibility, ..PhotoIonization, ..PhotoExcitation, ..PhotoRecombination, 
       ..AutoIonization, ..DoubleAutoIonization, ..DielectronicRecombination, ..ImpactExcitation, ..CoulombExcitation, ..CoulombIonization,
       ..PhotoDoubleIonization, ..PhotoExcitationFluores, ..PhotoExcitationAutoion, ..RayleighCompton, ..MultiPhotonDeExcitation,
       ..PhotoIonizationFluores, ..PhotoIonizationAutoion, ..ImpactExcitationAutoion, ..RadiativeAuger, 
       ..MultiPhotonIonization, ..MultiPhotonDoubleIon, ..InternalConversion  ==#


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
    + processSettings                ::Basics.AbstractProcessSettings  ... Provides the settings for the selected process.
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
    if  processSettings         == nothing  prsx                     = comp.processSettings         else  prsx                     = processSettings          end 
    
    
    cp = Computation(namex, nuclearModelx, gridx, propertySettingsx, configsx, asfSettingsx, initialConfigsx, initialAsfSettingsx,     
                        intermediateConfigsx,  intermediateAsfSettingsx, finalConfigsx, finalAsfSettingsx,        
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
    if  comp.processSettings != Nothing   sa = sa * "for the process (comp.process) and with the \n\ninitial configurations:       "
        for  config  in  comp.initialConfigs   sa = sa * string(config) * ",  "         end
        if  length(comp.intermediateConfigs) > 0
            sa = sa * "\nintermediate configurations:  "
            for  config  in  comp.intermediateConfigs     sa = sa * string(config) * ",  "     end
        end
        sa = sa * "\nfinal configurations:         "
        for  config  in  comp.finalConfigs     sa = sa * string(config) * ",  "         end
    else                                sa = sa * "for the properties $(comp.properties) and with the \nconfigurations:        "
        for  config  in  comp.configs   sa = sa * string(config) * ",  "                end
    end
    return( sa )
end


# `Base.show(io::IO, comp::Atomic.Computation)`  ... prepares a printout of comp::Atomic.Computation.
function Base.show(io::IO, comp::Atomic.Computation)
    sa = Base.string(comp);             print(io, sa, "\n\n")
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
    ##x if  Plasma in comp.properties     &&  comp.plasmaSettings   != PlasmaShift.Settings()
    ##x     println(io, "plasmaSettings:               \n$(comp.plasmaSettings)  ")             end
    if  Polarity in comp.properties   &&  comp.polaritySettings != MultipolePolarizibility.Settings()
        println(io, "polaritySettings:             \n$(comp.polaritySettings)  ")           end
    if  Yields in comp.properties     &&  comp.yieldSettings    != DecayYield.Settings()
        println(io, "yieldSettings:                \n$(comp.yieldSettings)  ")              end
    if  Zeeman in comp.properties     &&  comp.zeemanSettings   != LandeZeeman.Settings()
        println(io, "zeemanSettings:               \n$(comp.zeemanSettings)  ")             end
    end
end

#==
"""
`Basics.perform(computation::Atomic.Computation)`  
    ... to perform the computation as prescribed by comp. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.perform(computation::Atomic.Computation; output=true)`  
    ... to perform the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(computation::Atomic.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel = computation.nuclearModel

    # Distinguish between the computation of level energies and properties and the simulation of atomic processes
    if   length(computation.configs) != 0
        ##x basis     = Basics.performSCF(computation.configs, nModel, computation.grid, computation.asfSettings)
        ##x multiplet = Basics.performCI( basis, nModel, computation.grid, computation.asfSettings)
        multiplet = SelfConsistent.performSCF(computation.configs, nModel, computation.grid, computation.asfSettings, printout=true)
        LSjj.expandLevelsIntoLS(multiplet, computation.asfSettings.jjLS)
        #
        if output    results = Base.merge( results, Dict("multiplet:" => multiplet) ) 
                        results = Base.merge( results, Dict("grid:"      => computation.grid) )  end
        
        # Now compute all requested properties
        for settings  in computation.propertySettings
            if      typeof(settings) == Einstein.Settings
                outcome = Einstein.computeLines(multiplet,        computation.grid, settings)    
                if output    results = Base.merge( results, Dict("Einstein lines:" => outcome) )                  end
                #
            ## elseif  typeof(settings) == Hfs.Settings    && settings.calcIJFexpansion  
            ##     outcome = Hfs.computeHyperfineMultiplet(multiplet, nModel, computation.grid, settings)         
            ##     if output    results = Base.merge( results, Dict("IJF multiplet:" => outcome) )                   end
            ##     #
            elseif  typeof(settings) == Hfs.Settings
                outcome = Hfs.computeOutcomes(multiplet, nModel,  computation.grid, settings)         
                if output    results = Base.merge( results, Dict("HFS outcomes:" => outcome) )                    end
                #
            elseif  typeof(settings) == LandeZeeman.Settings 
                outcome = LandeZeeman.computeOutcomes(multiplet, nModel,  computation.grid, settings)      
                if output    results = Base.merge( results, Dict("Zeeman parameter outcomes:" => outcome) )       end
                #
            elseif  typeof(settings) == StarkShift.Settings 
                outcome = StarkShift.computeOutcomes(multiplet, nModel,  computation.grid, settings)      
                if output    results = Base.merge( results, Dict("Stark-shift outcomes:" => outcome) )            end
                #
            elseif  typeof(settings) == IsotopeShift.Settings 
                outcome = IsotopeShift.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Isotope parameter outcomes:" => outcome) )      end
                #
            elseif  typeof(settings) == AlphaVariation.Settings 
                outcome = AlphaVariation.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("alpha variation parameter outcomes:" => outcome) )      end
                #
            elseif  typeof(settings) == FormFactor.Settings 
                outcome = FormFactor.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Form factor outcomes:" => outcome) )            end
                #
            elseif  typeof(settings) == DecayYield.Settings 
                outcome = DecayYield.computeOutcomes(computation.configs, computation.asfSettings, 
                                                        multiplet, nModel, computation.grid, settings)     
                if output    results = Base.merge( results, Dict("Fluorescence and AutoIonization yield outcomes:" => outcome) )   end
                #
            elseif  typeof(settings) == MultipolePolarizibility.Settings 
                outcome = MultipolePolarizibility.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Polarizibility outcomes:" => outcome) )         end
                #
            elseif  typeof(settings) == ReducedDensityMatrix.Settings 
                outcome = ReducedDensityMatrix.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("RDM outcomes:" => outcome) )         end
                #
            #xx elseif  typeof(settings) == PlasmaShift.Settings 
            #xx     outcome = PlasmaShift.computeOutcomes(multiplet, nModel, computation.grid, computation.asfSettings, settings)         
            #xx     if output    results = Base.merge( results, Dict("Plasma shift outcomes:" => outcome) )            end
            end
        end
    
    #== Evaluate processes that need special SCF and CI procedures
    elseif  typeof(computation.processSettings)  in [AutoIonization.PlasmaSettings, PhotoIonization.PlasmaSettings]
        pSettings        = computation.processSettings
        plasmaSettings   = Plasma.Settings(pSettings.plasmaModel, pSettings.lambdaDebye, pSettings.ionSphereR0, pSettings.NoBoundElectrons)
      # initialBasis     = Basics.performSCF(computation.initialConfigs, nModel, computation.grid, computation.initialAsfSettings)
        initialMultiplet = perform("computation: CI for plasma",  initialBasis, nModel, computation.grid, computation.initialAsfSettings, plasmaSettings)
      # finalBasis       = Basics.performSCF(computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
        finalMultiplet   = perform("computation: CI for plasma",  finalBasis, nModel, computation.grid, computation.finalAsfSettings, plasmaSettings)
        #
        if      typeof(computation.processSettings)  == Plasma.AugerSettings
            outcome = AutoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("AutoIonization lines in plasma:" => outcome) )         end
        elseif  typeof(computation.processSettings)  == Plasma.PhotoSettings
            outcome = PhotoIonization.computeLinesPlasma(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Photo lines in plasma:" => outcome) )                  end
        else
            error("stop a")
        end  ==#
        
    else
        ##x initialBasis     = Basics.performSCF(computation.initialConfigs, nModel, computation.grid, computation.initialAsfSettings)
        ##x initialMultiplet = Basics.performCI( initialBasis, nModel, computation.grid, computation.initialAsfSettings)
        initialMultiplet = SelfConsistent.performSCF(computation.initialConfigs, nModel, computation.grid, computation.initialAsfSettings)
        LSjj.expandLevelsIntoLS(initialMultiplet, computation.initialAsfSettings.jjLS)
        ##x finalBasis       = Basics.performSCF(computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
        ##x finalMultiplet   = Basics.performCI(   finalBasis, nModel, computation.grid, computation.finalAsfSettings)
        finalMultiplet   = SelfConsistent.performSCF(computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
        LSjj.expandLevelsIntoLS(finalMultiplet, computation.finalAsfSettings.jjLS)
        if  output   results["initialMultiplet"] = initialMultiplet;   results["finalMultiplet"] = finalMultiplet    end 
        #
        #==
        if typeof(computation.processSettings) in [PhotoExcitationFluores.Settings, PhotoExcitationAutoion.Settings, PhotoIonizationFluores.Settings, 
                                                    PhotoIonizationAutoion.Settings, ImpactExcitationAutoion.Settings, 
                                                    DielectronicRecombination.Settings, ResonantInelastic.Setting]
            ##x intermediateBasis     = Basics.performSCF(computation.intermediateConfigs, nModel, computation.grid, computation.intermediateAsfSettings)
            ##x intermediateMultiplet = Basics.performCI( intermediateBasis, nModel, computation.grid, computation.intermediateAsfSettings)
            intermediateMultiplet = SelfConsistent.performSCF(computation.intermediateConfigs, nModel, computation.grid, computation.intermediateAsfSettings)
            if  output   results["intermediateMultiplet"] = intermediateMultiplet    end 
        end  ==#
        #
        if      typeof(computation.processSettings) == AutoIonization.Settings 
            outcome = AutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("AutoIonization lines:" => outcome) )                  end
        elseif  typeof(computation.processSettings) == RayleighCompton.Settings 
            outcome = RayleighCompton.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Rayleigh-Compton lines:" => outcome) )                end
        elseif  typeof(computation.processSettings) == DoubleAutoIonization.Settings   
            outcome = DoubleAutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Double-Auger lines:" => outcome) )                    end
        elseif  typeof(computation.processSettings) == DielectronicRecombination.Settings 
            outcome = DielectronicRecombination.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("dielectronic recombination pathways:" => outcome) )   end
        elseif  typeof(computation.processSettings) == MultiPhotonDeExcitation.Settings
            outcome = MultiPhotonDeExcitation.computeLines(computation.processSettings.process,
                                                            finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon-de-excitation lines:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoIonization.Settings   
            outcome = PhotoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photoionization lines:" => outcome) )                 end
        elseif  typeof(computation.processSettings) == PhotoDoubleIonization.Settings   
            outcome = PhotoDoubleIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Single-photon double-ionization lines:" => outcome) ) end
        elseif  typeof(computation.processSettings) == PhotoExcitation.Settings
            outcome = PhotoExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation lines:" => outcome) )                end
        elseif  typeof(computation.processSettings) == PhotoExcitationAutoion.Settings  
            outcome = PhotoExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-autoionization pathways:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoExcitationFluores.Settings
            outcome = PhotoExcitationFluores.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-fluorescence pathways:" => outcome) )        end
        elseif  typeof(computation.processSettings) == PhotoEmission.Settings
            outcome = PhotoEmission.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative lines:" => outcome) )                       end
        elseif  typeof(computation.processSettings) == ResonantInelastic.Settings 
            outcome = ResonantInelastic.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                        computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("RIXS pathways:" => outcome) )   end
        elseif  typeof(computation.processSettings) == RadiativeAuger.Settings
            outcome = RadiativeAuger.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative Auger sharings:" => outcome) )              end
        elseif  typeof(computation.processSettings) == PhotoRecombination.Settings 
            outcome = PhotoRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo recombination lines:" => outcome) )             end
        elseif  typeof(computation.processSettings) == ImpactExcitation.Settings 
            outcome = ImpactExcitation.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == InternalRecombination.Settings 
            outcome = InternalRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("internal-recombination lines:" => outcome) )          end
        elseif  typeof(computation.processSettings) == TwoElectronOnePhoton.Settings 
            outcome = TwoElectronOnePhoton.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("two-electron-one-photon lines:" => outcome) )         end
        elseif  typeof(computation.processSettings) == ParticleScattering.Settings 
            outcome = ParticleScattering.computeEvents(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("particle-scattering events:" => outcome) )         end
        elseif  typeof(computation.processSettings) == BeamPhotoExcitation.Settings 
            outcome = BeamPhotoExcitation.computeOutcomes(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("beam-assisted photo-excitation:" => outcome) )         end
        elseif  typeof(computation.processSettings) == HyperfineInduced.Settings 
            outcome = HyperfineInduced.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("hyperfine-induced transitions:" => outcome) )         end
            #
            #
        elseif  typeof(computation.processSettings) == PairA1P()  
            outcome = PairAnnihilation1Photon.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("pair-annihilation-1-photon lines:" => outcome) )       end
        elseif  typeof(computation.processSettings) == Coulex()
            outcome = CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == Coulion()
            outcome = CoulombIonization.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb ionization lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == PhotoIonAuto()   
            outcome = PhotoIonizationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-ionization-autoionization pathways:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoIonFluor()  
            outcome = PhotoIonizationFluores.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-ionizatiton-fluorescence pathways:" => outcome) )        end
        elseif  typeof(computation.processSettings) == ImpactExcAuto()  
            outcome = ImpactExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation-autoionization pathways:" => outcome) )     end
        elseif  typeof(computation.processSettings) == MultiPI()   
            outcome = MultiPhotonIonization.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon single ionization:" => outcome) )        end
        elseif  typeof(computation.processSettings) == MultiPDI()   
            outcome = MultiPhotonDoubleIon.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon double ionization:" => outcome) )        end
        elseif  typeof(computation.processSettings) == InternalConv()   
            outcome = InternalConversion.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("internal conversion lines:" => outcome) )        end
        else
            error("stop b")
        end
    end
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
    end  ==#


end # module
