
#
println("Dl) Apply & test the ImpactExcitation module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-ImpactExcitation.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, Alok")           ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("continuum: SCF potential", Basics.DFSField(0.50))  ## use GBL_CONT_POTENTIAL  ... to access this SCF potential internally, please

grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 2.5e-2, rbox = 30.0)
## nm = nucModel = Nuclear.Model(12., "Fermi")

if  false
    # Last successful:  xxMay2024
    # Calculation of EIE rate coefficients for ...
    rateSettings = ImpactExcitation.RateSettings(20, 6, [10^4, 10^5, 10^6.0, 10^7, 10^7.3, 10^7.5, 10^8])
    settings = ImpactExcitation.Settings(ImpactExcitation.Settings(); lineSelection = LineSelection(true, indexPairs=[(1,2)]), 
                                         calcRateCoefficient = true, rateSettings = rateSettings)

    wa = Atomic.Computation(Atomic.Computation(), name="impact excitation", grid=grid, nuclearModel=Nuclear.Model(12.), 
                                                  initialConfigs  = [Configuration("1s^2 2s^2 2p^6 3s")],
                                                  finalConfigs    = [Configuration("1s^2 2s^2 2p^6 3p")], 
                                                  processSettings = settings)
    @show wa

    perform(wa, output=true)
    #
elseif false
    # Last successful:  xxMay2024
    # Calculation of EIE rate coefficients for ...
    # ... please, improve this branch and tell briefly how Julia/JAC need to be invoked in order to make use
    #     of such serial computations, if any.
    # Defining a Fermi distributed nuclear model
    nucModel  = Nuclear.Model(12., "Fermi")
    basis     = Basics.performSCF([Configuration("1s^2 2s^2 2p^6 3s"), Configuration("1s^2 2s^2 2p^6 3p")], nucModel, grid, AsfSettings())
    multiplet = Basics.performCI(basis, nucModel, grid, AsfSettings());
    # Settings for Rate coefficient / effective collision strength calculations
    temperatures = [5e4, 1e5, 2e5, 5e5, 1e6, 2e6]
    rateSettings = ImpactExcitation.RateSettings(30, 8, temperatures)
    settings     = ImpactExcitation.Settings(ImpactExcitation.Settings(); lineSelection = LineSelection(true, indexPairs=[(1,2),(1,3)]),
                                                                          calcRateCoefficient = true, rateSettings = rateSettings )
    lines, rates = ImpactExcitation.computeLines(multiplet, multiplet, nucModel, grid, settings;  output=true)
    #
elseif false
    # Last successful:  xxMay2024
    # Calculation of EIE rate coefficients for ...
    # ... please, improve this branch and tell briefly how Julia/JAC need to be invoked in order to make use
    #     of such distributed computations.
    using Distributed
    @everywhere using JAC
    # Atomic structure calculations
    grid      = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, hp = 5e-2, rbox = 30.0)
    nucModel  = Nuclear.Model(12., "Fermi")
    basis     = Basics.performSCF([Configuration("1s^2 2s^2 2p^6 3s"), Configuration("1s^2 2s^2 2p^6 3p")], nucModel, grid, AsfSettings())
    multiplet = Basics.performCI(basis, nucModel, grid, AsfSettings());
    #
    @everywhere begin
    # Copying data to all worker processes
    grid      = $(grid)
    nucModel  = $(nucModel)
    multiplet = $(multiplet)
    # Cross section / rate coefficient / effective collision strength calculations
    temperatures = [5e4, 1e5, 2e5, 5e5, 1e6, 2e6]
    rateSettings = ImpactExcitation.RateSettings(30, 8, temperatures)
    settings = ImpactExcitation.Settings(ImpactExcitation.Settings(); lineSelection = LineSelection(true, indexPairs=[(1,2),(1,3)]),
                                                                      calcRateCoefficient = true, rateSettings = rateSettings)
    end
    #
    lines, rates = ImpactExcitation.computeLines(multiplet, multiplet,nucModel, grid, settings; output=true);
    #
end
#
setDefaults("print summary: close", "")


