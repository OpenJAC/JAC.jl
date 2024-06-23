#
println("Jb) Apply & test the line-shift computations.")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

if  true
    # Last successful:  15Jun2024
    # Compute the level energy shifts in a plasma for Fe X ions
    nm          = Nuclear.Model(26.0)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 10.0 * 1)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    settings    = Plasma.Settings(temperature, rho)
    scheme      = Plasma.LineShiftScheme(Basics.DebyeHueckel(), 0.25, 0., 0)
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=settings,                            
                                     refConfigs = [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")] )
    @show wa
    wb          = perform(wa, output=true)
    #
elseif  false
    # Last successful:  xxxx
    # It still need to be worked out how the line-shift rates are to be calculated.
    # Compute the 
    nm          = Nuclear.Model(6.0)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 10.0 * 1)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    settings    = Plasma.Settings(temperature, rho)
    scheme      = Plasma.LineShiftScheme(Basics.DebyeHueckel(), 0.25, 0., 0)
    # 
    wa   = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs=[Configuration("1s 2s^2 2p^6")],
                            finalConfigs  =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")], 
                            processSettings = AutoIonization.PlasmaSettings(Basics.DebyeHueckel(), 0.25, 0., 0, true, LineSelection() ) )

    wb = perform(wa)
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=settings,                            
                                     refConfigs = [Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")] )
    @show wa
    wb          = perform(wa, output=true)
end
