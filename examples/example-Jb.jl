#
println("Jb) Apply & test the line-shift computations.")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)

if  true
    # Last successful:  17Jul2024
    # Compute the Ne K-LL Auger energies and rates within a Debye-Hueckel plasma
    # Note that the Debye-Hueckel potential is not included into the Auger transition operator ... and that further tests are needed.
    nm              = Nuclear.Model(10.0)
    settings        = Plasma.Settings()
    initialConfigs  = [Configuration("1s 2s^2 2p^6")] 
    finalConfigs    = [Configuration("1s^2 2p^6")]   ## , Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2s^2 2p^4")
    lineSettings    = AutoIonization.PlasmaSettings()
    scheme          = Plasma.LineShiftScheme(Basics.DebyeHueckel(), initialConfigs, finalConfigs, lineSettings)
    
    comp            = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=nm, grid=grid, settings=Plasma.Settings() )
    @show comp
    wb              = perform(comp, output=true)
    #
elseif  false
    #
    # Last successful:  unknown ...
    # Prepare an example from Saha (2005++) to be later shown in some paper.
    #
elseif  false
    #
    # Last successful:  unknown ...
    # The role of these plasma-shift computations are not fully clear ... but have been removed from the Atomic.Computations()
    # This branch need to be worked out.
    setDefaults("print summary: open", "zzz-PlasmaShift.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=JAC.Radial.Grid(true), nuclearModel=Nuclear.Model(26.), 
                            configs=[Configuration("[Ne] 3s^2 3p^5"), Configuration("[Ne] 3s 3p^6")],
                            propertySettings=[ PlasmaShift.Settings(PlasmaShift.DebyeHueckel(), 0.25, 0., 0) ] )

    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end
