
#
println("Eh) Photoabsorption cascade computations.")
using JLD
#
setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure sine")
setDefaults("unit: energy", "eV")
grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
#
if  false
    # Calculation of a photoabsorption cascade that is purely based on the photoionization of the target
    setDefaults("print summary: open", "photoabsorption-cascade.sum")
    name = "Photoabsorption of neon"
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(12.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=Cascade.PhotoAbsorptionScheme([Photo(), PhotoExc()], [E1, M1], [en for en = 25.0:25.0:300.0], 
                                                                     [Shell("2s"), Shell("2p")], Shell[Shell("3s"), Shell("3p")], LevelSelection(false),
                                                                     [0, 1, 2], 0. ),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Simulation of the (total) photoabsorption cross sections
    setDefaults("print summary: open", "photoabsorption-simulation.sum")

    data     = [JLD.load("zzz-cascade-excitation-computations-2023-03-11T21.jld")]
    name     = "Simulation of total photoabsorption CS for Ne"
    energies = [en for en = 25.0:1.0:101.0 ]
    property = Cascade.PhotoAbsorptionCS(false, energies, [(1,1.0)], Configuration[])
    settings = Cascade.SimulationSettings(false, false, 0.)
    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, 
                              property=property, settings=settings, computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
end
