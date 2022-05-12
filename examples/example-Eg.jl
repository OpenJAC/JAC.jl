
#
println("Eg) Expansion opacity calculations.")
using JLD
#
setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure sine")
grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
#
if  false
    # Calculation of the expansion opacity cascade for P^2+
    setDefaults("print summary: open", "opacity-cascade.sum")
    name = "Expansion opacities for P^2+"
    # asfSettings = AsfSettings(AsfSettings(), maxIterationsScf = 36)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(15.), grid=grid, approach=Cascade.SCA(),
                               scheme=Cascade.ExpansionOpacityScheme([E1], 3.0, 0., 0., 1, [Shell("3s"), Shell("3p")], [Shell("3s"), Shell("3p"), Shell("3d")]),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6 3s^2 3p")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Simulation of the expansion opacity cascade for P^2+
    setDefaults("print summary: open", "opacity-simulation.sum")

    data = [JLD.load("zzz-cascade-expansion-opacity-computations-2022-05-10T20.jld")]
    name = "Simulation of expansion opacities for P^2+"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.ExpansionOpacities(1.0, 1000., 1.0e5, 0., [1.0e3, 1.0e4, 1.0e5]), 
                              settings=Cascade.SimulationSettings(true, false, 0.), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")

end
