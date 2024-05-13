
println("Fa) Compute & test a three-step cascade model following the 1s-3p photo-excitation of Si^+.")

using JLD
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


if  true   
    # Last successful:  unknown ...
    # Cascade computation for Si^- with a 1s --> 3p excitation
    setDefaults("print summary: open", "zzz-Cascade-computation.sum")

    name = "Cascade after Si- 1s-3p excitation"
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(14.), grid=grid, approach=Cascade.AverageSCA(),
                            scheme=Cascade.StepwiseDecayScheme([Auger, Radiative], 2, 0, Shell[], Shell[]),
                            initialConfigs = [Configuration("1s^1 2s^2 2p^6 3s^2 3p^4")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
else
    # Last successful:  unknown ...
    # Compute 
    setDefaults("print summary: open", "zzz-Cascade-simulation.sum")

    ## data = [JLD.load("zzz-cascade-decay-computations-2020-01-25T11.jld"), JLD.load("zzz-cascade-decay-computations-2020-01-25T11.jld")]
    data = [JLD.load("zzz-cascade-decay-computations-2020-01-25T21.jld")]
    name = "Simulation after Si- 1s-3p excitation"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, properties=Cascade.AbstractSimulationProperty[Cascade.IonDistribution()], ## , Cascade.FinalLevelDistribution()], 
                            settings=Cascade.SimulationSettings(0., 0., 0., 0., [(1, 2.0), (2, 1.0), (3, 0.5)]), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
end
