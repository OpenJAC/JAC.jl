
#
#=
Example: Simulation of the four-step cascade after 1s-3p photo-excitation of Si^+: Configuration model
=#
println("E2) Simulation of the four-step cascade after 1s-3p photo-excitation of Si^+: Configuration model.")
#
using JLD
## wr = load("xxx-cascade-2018-01-20T19.jld");   data = wr["cascade data:"];
wa = Cascade.Simulation([Cascade.IonDist, Cascade.FinalDist], Cascade.ProbPropagation, Cascade.SimulationSettings() )
wb = perform(wa, data; output=false)

nothing

