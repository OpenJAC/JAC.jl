#
println("Ha) Tests of the HighHarmonic module to calculate a high-harmonic spectrum in single-electron approximation.")
#
setDefaults("print summary: open", "zzz-HHG.sum")
grid       = JAC.Radial.Grid("grid: exponential")
setDefaults("standard grid", grid::Radial.Grid; printout=true)

observable = HighHarmonic.HhgSpectrum()
approach   = HighHarmonic.HhgHydrogenicSaddlePoint()
pulse      = HighHarmonic.HhgLinearPlaneWaveLaser(2.0, intensity=3.51e16)   ## omega = 2.0 eV; intensity = 10^14 W/cm^2
timeMesh   = HighHarmonic.computeTimeMesh(2, 50)                                    ## 2 cycles with 50 time points / cycle
initialOrb = Radial.Orbital(Subshell("1s_1/2"), 0.5)   
target     = HighHarmonic.TargetCloud( HighHarmonic.LocalizedTargetModel( CartesianPoint(0., 0., 0.)), WeightedCartesian[])
detector   = HighHarmonic.DetectorScreen( HighHarmonic.LocalizedDetectorModel( CartesianPoint(0., 0., 0.)), WeightedCartesian[])

wa = HighHarmonic.Computation( observable, approach, pulse, timeMesh, Nuclear.Model(2.0), initialOrb, target, detector)

## wb = HighHarmonic.perform(wa; output=true)
## setDefaults("print summary: close", "")
