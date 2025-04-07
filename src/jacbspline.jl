
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
#==
include("module-BsplinesN.jl");        using ..BsplinesN
include("module-Hamiltonian.jl");      using ..Hamiltonian
include("module-SelfConsistent.jl");   using ..SelfConsistent  ==#

grid = Radial.Grid(true)
primitives = BsplinesN.generatePrimitives(grid)

#==
wbL = primitives.bsplinesL;    wbS = primitives.bsplinesS
@show BsplinesN.computeOverlap(wbL[55], wbL[55],   grid::Radial.Grid)
@show BsplinesN.computeNondiagonalD( 1, -1, wbL[2], wbL[5], grid::Radial.Grid)
@show BsplinesN.computeNondiagonalD(-1, -1, wbL[2], wbL[5], grid::Radial.Grid)


nm = Nuclear.Model(10.)
if       nm.model == "point"
    pot = Nuclear.pointNucleus(nm.Z, primitives.grid)
elseif   nm.model == "Fermi"
    pot = Nuclear.fermiDistributedNucleus(nm.radius, nm.Z, primitives.grid) 
elseif   nm.model == "uniform"
    pot = Nuclear.uniformNucleus(nm.radius, nm.Z, primitives.grid)
else
    error("stop a")
end

shells    = Basics.generateShellList(1, 2, 1)
subshells = Basics.generateSubshellList(shells)
BsplinesN.generateOrbitals(subshells, pot::Radial.Potential, nm, primitives; printout=true)  ==#


wa  = Atomic.Computation(Atomic.Computation(), name="Ne", nuclearModel=Nuclear.Model(10.00); grid=grid,
                         configs = [Configuration("[He] 2s 2p^4")] )
wb  =perform(wa; output=true)
bs1 = wb["multiplet:"].levels[1].basis

asfSettings = AsfSettings();   nm=Nuclear.Model(10.00)
bs2 = SelfConsistent.solveMeanFieldBasis(bs1, nm, primitives, asfSettings; printout=true)


configs = [Configuration("[He] 2s 2p^4")]
asfSettings = AsfSettings()

SelfConsistent.performSCF(configs, nm, grid, asfSettings, printout=true)
