
# Compile and print the examples (files) for developing the JAC tools  
println("The examples are ordered into several main branches: \n")
println("A)  Tests & development of the electronic structure part.")
println("B)  Tests & development of atomic amplitudes")
println("C)  Tests & development of atomic properties")
println("D)  Tests & development of basic atomic processes")
println("E)  Tests & development of composed atomic processes")
println("F)  Tests & development of atomic cascades")
println("G)  Tests & development of symbolic evaluations of Racah expressions")
println("H)  Tests & development of semiempirical computations")
println("I)  Tests & development of atomic response computations")
println("J)  Tests & development of plasma computations \n")
print(  "Enter a Letter (A, B, ...) to select the associated examples; Q ... to quit:  ")
char = read(stdin, Char) 

if      char == 'A'
    println("A)  Tests & development of the electronic structure part.")
    println("---------------------------------------------------------")
    println("Aa) Test of several radial e-e potentials.")
    println("Ab) Tests towards a SCF field: B-spline primitives and one-particle spectra in a local potential.")
    println("Ac) Test of the CI part for an internally generated neon multiplet without Breit interaction.")
    println("Ad) Test of the frequency-independent Breit interaction for an internally generated neon multiplet.")
    println("Ae) Test of the QED model corrections to the level structure of atoms and ions.")
    println("Af) Generate and normalize continuum orbitals in a local potential.")
    println("Ag) Test of the jj-LS transformation of levels from a given multiplet.")
    println("Ah) Test of a mean-field basis, mean-field multiplet and configuration-interaction (CI) expansion.")
    println("Ai) Test of restricted-active-space (RAS) expansion.")
    println("Aj) Test of the Green(function) expansion.")
    println("Ak) Test of spin-angular coefficients.")
    println("Al) Test of parallel computating with Julia.")
    #
elseif  char == 'B'
    println("B)  Tests & development of atomic amplitudes")
    println("--------------------------------------------")
    println("Ba) Tests of the dipole, em multipole and momentum-transfer amplitudes.")
    println("Bb) Tests of the parity non-conservation, Schiff moment and anapole moment amplitudes.")
    #
elseif  char == 'C'
    println("C)  Tests & development of atomic properties")
    println("---------------------------------------------")
    println("Ca) Test of the Einstein module with ASF from an internally generated multiplet.")
    println("Cb) Test of the Hfs module for HFS A,B parameters and hyperfine representation with ASF from an internally generated multiplet.")
    println("Cc) Test of the IsotopeShift module with ASF from an internally generated multiplet.")
    println("Cd) Test of the LandeZeeman module with ASF from an internally generated multiplet.")
    println("Ce) Test of the FormFactor module with ASF from an internally generated multiplet.")
    println("Cf) Test of the reduced 1- and 2-particle density matrices & natural orbitals.")
    println("Cg) Test of the PlasmaShift module with ASF from an internally generated multiplet.")
    println("Ch) Test of the RadiativeOpacity module with ASF from an internally generated multiplet.")
    println("Ci) Test of the AlphaVariation module with ASF from an internally generated multiplet.")
    println("Cj) Test of the MultipolePolarizibility module with ASF from an internally generated multiplet.")
    println("Ck) Test of the DecayYield module with ASF from an internally generated multiplet.")
    #
elseif  char == 'D'
    println("D)  Tests & development of basic atomic processes")
    println("--------------------------------------------------")
    println("Da) Test of the PhotoEmission module with ASF from an internally generated initial- and final-state multiplet.")
    println("Db) Test of the PhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dc) Test of the Photoionization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dd) Test of the PhotoRecombination module with ASF from an internally generated initial- and final-state multiplet.")
    println("De) Test of the Auger module with ASF from an internally generated initial- and final-state multiplet.")
    println("Df) Test of the Dielectronic module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Dg) Test of the RayleighCompton module with ASF from an internally generated initial- and final-state multiplets.")
    println("Dh) Test of the MultiPhotonDeExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    
    println("Di) Test of the DoubleAutoIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dj) Test of the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dk) Test of the RadiativeAuger module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dl) Test of the ImpactExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dm) Test of the InternalRecombination module with ASF from an internally generated initial and final-state multiplet.")
    println("Dn) Test of the TwoElectronOnePhoton module with ASF from an internally generated initial and final-state multiplet.")
    println("Dp) Test of the ParticleScattering module with ASF from an internally generated initial and final-state multiplet.")
    println("Dq) Test of the InternalConversion module with ASF from an internally generated initial- and final-state multiplet.")
    println("Dr) Test of the CoulombIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ds) Test of the CoulombExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    #
elseif  char == 'E'
    println("E)  Tests & development of composed atomic processes")
    println("-----------------------------------------------------")
    println("Ea) Test of the BeamPhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")
    println("Eb) Test of the PhotoExcitationFluores module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Ec) Test of the PhotoExcitationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Ed) Test of the Photoionization and PlasmaShift module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ef) Test of the Auger and PlasmaShift modules with ASF from an internally generated initial- and final-state multiplet.")
    #
    println("Eg) Test of the ImpactExcitationAutoion module with ASF from an internally generated initial-, intermediate and final-state multiplet.")
    println("Ei) Test of the MultiPhotonIonization module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ej) Test of the MultiPhotonDoubleIon module with ASF from an internally generated initial- and final-state multiplet.")
    println("Ek) Test of the PairAnnihilation1Photon.jl module with ASF from an internally generated initial- and final-state multiplet.")
    println("El) Test of the PhotoIonizationFluores module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    println("Em) Test of the PhotoIonizationAutoIon module with ASF from an internally generated initial-, intermediate and final-state multiplets.")
    #
elseif  char == 'F'
    println("F)  Tests & development of atomic cascades")
    println("------------------------------------------")
    println("Fa) Three-step cascade after 1s-3p photo-excitation of Si^+: Configuration model.")
    println("Fb) Simulation of the four-step cascade after 1s-3p photo-excitation of Si^+: Configuration model.")
    println("Fc) Three-step cascade computations and simulations for the decay of the neon 1s^-1 3p hole states: AverageSCA model.")
    println("Fd) Photoabsorption of Ne^+: AverageSCA model.")
    println("Fe) Dielectronic recombination cascade.")
    println("Ff) Radiative recombination cascade.")
    println("Fg) Expansion opacity (cascade) calculations.")
    println("Fh) Photoabsorption cascade computations.")
    println("Fi) Electron-impact excitation (cascade) computations for lithium-like Ne.")
    #
elseif  char == 'G'
    println("G)  Tests & development of symbolic evaluations of Racah expressions")
    println("--------------------------------------------------------------------")
    println("Ga) Symbolic evaluation by means of special values & recursion relations of the Wigner n-j symbols.")
    println("Gb) Symbolic simplification of Racah expressions by means of sum rules.")
    println("Gc) Symbolic simplification of recoupling coefficients by means of sum rules.")
    println("Gd) Symbolic simplification of Racah expressions including spherical harmonics and Wigner rotation matrices.")
    println("Ge) Symbolic simplification of spherical tensor operators, matrix elements and spherical amplitudes.")
    #
elseif  char == 'H'
    println("H)  Tests & development of semiempirical computations")
    println("-----------------------------------------------------")
    println("Ha) Tests of the semi-empirical computation of ADK rates.")
    println("Hb) Tests of (empirical) impact-ionization cross sections.")
    #
elseif  char == 'I'
    println("I)  Tests & development of atomic response computations")
    println("-------------------------------------------------------")
    println("Ia) Tests of the HighHarmonic module to calculate a high-harmonic spectrum in single-electron approximation.")
    println("Ib) Tests of the StrongField module to calculate energy and momentum distributions of photoelectron in SFA.")
    println("Ic) Tests of the StrongField module to demonstrate the Coulomb asymmetry in ATI azimuthal angular distributions.")
    println("Id) Tests of the StrongField2 module to compute in ATI photoelectron momentum distributions.")
    #
elseif  char == 'J'
    println("J)  Tests & development of plasma computations")
    println("----------------------------------------------")
    println("Ja) Tests of average-atom computations.")
    #
elseif  char == 'Q'
    return( nothing )
else    
    println("Unknown Letter; redo ...")
end

