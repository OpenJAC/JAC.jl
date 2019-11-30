

[![Build Status](https://travis-ci.com/OpenJAC/JAC.jl.svg?branch=master)](https://travis-ci.com/OpenJAC/JAC.jl)
[![codecov](https://codecov.io/gh/OpenJAC/JAC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenJAC/JAC.jl)



# Jena Atomic Calculator (JAC) for the computation of atomic representations, processes and cascades


## What is JAC?

We here provide a first public version of **JAC**, the **Jena Atomic Calculator** and an open-source Julia package for 
doing atomic computations. JAC is a (relativistic) electronic structure code for the computation of (atomic many-electron) 
interaction amplitudes, properties as well as a large number of excitation and decay processes for open-shell 
atoms and ions across the whole periodic table. In forthcoming years, moreover, JAC will -- more and more -- facilitate 
also studies on atomic cascades, responses to external fields and particles, the time-evolution of atoms and ions 
as well as selected symbolic computations of expressions from Racah's algebra. 

A primary guiding philosophy of JAC was to develop a **general and easy-to-use toolbox for the atomic physics 
community**, including an interface that is equally accessible for working spectroscopiests, theoreticians and 
code developers. Beside of its simple use, however, I also wish to provide a modern code design, a reasonable 
detailed documentation of the code as well as features for integrated testing. In particular, many typical computations
and the handling of atomic data should appear within the code similar to how they would appear in spoken or written 
language. Shortly speaking, JAC aims to provide a powerful **platform for daily use and to extent atomic theory 
towards new applications** or, in short, a **community platform for Just Atomic Computations**.

**Remark**: Although major efforts have been undertaken during the past two years, JAC is still in a very early state 
of its development and includes features that are only partly implemented or not yet tested in all detail. 
Despite of possible failures and deficiencies of the present code, however, I here annouce JAC and kindly ask potential
users and developers for response, support and encouragement.




## *Kinds* of computations

In some more detail, JAC distinguishes and aims to support (partly still within the future) **nine kinds of 
computations** which can be summarized as follows [(Figure)](Overview-Jac.pdf):

1. **Atomic computations**, based on explicitly specified electron *configurations*: This kind refers to the 
    computation of level energies, atomic state representations and to either one or several atomic properties 
    for selected levels from a given multiplet. It also help compute **one** selected process at a time, if atomic 
    levels from two or more multiplets are involved in some atomic transition.
2. **Atomic representations**: This kind concerns different representations of atomic wave functions; in particular,
    it includes systematically-enlarged restricted active-space (RAS) computations of atomic states and level 
    energies due to a pre-specified active space of orbitals as well as due to the (number and/or kind of) virtual 
    excitations that are taken to be into account. Such RAS computations are normally performed stepwise by making 
    use of the (one-electron) orbital functions from some prior step. Other atomic representations refer 
    to approximate atomic Green functions and, in the future, combined techniques with concepts from close-coupling, 
    (exterior) complex scaling, DMRG or perturbation theory.
3. **Interactive computations**: Here, the (large set of) methods of the JAC program are applied interactively,
    either directly from the REPL or by using some short Julia script in order to compute and evaluate
    the desired observables (atomic parameters), such as energies, expansion coefficients, transition matrices
    and amplitudes, rates, cross sections, etc. An interactive computation typically first prepares and applies 
    (certain instances of) JAC’s data types, such as orbitals, configuration-state functions (CSF), atomic 
    bases, levels, multiplets, and others. And like Julia, that is built on many (high-level) functions and 
    methods, JAC then provides the required language elements for performing specific atomic computations 
    at different degree of complexity and sophistication.
4. **Atomic cascade computations**: A cascade typically includes ions of an element in three or more charge 
    states that are connected to each other by different atomic processes, such as photoionization, dielectronic 
    recombination, Auger decay, radiative transitions, and where the relative level population of these charge
    states is determined by the set-up and geometry of the given experiment. Cascade computations are usually 
    based on some predefined *(cascade) approach* that enables one to automatically select the state-space of 
    the ions, to choose the atomic processes to be considered for the various steps of the cascade, and to 
    specify perhaps additional     restrictions in order to keep the computations feasible.
5. **Atomic responses**: With this kind, I wish to support in the future computations that help analyze 
    the response of atoms to incident beams of light pulses and particles, such as field-induced ionization 
    processes, high-harmonic generation and several others. For these responses, the detailed structure of 
    the atoms and ions has often not yet been considered until today but will become relevant as more elaborate 
    and accurate measurements will become feasible.
6. **Atomic time-evolution of statistical tensors**: We here wish to simulate the population and coherences
    of (atomic) levels using the *Liouville equation*, when atoms and ions are irradiated by (intense) light
    pulses. For these computations, however, we always assume that the level structure of the atoms is kept 
    intact. Further (decay) processes of the excited atoms and ions can be takent into account by some *loss 
    rate*, but without that the atoms can leave the *pre-specified space of sublevels*. In particular, I here plan 
    to consider the interaction of atoms and ions with pulses of different shape, polarization strength 
    and duration.
7. **Atomic descriptors**: In building models for machine learning, the definition of proper (atomic) descriptors 
    or *feature transformations* is central for predicting the physical properties and behaviour of atoms and ions; 
    here, we shall provide various of such descriptors useful for atomic physics.
8. **Semi-empirical estimates** of atomic properties, such as cross sections, stopping powers, asymptotic
    behaviour, etc. These *estimates* often refer to simple model computations or the use of fit functions.
    They are only implemented when data are needed but no *ab-initio* computations of the properties are 
    feasible otherwise.
9. **Symbolic evaluation of expressions from Racah's algebra**: This kind refers to the algebraic transformation
    and simplification of (Racah) expressions, which may generally include any number of Wigner n-j symbols 
    of different kind as well as (various integrals over) the spherical harmonics, the Wigner rotation matrices
    and the Kronecker and triangular deltas. Of course, the complexity of such *Racah expressions* increases 
    very rapidly as more Wigner symbols are involved. A symbolic evaluation of these expressions is naturally 
    based on the knowledge of a large set of special values, orthogonality relations and *sum rules* 
    that may include rules with a (multiple) summation over dummy indices, cf. the monography by 
    Varshalovich *et al* (1988).


       

## Documentation & News

A detailed [User Guide, Compendium & Theoretical Background to JAC](UserGuide-Jac.pdf)  is available that
describes the **use and underlying atomic theory** of the JAC code. News about recent developments of JAC
are summarized [here](NEWS.md).



## Licence & Reference

The code in this repository is distributed under the [MIT licence](LICENSE.md). The associated 
[User Guide, Compendium & Theoretical Background to JAC](UserGuide-Jac.pdf) is distributed under the Creative 
Commons Attribution 4.0 International (CC BY 4.0) license.

For reference to (using) this code, please, use the Computer Physics Communications publication on JAC:

* S. Fritzsche, [Computer Physics Communications 240, 1 (2019)](https://doi.org/10.1016/j.cpc.2019.01.012)



## Installation

For the installation of the JAC package, you need first to install **gfortran** to compile the Fortran code
of the angular coeffients. Afterwards, you can install JAC by just entering the package manager (with <Alt> ]) 
and by typing
```
pkg> add https://github.com/OpenJAC/JAC.jl
```
Moreover, to run the tutorials, yo need to have Jupyter notebooks running.

If problems occur with the installation, please, see the [Installation Guide](InstallationGuide.md).



## Dependencies and external code used in JAC

The JAC code makes use of:
* standard Julia packages, such as SpecialFunctions, FortranFiles, GaussQuadrature, GSL and QuadGK.
* ANCO library, G. Gaigalas, S. Fritzsche and I. P. Grant, Comp. Phys. Commun. 139, 269 (2001).
* Matrix elements from G. Gaigalas, T. Zalandauskas and S. Fritzsche, Comp. Phys. Commun. 157, 239 (2004).

    
    
## Quickstart

The numerous features of JAC can be easily understood by (first) following the tutorials that are distributed 
together with the code. Further details can then be found from the 
[User Guide, Compendium & Theoretical Background to JAC](UserGuide-Jac.pdf). Make use the index or a
full-text search to find selected items in this (.pdf) User Guide.

A very **simple example** has been discussed in the [CPC reference](https://doi.org/10.1016/j.cpc.2019.01.012)
above and just refers to the low-lying level structure and the Einstein A and B coefficients of the 
3s 3p^6 + 3s^2 3p^4 3d -> 3s^2 3p^5  transition array for Fe^{9+} ions, also known as the spectrum Fe X. 
To perform such a computation within the framework of JAC, one needs to specify the initial- and final-state 
configurations by an instance of an `Atomic.Computation`, together with the specifier `process=Radiative`. 
We here also provide a title (line), the multipoles (default E1) and the gauge forms for the coupling of the 
radiation field that are to be applied in these calculations:


```julia
    comp = Atomic.Computation("Energies and Einstein coefficients for the spectrum Fe X",  Nuclear.Model(26.);
                    initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
                    finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")], 
                    process = Radiative, 
                    processSettings = Radiative.Settings([E1, M1, E2, M2], [UseCoulomb, UseBabushkin] )
    perform(comp::Atomic.Computation)
```

This example is discussed also in one of the [tutorials](tutorials/09-compute-Fe-X-spectrum.ipynb) below.
    
    
## Tutorials

The following IJulia/jupyter notebooks introduce the reader to JAC and demonstrate several features of this toolbox.  
They can be explored statically at GitHub or can be run locally after the software repository has been cloned and installed.
In order to modify the cell-output of the notebooks and to better print *wide tables*, you can create or modify the file
~/.jupyter/custom/custom.css in your home directory and add the line:  div.output_area pre { font-size: 7pt;} .

* [Getting started](tutorials/01-getting-started.ipynb)
* [Simple estimates for hydrogenic atoms and ions](tutorials/02-obtain-hydrogenic-orbitals.ipynb)
* [Specifying nuclear models and potentials](tutorials/03-define-nuclear-model-parameters.ipynb)
* [Selection and use of atomic potentials](tutorials/04-compare-radial-potentials.ipynb)
* [Self-Consistent-Field (and CI) computations for carbon](tutorials/05-compute-SCF+CI-carbon-III.ipynb)
* [Determine LS notation for atomic levels](tutorials/06-expand-levels-into-LS-basis.ipynb)
* [Estimate QED corrections for beryllium-like ions](tutorials/07-estimate-QED-for-beryllium-like.ipynb)
* [Compute the atomic level structure in a Debye-Hückel plasma](tutorials/08-compute-atoms-in-DH-plasma.ipynb)
* [Compute transition probabilities for Fe X](tutorials/09-compute-Fe-X-spectrum.ipynb)
* [Compute the K-LL Auger rates of atomic neon](tutorials/10-compute-Ne-KLL-Auger-spectrum.ipynb)
* [Compute K-LL Auger rates in a Debye-Hückel plasma](tutorials/11-compute-Auger-rates-in-DH-plasma.ipynb)
* [Compute the 2s, 2p photoionization of argon](tutorials/12-compute-Ar-2s-2p-photoionization.ipynb)
* Several [other tutorials](tutorials/) are available, and this list will be extended with the further development of JAC.



## Current limitations of JAC

Although JAC has been designed for all atoms and ions across the periodic table, a number of limitations
occur:
* All self-consistent-field computations are based on a local potential (e.g. core-Hartree, Kohn-Sham, 
  Dirac-Hartree-Slater, ...) that can be controlled by the user.



## Encouragement & Contribution

The scope of JAC is much wider than what I can (and plan to) implement myself here in Jena. 
With JAC's upload to GitHub, I therefore wish to **encourage the users to fork the code and to report improvements,
failures, bugs, etc.** Non-trivial changes to the code can be made available via pull requests, i.e. 
by submitting code for review (by other users) prior to their merger with the master code. 

In particular, **I like to encourage contributions from the atomic physics community** if the overall style of the 
package is maintained and if consensus exists how to add new features to the code. The goal should be to *avoid*
duplication and inhomogeneity across the package as well as to implement (too) specific features that may cause 
issues in the future. External support by developers may include incremental improvements as well as multiple 
approaches for algorithms and modules in order to provide well-tested alternatives, for instance, if some particular 
approach does not work properly in all applications. Moreover, emphasis will be placed first on all those 
applications that receive enough attention by the community. 

In contrast, I shall not support those developments which appears too sophisticated or detrimental to a 
long-term maintenance of the code. Other specialized parts might be incorporated later if the code has left its 
early stage of development and becomes robust enough.

Although a good number of tests have been made on JAC, this is still a very first implementation, and no code is
error free. I shall therefore appreciate reports from the users if problems are encountered or, more helpful, 
if solutions are provided. One of the simplest way to start contributing to JAC is writing a tutorial, in addition 
to those provided above, in order to navigate others to the task of a new user. Also, new graphical user 
interface and plotting features on different outcomes of atomic computations will be very helpful for the community. 
A few further suggestions for extending and improving JAC can be found in section 1.7 in the 
[User Guide, Compendium & Theoretical Background to JAC](UserGuide-Jac.pdf).



## Developers:

* Stephan Fritzsche,  s.fritzsche@gsi.de (U Jena, Germany)



## Supporters:

