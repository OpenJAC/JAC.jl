

[![codecov](https://codecov.io/gh/OpenJAC/JAC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenJAC/JAC.jl)
[![Build Status](https://github.com/OpenJAC/JAC.jl/workflows/CI/badge.svg)](https://github.com/OpenJAC/JAC.jl/actions)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/OpenJAC/JAC.jl/master)  ... use JAC immediately with Binder in the cloud.


# Jena Atomic Calculator (JAC) for the computation of atomic representations, processes and cascades

*Last update:* April, 15th, 2025


## What is JAC?

**JAC**, the **Jena Atomic Calculator**, provides an open-source Julia package for doing atomic computations of various kind
and complexity. In particular, JAC is a (relativistic) electronic structure code for the computation of (atomic many-electron) 
interaction amplitudes, properties as well as a good number of excitation and decay processes for open-shell atoms and ions 
across the whole periodic table. In recent years, moreover, emphasis has been placed support atomic cascade computations
in different physical contexts as well as symbolic analysis (simplifications) of expressions from Racah's algebra.
Some further work is done (or planned) to incorporate central features for studying atomic -- strong-field -- responses
to external fields and particles, or the time-evolution of atoms and ions in the framework of (time-dependent) density
matrix and Liouville equations.

A primary guiding philosophy of JAC was to develop a **general and easy-to-use toolbox for the atomic physics community**, 
including an interface that is equally accessible for scientists working in astro and plasma physics, atomic spectroscopy,
theoretical physics as well as for code developers. Beside of its simple use, however, we wish to provide and support a modern
code design, a reasonable detailed documentation of the code and features for integrated testing. Indeed, the JAC toolbox 
facilitates many typical computations and the handling of atomic data by providing input interfaces similar to what one 
uses in a in spoken or written language. Shortly speaking, JAC aims to provide a powerful **platform for daily use and 
to extent atomic theory towards new applications** or, eventually, a **community platform for Just Atomic Computations**.

**Remark**: Although major efforts have been undertaken during the past eight years, JAC is still in an *early state* 
of development and includes features that have been implemented only partly or which are not yet tested well. JAC is an
open-source *physics code* (in its original sense) for which, despite of possible failures and deficiencies, we kindly ask 
potential users and developers for response, support and encouragement.




## *Kinds* of computations

In some more detail, JAC distinguishes and aims to support **different kinds of computations** which can be summarized 
as follows:

1. **Atomic computations**, based on explicitly specified electron *configurations*: This kind refers to the frequently
    applied computation of level energies and atomic state representations. It also supports the computation of either 
    (one or) several atomic properties for selected atomic levels from a given multiplet or of **one** selected process 
    at a time, if atomic levels from two or more multiplets and/or charge states are involved in some atomic transition,
    such as the photo or autoionization of atoms.
2. **Atomic representations**: This kind concerns the generation of different representations of atomic wave functions; 
    in particular, it includes systematically-enlarged restricted active-space (RAS) computations of atomic states 
    and level energies due to a pre-specified active space of orbitals as well as due to the (number and/or kind of) 
    virtual excitations to be taken to be into account. Such RAS computations are normally performed stepwise by making 
    use of the (one-electron) orbital functions from some prior step. Other atomic representations refer to approximate 
    atomic Green functions and, in the future, could be combined with concepts from close-coupling, (exterior) complex 
    scaling, DMRG or perturbation theory.
3. **Interactive computations**: In practice, the (large set of) methods of the JAC toolbox can always be applied also 
    interactively, either directly at Julia's REPL or by using some -- more or less -- short Julia scripts in order to 
    compute and evaluate the desired observables (atomic parameters), such as energies, expansion coefficients, transition
    matrices and amplitudes, rates, cross sections, etc. An interactive computation typically first prepares and applies 
    (certain instances of) JAC’s data types, such as orbitals, configuration-state functions (CSF), atomic bases, levels, 
    multiplets, and others. And like Julia, that is built upon many (high-level) functions and methods, JAC then provides 
    the required language elements for performing specific atomic computations at different degree of complexity and 
    sophistication.
4. **Atomic cascade computations**: An atomic cascade typically includes ions (of one element) in three or more 
    different charge states. These charge states are connected to each other by different atomic processes, such as 
    photo ionization, radiative and dielectronic recombination, Auger decay, the photo excitation and emission, or various 
    others. In practice, moreover, the relative level population of these charge states is usually determined by the
    specific set-up and geometry of the considered experiment. These cascade computations are usually based on some 
    predefined *(cascade) approach* that enables one to automatically select the state-space of the ions, to choose the 
    atomic processes to be considered for the various steps of the cascade, and to specify perhaps additional limitations
    in order to keep the computations feasible. In addition, these cascade computations are generally divided into two
    parts, the (cascade) computation for determining all necessary many-electron amplitudes and the (so-called)
    simulations to combine the amplitudes due to the experimental scenario of interest.
5. **Empirical Computations**: Not all atomic  properties and processes, such as the -- single and multiple -- 
    electron-impact ionization, stopping powers or tunnel ionization rates, can be efficiently described by ab-initio 
    many-body techniques. If needed, they are often easier computed by using empirical formulas and models. 
    A number of empirical computations are now supported to deal with such models; are often based on simple
    electronic structure calculations, together with empirically obtained parameters. These computations are only 
    implemented when data are needed but no *ab-initio* computations of the involved processes appears to be feasible. 
6. **Plasma Computations**: The notion of atoms and ions in plasma has been frequently applied to analyze the 
    behaviour of plasma in situations where the level structure and effective single-electron states remain partly 
    intact. Useful plasma computations refer to shifted photo lines, ionic mixtures in local and non-local 
    thermodynamic equilibria, or applications of the average-atom model.
7. **Symbolic evaluation of expressions from Racah's algebra**: This kind refers to the algebraic transformation
    and simplification of (Racah) expressions, which may generally include any number of Wigner n-j symbols 
    of different kind as well as (various integrals over) the spherical harmonics, the Wigner rotation matrices
    and the Kronecker and triangular deltas. Of course, the complexity of such *Racah expressions* increases 
    very rapidly as more Wigner symbols are involved. A symbolic evaluation of these expressions is naturally 
    based on the knowledge of a large set of special values, orthogonality relations and *sum rules* that may include 
    rules with a (multiple) summation over dummy indices, cf. the monography by Varshalovich *et al* (1988).

    
Several other *different kinds of computations* have been prepared and will support the applications of the 
JAC toolbox but come with a rather limited implementation so far.
    
8. **Atomic responses**: With this kind, we partly support computations in intense laser field; they also help 
    analyze the response of atoms to incident beams of light pulses and particles, such as field-induced 
    ionization processes, high-harmonic generation and several others. For these responses, the detailed 
    structure of atoms and ions has not been considered much until today. A partial-wave formulation of these 
    strong-field processes enables one to clearly distinguish between contributions due to the atomic target,
    the Volkov states, or the shape and phase of the incident light.
9. **Atomic time-evolution of statistical tensors**: We here wish to simulate the population and coherences
    of (atomic) levels using the *Liouville equation*, when atoms and ions are irradiated by (intense) light
    pulses. For these computations, however, we shall assume always that the level structure of the atoms is kept 
    intact. Further (decay) processes of the excited atoms and ions can be taken into account by some *loss 
    rate*, but without that the atoms can leave the *pre-specified space of sublevels*. In particular, we here 
    plan to consider the interaction of atoms and ions with pulses of different shape, polarization strength 
    and duration.

       

## Documentation & News

A detailed [User Guide, Compendium & Theoretical Background to JAC](docs/UserGuide-Jac.pdf)  is available that
describes the **use and underlying atomic theory** of the JAC code. News about recent developments of JAC
are summarized [here](NEWS.md).



## Licence & Reference

The code in this repository is distributed under the [MIT licence](LICENSE.md). The associated 
[User Guide, Compendium & Theoretical Background to JAC](docs/UserGuide-Jac.pdf) is distributed under the Creative 
Commons Attribution 4.0 International (CC BY 4.0) license.

For reference to (using) this code, please, use the Computer Physics Communications publication on JAC:

* S. Fritzsche: A fresh computational approach to atomic structures, processes and cascades 
     [Computer Physics Communications 240, 1 (2019)](https://doi.org/10.1016/j.cpc.2019.01.012)
* G. Gaigalas & S. Fritzsche: Angular coefficients for symmetry-adapted configuration states in jj-coupling.
     [Comp. Phys. Commun. 267, 108086 (2021)](https://doi.org/10.1016/j.cpc.2021.108086)
* S. Fritzsche, P. Palmeri & S. Schippers: Atomic cascade computations. [Symmetry 13, 520 (2021)](https://doi.org/10.3390/sym13030520)
* S. Fritzsche: Symbolic evaluation of expressions from Racah’s algebra. [Symmetry 13, 1558 (2021)](https://doi.org/10.3390/sym13091558)
* S. Fritzsche & A. Surzhykov: Approximate atomic Green functions. [Molecules 26, 2660 (2021)](https://doi.org/10.3390/molecules26092660)
* S. Fritzsche: Dielectronic recombination strengths and plasma rate coefficients of multiply-charged ions.
     [A&A 656, A163 (2021)](https://doi.org/10.1051/0004-6361/202141673)
* S. Fritzsche: Level structure and properties of open f-shell elements. [Atoms 10, 7 (2022)](https://doi.org/10.3390/atoms10010007)
* S. Fritzsche: Photon emission from hollow ions near surfaces. [Atoms 10, 37 (2022)](https://doi.org/10.3390/atoms10020037)
* S. Fritzsche, B. Böning: Strong-field ionization amplitudes for atomic many-electron targets. [Atoms 10, 70 (2022)](https://doi.org/10.3390/atoms10030070)
* S. Fritzsche: Application of symmetry-adapted atomic amplitudes. [Atoms 10, 127 (2022)](https://doi.org/10.3390/atoms10040127)
* S. Fritzsche, A.V. Maiorova & Z.W. Wu: Radiative recombination plasma rate coefficients of multiply-charged ions.
     [Atoms 11, 50 (2023)](https://doi.org/10.3390/atoms11030050)
* S. Fritzsche, L.G. Jiao, Y.C. Wang & J.E. Sienkiewicz; Collision strengths of astrophysical interest for multiply charged ions.
     [Atoms 11, 80 (2023)](https://doi.org/10.3390/atoms11050080)
* S. Fritzsche, A.K. Sahoo, L.Sharma, Z.W. Wu & S. Schippers; Merits of atomic cascade computations.
     [European Physical Journal D 78, 75 (2024)](https://link.springer.com/article/10.1140/epjd/s10053-024-00865-z)
* S. Fritzsche, H.K. Huang, Z.K. Huang, S. Schippers, W.Q. Wen and Z.W. Wu; Dielectronic recombination into high-n Rydberg shells.
     [European Physical Journal D 79, 22 (2025)](https://link.springer.com/article/10.1140/epjd/s10053-025-00969-0)
* S. Fritzsche; Atomic input for modeling ionic mixtures in astrophysical plasma.
     [European Physical Journal A 61, 63 (2025)](https://doi.org/10.1140/epja/s10050-025-01535-3)


See also [`CITATION.bib`](CITATION.bib) for the relevant references(s).




## Installation

In Julia, you can install the JAC package like any other package by by just entering the package manager (with <Alt> ]) 
and by typing
```
pkg> add https://github.com/OpenJAC/JAC.jl
```
Further information is given in the demo [Install Julia ...](demos/xxx.jl)


## Dependencies and external code used in JAC

The JAC code makes use of:
* standard Julia packages, such as BSplineKit, SpecialFunctions, GaussQuadrature, GSL and QuadGK.
* Matrix elements from G. Gaigalas and S. Fritzsche, Comp. Phys. Commun. 267, 108086 (2021).

    
    
## Quickstart

The 'simplest access' to the JAC toolbox is by using Binder in the cloud. If you click here:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/OpenJAC/JAC.jl/master) 

you will get a Jupyter notebook where you can call 'using JAC' in order to have Julia and JAC (completely) 
installed. -- Then you can run all examples and calls like on your own computer, just a bit slower 
(say, by a factor 3..5). This will help you to run a few first examples (as shown in the example folder above) 
and in order to decide of whether you wish to install the code locally.


You can also directly access the 
[Getting started with JAC](demos/...)
tutorial in the cloud, and similar for other tutorials that are distributed together with the code. Further details 
can then be found from the [User Guide, Compendium & Theoretical Background to JAC](docs/UserGuide-Jac.pdf). 
Make use of the index or a full-text search to find selected items in this (.pdf) User Guide.


A very **simple example** has been discussed in the [CPC reference](https://doi.org/10.1016/j.cpc.2019.01.012)
above and just refers to the low-lying level structure and the Einstein A and B coefficients of the 
3s 3p^6 + 3s^2 3p^4 3d -> 3s^2 3p^5  transition array for Fe^{9+} ions, also known as the spectrum Fe X. To perform such 
a computation within the framework of JAC, one needs to specify the initial- and final-state configurations by an instance
of an `Atomic.Computation`. We here also provide a title (line), the multipoles (default E1) and the gauge forms for the 
coupling of the radiation field that are to be applied in these calculations:


```julia
    grid = Radial.Grid(true);   setDefaults("standard grid", grid)
    defaultsSettings = PhotoEmission.Settings()
    photoSettings    = PhotoEmission.Settings(defaultsSettings, multipoles=[E1, M1], gauges=[UseCoulomb, UseBabushkin], printBefore=true)
    
    comp = Atomic.Computation(Atomic.Computation(), name="Energies and Einstein coefficients for the spectrum Fe X",  
                              grid = grid, nuclearModel = Nuclear.Model(26.);
                              initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
                              finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")], 
                              processSettings = photoSettings ); 
    perform(comp::Atomic.Computation)
```

This example is discussed also in one of the [tutorials](tutorials/51-compute-Fe-X-spectrum.ipynb) below or directly in the
[cloud](https://mybinder.org/v2/gh/openjac/JAC.jl/master?labpath=tutorials%2F51-compute-Fe-X-spectrum.ipynb).
    
    
## Demos

The following Pluto.jl notebooks introduce the reader to JAC and demonstrate several features of this toolbox. ---
They can be explored statically at GitHub or can be run locally after the software repository has been cloned and installed.

* [Getting started with Julia](demos/11-getting-started-with-Julia.ipynb)
* [Getting started with JAC](demos/12-getting-started-with-JAC.ipynb)
* [Simple estimates for hydrogenic atoms and ions](demos/13-compute-hydrogenic-orbitals.ipynb)
* [Specifying nuclear models and potentials](demos/16-define-nuclear-model-parameters.ipynb)
* [Selection and use of atomic potentials](demos/21-compare-radial-atomic-potentials.ipynb)
* [Self-Consistent-Field (and CI) computations for carbon](demos/22-compute-SCF+CI-carbon-III.ipynb)
* [Generate extended configuration lists](demos/23-generate-configuration-lists.ipynb)
* [Determine LS notation for atomic levels](demos/25-expand-levels-into-LS-basis.ipynb)
* [Estimate QED corrections for beryllium-like ions](demos/26-estimate-QED-for-beryllium-like.ipynb)
* [Compute the atomic level structure in a Debye-Hückel plasma](tutorials/28-compute-atoms-in-plasma.ipynb)
* [Generate an atomic mean field and apply it for CI computations](tutorials/31-generate-mean-field+ci-expansion.ipynb)
* [Compute transition probabilities for Fe X](tutorials/51-compute-Fe-X-spectrum.ipynb)
* [Compute the 2s, 2p photoionization of argon](tutorials/53-compute-Ar-2s-2p-photoionization.ipynb)
* [Compute the K-LL Auger rates of atomic neon](tutorials/54-compute-Ne-KLL-Auger-spectrum.ipynb)
* [Compute K-LL Auger rates in a Debye-Hückel plasma](tutorials/55-compute-Auger-rates-in-DH-plasma.ipynb)
* Several [other tutorials](tutorials/) are available, and this list will be extended with the further development of JAC.



## Current limitations of JAC

Although JAC has been designed for all atoms and ions across the periodic table, a number of limitations occur:
* All self-consistent-field computations are based so far on a local potential (e.g. core-Hartree, Kohn-Sham, 
  Dirac-Hartree-Slater, ...) that can be controlled by the user.
* Until the present, no serious optimization has been done for the code; this restricts most computations
  to CSF expansion with several hundred CSF.
* All continuum orbitals are generated in a local potential (Dirac-Hartree-Slater or others, as above) of the ionic core, 
  and without the explicit treatment of the exchange interaction.



## Encouragement & Contribution

The scope of JAC is much wider than what we can (and plan to) implement ourselves here in Jena. 
With making JAC an official Julia package (at Github), we wish to **encourage the users to fork the code and to 
make and announce improvements, to report failures, bugs, etc.** Non-trivial changes to the code can be made available 
via pull requests, i.e. by submitting code for review (by us and other users), and with the goal to merge these advancements
with the master code. However, since JAC is still a *physics code*, this merging may need some time to enable us to
understand and test for side-effects upon other parts of the code.

In particular, **we like to encourage contributions from the atomic physics community** to contribute to the code, provided
that the overall style of the package is maintained and if consensus exists how to add new features to the code. The goal 
should be to *avoid* duplication and inhomogeneity across the package as well as to implement (too) specific features 
that will cause issues in the future. External support by developers may include incremental improvements as well as 
multiple approaches for algorithms and modules in order to provide well-tested alternatives, for instance, if some 
particular approach does not work properly in all applications. Moreover, emphasis will be placed first on all those 
applications that receive enough attention by the community. --- In contrast, we shall not support developments that 
are highly sophisticated or detrimental to a long-term maintenance of the code. 

Although a good number of tests and applications have been performed by using JAC, this code still in an early stage, 
and no code is error free. We shall therefore appreciate reports from the users if problems are encountered or, 
more helpful, if solutions are suggesteg and/or provided. One of the simplest way to start contributing to JAC is 
writing some new *demos* (Pluto notebooks), in addition to those provided above, in order to navigate others to the 
task of a new user. Also, new graphical user interface and plotting features on different outcomes of atomic computations 
will be very helpful for the community. A few further suggestions for extending and improving JAC can be found in 
section 1.7 in the [User Guide, Compendium & Theoretical Background to JAC](docs/UserGuide-Jac.pdf).



## Developers:

* Fritzsche, Stephan, s.fritzsche@gsi.de (U Jena, Germany)
* Sahoo, Aloka Kumar, aloka_s@ph.iitr.ac.in (HI Jena, Germany)
* Huang, Houke (Institute of Modern Physics, Lanzhou, China)
* Wang, Wu (HI Jena, Germany & U Haikou, China)
* Li, Bowen (U Lanzhou, China)



## (Former) Supporters:

* Böning, Birger (HI Jena, Germany)
* Dar, Danish F. (U Jena, Germany)
* Gaigalas, Gediminas (U Vilnius, Lithuania)
* Gilles, Jan (PTB Braunschweig, Germany)
* Hofbrucker, Jiri (U Jena, Germany)
* Jiao, Li-Guang (HI Jena, Germany & Jilin U Changchun, China)
* Schippers, Stefan (U Giessen, Germany)
* Sienkiwicz, Joseph (TU Gdansk, Poland)
* Surzhykov, Andrey (U Braunschweig, Germany)
* Volotka, Andrey (HI Jena, Germany & St. Petersburg, Russia)
* Wang, Yuan-Cheng (HI Jena, Germany & U Shenyang, China)
* Wu, Zhongwen (HI Jena, Germany & U Lanzhou, China)
