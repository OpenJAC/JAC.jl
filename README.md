

.. image:: https://travis-ci.com/OpenJAC/JAC.jl.svg?branch=master   
    :target: https://travis-ci.com/OpenJAC/JAC.jl


# Jena Atomic Calculator (JAC)* for the computation of atomic structures, processes and cascades



## What is `JAC`?

We here provide a first public version of *JAC*, the Jena Atomic Calculator and an open-source Julia package for 
doing atomic computations. JAC is a (relativistic) electronic structure code for the computation of (many-electron) 
interaction amplitudes, properties as well as a large number of excitation and decay processes for open-shell 
atoms and ions across the whole periodic table. In the future, moreover, JAC will -- more and more -- facilitate 
also studies on atomic cascades, responses as well as the time-evolution of atoms and ions. 

A primary guiding philosophy of JAC was to develop a **general and easy-to-use toolbox for the atomic physics 
community**, including an interface that is equally accessible for working spectroscopiest, theoreticians and 
code developers. Beside of its simple use, however, I also wish to provide a modern code design, a reasonable 
detailed documentation of the code and features for integrated testing. In particular, most typical calculations and 
the handling of atomic data should appear within the code similar to how they would appear in spoken or written 
language. Shortly speaking, JAC aims to provide a powerful **platform for daily use and to extent atomic theory 
towards new applications**.



## *Kinds* of computations

In some more detail, JAC distinguishes and aims to support (partly still in the future) **seven kinds of 
computations** which can be summarized as follows:

    1. **Atomic computations**, based on explicitly specified electron configurations: This kind refers to the 
       computation of level energies, atomic state representations and to either one or several atomic properties 
       for selected levels of a given multiplet. It also help compute *one* selected process at a time, if atomic 
       levels from two or more multiplets are involved in atomic transitions.
    2. **Restricted active-space computations (RAS)**: This kind concerns systematically-enlarged calculations
       of atomic states and level energies due to a specified active space of orbitals as well as due to the
       (number and/or kind of) virtual excitations that are taken to be into account. Such RAS computations are 
       normally performed stepwise by utilizing the (one-electron) orbital functions from some prior step.
    3. **Interactive computations**: Here, the (large set of) methods of the JAC program are applied interactively,
       either directly (from the REPL) or by using some short Julia script in order to compute and evaluate
       the desired observables (parameters), such as energies, expansion coefficients, transition matrices, rates,
       cross sections, etc. An interactive computation typically first prepares and applies (instances of) JAC’s
       data types, such as orbitals, configuration-state functions (CSF), bases, multiplets, and others. And like 
       Julia is built on many (high-level) functions and methods, JAC then provides the required language elements 
       for performing specific atomic computations at various degree of complexity and sophistication.
    4. **Atomic cascade computations**: A cascade typically includes ions of an element in three or more charge 
       states that are connected to each other by different atomic processes, such as photoionization, dielectronic 
       recombination, Auger decay, radiative transitions, and where the level population is determined by the 
       set-up and geometry of the given experiment. Cascade computations are usually based on some predefined 
       (cascade) approach that enables one to select automatically the state-space of the ions, the atomic 
       processes to be considered for the various steps of the cascade, and to specify perhaps additional 
       restrictions in order to keep the computations feasible.
    5. **Atomic responses**: With this kind, I wish to support in the future computations that help analyze 
       the response of atoms to incident beams of light pulses and particles, such as field-induced ionization 
       processes, high-harmonic generation and several others. For these responses, the detailed atomic 
       structure has often not yet been considered until the present but will become relevant as more elaborate 
       and accurate measurements are becoming feasible.
    6. **Atomic time-evolution of statistical tensors**: We here wish to simulate the population and coherences
       of (atomic) levels using the *Liouville equation*, when atoms and ions are irradiated by (intense) light
       pulses. For these computations, however, we always assume that the level structure of the atoms is kept 
       intact. Further (decay) processes of the excited atoms and ions can be takent into account by some loss 
       rate, but without that the atoms leave the *pre-specified space of sublevels*. In particular, I here plan 
       to consider the interaction of atoms and ions with pulses of different shape, polarization strength 
       and duration.
    7. **Semi-empirical estimates** of atomic properties, such as cross sections, stopping powers, asymptotic
       behaviour, etc. These *estimates* often refer to simple model computations or the use of fit functions.
       They are only implemented when data are needed but no *ab-initio* computations of the properties are 
       feasible otherwise.

       

## Documentation & News

A detailed `Manual, Compendium & Theoretical Background to JAC <Manual-Jac-dist.pdf>`_  is available that
describes the **use and underlying atomic theory** of the JAC code. News about recent developments of JAC
are summarized `here <NEWS.rst>`_.



## Licence & Reference

The code in this repository is distributed under the `MIT licence <LICENSE.md>`_. The associated manual 
`Manual, Compendium & Theoretical Background to JAC <Manual-Jac-dist.pdf>`_ is distributed under the Creative 
Commons Attribution 4.0 International (CC BY 4.0) license.

For reference, please, use the Computer Physics Communications publication on JAC:

* S. Fritzsche, `Computer Physics Communications xx, yy (2019); <https://doi.org/10.1016/j.cpc.2019.01.012>`_



## Dependencies and external code used in JAC

The JAC code makes use of:
* standard Julia packages, such as SpecialFunctions, FortranFiles, GaussQuadrature, GSL and QuadGK.
* ANCO library, G. Gaigalas, S. Fritzsche and I. P. Grant, Comp. Phys. Commun. 139, 269 (2001).
* Matrix elements from G. Gaigalas, T. Zalandauskas and S. Fritzsche, Comp. Phys. Commun. 157, 239 (2004).

    
    
## Quickstart

The numerous features of JAC can be easily understood by following the tutorials that are distributed together
with the code. Further details can then be found from the 
`Manual, Compendium & Theoretical Background to JAC <Manual-Jac-dist.pdf>`_. Make use the index or a
full-text search to find selected items in this (.pdf) manual.

A very first **simple example** has been discussed in the reference above and refers to the low-lying level 
structure and the Einstein A and B coefficients of the 3s 3p^6 + 3s^2 3p^4 3d → 3s^2 3p^5 transition array 
for Fe^9+ ions, also known as the spectrum Fe X. To perform such a computation within the framework of JAC, 
one needs to specify the initial- and final-state configurations in an instance of an `Atomic.Computation`, 
together with the specifier `process=RadiativeX`. We here also provide a title (line), the multipoles 
(default E1) and the gauge forms for the coupling of the radiation field that are to be applied in these 
calculations:


```julia
    comp = Atomic.Computation("Energies and Einstein coefficients for the spectrum Fe X",  Nuclear.Model(26.);
                    initialConfigs = [Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 3d")],
                    finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")], 
                    process = RadiativeX, 
                    processSettings = Radiative.Settings([E1, M1, E2, M2], [UseCoulomb, UseBabushkin] )
    perform(comp::Atomic.Computation)
```

This example is discussed also in the `tutorial <06-compute-Fe-X-spectrum.ipynb>`_.
    
    
## Tutorials

The following IJulia/jupyter notebooks introduce the reader to JAC and demonstrate various features of this toolbox.  
They can be explored statically at GitHub or can be run locally after the software repository has been cloned and installed.
In order to modify the cell-output of the notebooks and to better print the *wide tables*, you can create or modify the file
~/.jupyter/custom/custom.css in your home directory and add the line:  div.output_area pre { font-size: 7pt;} .

* `Getting started <tutorials/01-getting-started.ipynb>`_ 

* `Simple hydrogenic estimates <tutorials/02-hydrogenic-computations.ipynb>`_

* `Nuclear models and potentials <tutorials/03-setting-the-nucleus.ipynb>`_

* `Atomic potentials <tutorials/04-compare-radial-potentials.ipynb>`_

* `SCF + CI computations for carbon <tutorials/05-compute-SCF+CI-carbon-III.ipynb>`_

* `Einstein coefficients for Fe X <tutorials/06-compute-Fe-X-spectrum.ipynb>`_

* `Neon K-LL Auger rates <tutorials/07-compute-Ne-KLL-Auger-spectrum.ipynb>`_

* `Argon 2s, 2p photoionization <tutorials/09-compute-Ar-2s-2p-photoionization.ipynb>`_

* Several further tutorials are also available.



## Encouragement & Contributions

The scope of JAC is much larger than what I can (and plan to) implement myself here in Jena. 
With JAC's upload to GitHub, I therefore wish to encourage the users to fork the code and to report improvements,
failures, bugs, etc. Non-trivial changes to the code can be made available via pull requests, i.e. 
by submitting code for review (by other users) prior to their merger with the master code. 

In particular, **I like to encourage contributions from the atomic physics community** if the overall style of the 
program is maintained and if consensus exists how to add new features to the code. The goal should be to avoid 
duplication and inhomogeneity across the package as well as to implement (too) specific features that may cause 
issues in the future. External support by developers may include incremental improvements as well as multiple 
approaches for algorithms and modules in order to provide well-tested alternatives, for instance, if some particular 
approach does not work properly in all applications. Moreover, emphasis will be placed first on all those 
applications that receive enough attention by the community. 

In contrast, we shall not support those developments which appears too sophisticated or detrimental to a 
long-term maintenance of the code. Other specialized parts might be incorporated later if the code has left its 
early stage of development and becomes robust enough.

Although a good number of tests have been made on JAC, this is still a very first implementation, and no code is
error free. I shall therefore appreciate reports from the users if problems are encountered or, more helpful, 
if solutions are provided. One of the simplest way to start contributing to JAC is writing a tutorial, in addition 
to those provided above, in order to navigate others to the task of a new user. Also, new graphical user 
interface and plotting features on different outcomes of atomic computations will be very helpful for the community. 
A few further suggestions can be found in section 1.7 in the 
`Manual, Compendium & Theoretical Background to JAC <Manual-Jac-dist.pdf>`_.



## Developers:

* Stephan Fritzsche,  `s.fritzsche@gsi.de`



## Supporters:

