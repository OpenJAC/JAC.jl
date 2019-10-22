
"""
`module  JAC.Details`  
    ... a submodel to collect Doc string information about further details about the JAC program.
"""
module Details

    """
    **`Atomic amplitudes`**

    Furthermore, JAC supports the computation of the following (transition) amplitudes:

    + Auger.amplitude()                 ...... Amplitude for the Auger electron emission of an electron.
    + FormFactor.amplitude()            ...... Amplitude for the momentum transfer q.
    + Hfs.amplitude()                   ...... Amplitude for the hyperfine interaction with the magnetic-dipole or
                                               electric-quadrupole field of the nucleus.
    + LandeZeeman.amplitude()           ...... Amplitude for the interaction with an external magnetic field.
    + ParityNonConservation.amplitude() ...... PNC, anapole and EDM amplitudes.
    + PhotoIonization.amplitude()       ...... Photoionization amplitudes for the absorption of a photon and the release of a 
                                               free electron.
    + PhotoRecombination.amplitude()    ...... Photorecombination amplitudes for the emission of a photon due to the capture of a 
                                               free electron.
    + PhotoEmission.amplitude()         ...... Transition amplitude for the emission or absorption of an multipole photon.

    """
    function amplitudes()  end




    """
    #### `Atomic computations, based on explicitly specified electron configurations`

        > computation = Atomic.Computation(..)      ... specify the computation, either by Atomic.Computation(
                                                        "interactive") or by the use of explicit constructors.
        > perform(computation)                      ... run the computation and obtain the requested results.

    + A typical computation, that is based on explicitly specified electron configurations, refers to level 
        energies, atomic states and to either one or several atomic properties of a given multiplet or to just 
        one selected process or amplitude; cf. the supported properties, amplitudes and processes below.
    + Such a computation is always based on a set of explicitly specified (non-relativistic) configurations, 
        to describe either a single multiplet of interest, or the initial, intermediate and final multiplets for 
        describing excitation and decay processes.
    + Although all configurations need to be specified explicitly, the size of these computations may grow 
        rapidly if open-shell configurations are involved.


    #### `Complete active-space computations (CAS)`   (not yet properly implemented)

        > computation = Atomic.Cas(..)              ... specify the computation, either by Atomic.Cas("interactive") 
                                                        or by the use of explicit constructors.
        > perform(computation)                      ... to run the CAS computation and to obtain all specified results.

    + A CAS computation refers to the systematically-extended generation of atomic states and level energies due to a specified 
        active space of orbitals as well as the (number and/or kind of the) virtual excitations to be concluded.
    + A CAS computations is normally performed stepwise by using the self-consistent field from some prior step.
    + If a set of systematically enlarged atomic states (ASF) is generated, all atomic properties and processes below 
        can be still computed with these ASF by (so-called) `interactive computations.`


    #### `Interactive computations`

    + In an interactive computation, the methods of the JAC program are applied interactively, either within REPL
        or by just a short script, in order to compute energies, expansion coefficients, transition matrices, rates, 
        cross sections, etc. with explicitly specified entities from the JAC program such as orbitals, bases, multiplets and others.
    + All methods from the JAC module and its submodules can be utilized, although special methods are often 
        available to facilitate the computations. Like for all Julia methods, the JAC methods are here just 
        (high-level) language elements for performing computations at various levels of sophistication.
    + This type of computation has not yet been fully implemented; further features will be added and extended; 
        see below.


    #### `Atomic cascade computations`

        > computation = Atomic.Cascade(..)    ... specify the computation, either by Atomic.Cascade("interactive") 
                                                  or by the use of explicit constructors.
        > perform(computation)                ... run the cascade computation, i.e. to determine all specifiec 
                                                  transition amplitudes between levels from two or more charge states.
        > evaluate(data::CascadeData)         ... evaluate the set of lines and amplitudes due to different criteria 
                                                  in order to obtain information about the cascade, such as ion yield,
                                                  branching fractions, major decay pathes, etc.
                                                    
    + A cascade computation typically refers to three or more charge states of an atom that are connected by various processes, 
        such as photoionization, dielectronic recombination, Auger decay, radiative transitions, etc.
    + The particular atomic processes that are taken into account for the inidividual steps of the cascade must be specified 
        explicitly.
    + Different (cascade) approaches have been predefined in order to deal with atomic cascades, for instance, within a 
        configuration-average approach, by including the fine-structure of all levels of interest or by incorporating even shake-up 
        and shake-off processes.


    #### `Atomic responses` (not yet properly implemented)

    + We here hope for further support from outside-users.                                 


    #### `Time evolution of statistical tensors in (intense) light pusles` (not yet implemented)

        > pulse1 = Pulse.ExperimentalCharacterization(..)      ... specify one or several light pulses for the 
                                                                   interaction of the atom with the radiation field.
        > adapt(pulse1, pulse2, ...)                           ... to adapt the internal representation of the light 
                                                                   pulses to a form that is suitable for the time 
                                                                   evolution of the statistical tensors of the atom.
        > frame = Atomic.LevelFrame(..)                        ... specify the level frame for the time evolution, 
                                                                   either by Atomic.LevelFrame("interactive") or by 
                                                                   the use of explicit constructors.
        
        > evolve(tensors::Array{Statistical.Tensor,1}, ..)     ... perform the time evolution of the statistical 
        > evaluate(tensors::Array{Statistical.Tensor,1}, ..)       tensor due to the given information.
                                                    
    + A time evolution of statistical tensors always proceeds within a pre-specified set of sublevels; all 
        (decay) processes that lead the system out of these sublevels must be treated by loss rates.
    + Although such a time evolution can deal with pulses of different shape, strength and length, it is
        assumed that they are `weak enough` to not disturb the level structure and level sequence substantially, 
        i.e. that every sublevel can still be characterized by its given (total) energy and symmetry.
    + At present, however, very little is yet prepared to deal with pulses of different shape or such 
        time evolutions.
    + We here hope for further support from outside-users.                                 


    #### `Semi-empirical estimates of cross sections, etc.` (not yet properly implemented)

        > estimation = Semiempirical.Estimation(..)       ... specify some properties and parameters that are to 
                                                              be estimated.
        > perform(estimation, ...)                        ... perform this estimation explicitly, often by 
                                                              providing some additional information and parameters.                  
                    

    + A semi-empirical `estimation' refers to simple model computation or the evaluation of some fit functions 
        in order to provide atomic data that cannot be generated so easily by ab-initio computations. Estimations 
        are typically built on -- more or less -- sophisticated models and external parameter optimzations.
    + Three simple examples refer to the binding energy of an inner-shell electron, which is taken from some 
        tabulation, the total or shell-dependent electron-impact ionization cross sections as well as ...
    + Estimations are obtained without the explicit computation of matrix elements or quantum amplitudes and 
        are typically derived much faster than any ab-initio computation.


    #### `Symbolic evaluation of expressions from Racah's algebra` (not yet properly implemented)

    + We here hope for further support from outside-users.                                 


    """
    function kindsOfComputation()  end




    """
    ### `Atomic cascade computations and approximations`

    There has been much recent interest in studying and analyzing the multiple electron emission from atoms and ions following their inner-shell
    excitation or ionization. Such multiple Auger emission studies have been performed for both, noble gases after inner-shell excitation and 
    ionization as well as for complex ions, such as Fe^+ (Schippers et al. 2017), Fe^+, Cd (Andersson et al. 2016), and others. Because of the 
    rapidly increasing complexity in dealing with multiple open shell, these atoms and ions can often not be treated at the same level of 
    sophistication like single autoionization processes. Therefore, a proper hierarchy of (cascade) approaches is desirable that can be utilized 
    in order to adapt the simulations to the different requests from the side of experiments. 

    In general, the atomic processes that can be included into a cascade computations are the same as listed above. In practice, however,
    the allowed atomic processes for cascade simulations only include: 
    a) Radiative transitions; b) Auger transitions; c) Photoionization lines.

    To deal with atomic excitation and decay cascades, JAC supports **four approaches**, i.e. levels of approximations:


    ### `Average singe-configuration approach (averageSCA)`  (not yet implemented)

    + This approach makes use of only one 'common set of orbitals' for all ionization stages of the atom and of *configuration-averaged*
        data throughout the simulations; at present, it simply applies the pre-defined radial orbitals of Bunge & Bunge (1993) or others 
        to set-up all necessary wavefunctions.
    + Here, only approximate configuration averaged and decay rates are applied by simply averaging over all fine-structure levels, rates and
        cross sections.
    + For dealing with continuum electrons, we use approximate orbitals without exchange and without any relaxation of the electron density 
        and by only using the first 3-4 (non-vanishing) partial waves in the computation of the photo- and autoionization amplitudes.
    + This averageSCA approach is naturally restricted also to only electric-dipole allowed (E1) transitions.
    + For complex configurations, a maximum of typical 3 total J symmetries are utilized in order to estimate the averaged rates and 
        cross sections.
    + All further improvements (with regard to this approach) CAN then occur only in a more sophisticated but also much more expensive approach.
    + **Properties** that can be considered in this averagedSCA approach: 1) probability flux between configurations; 2) relative ion yields; ...


    ### `Singe-configuration approach (SCA)`  (not yet implemented)

    + This approach applies either a common set of orbitals for all ionization stages of the atom; i.e. the pre-defined radial orbitals 
        of Bunge & Bunge (1993) or an individually-optimized orbital set for each configuration.
    + All fine-structure (many-electron) transitions amplitudes are calculated explicitly for the selected atomic processes.
    + The exchange interaction of the emitted electrons with the bound-state density and the use of different approximations for dealing with
        the continuum orbitals can be controlled to some extent; the orbital relaxation between different orbitals can be incorporated into
        the amplitudes as far this is supported by the JAC program.
    + Different multipole components can be selected for the interaction with the radiation field.
    + This single-configuration approach (SCA) is expected to already provide a reasonable description of all strong decay pathes.
    + The display and discussion of all results should be configuration 'and' level-wise; various tools for displaying the results are
        desirable and need to be developed in the future.
    + **Properties** that can be considered in this SCA approach: 1) probability flux between configurations; 2) relative ion yields; 
        3) list of configurations with corresponding level (splitting); 4) analysis of fluorescence lines; ...


    ### `Multiple-configuration approach (MCA)`  (not yet implemented)

    + When compared with the SCA approach, the central idea of the MCA approach is to incorporate also electron-electron correlations by 
        configuration mixing contributions with closeby-lying configurations and, hence, to 'group' the configurations together by using
        physical insight into the cascade process.
    + This separation of the various electron configurations of each stage of ionization into different groups is typically based on the 
        energy, a maximum size (No. of CSF) of any individual groups or on some prior knowledge about some strongly-interacting configurations; 
        this approach is often adapted manually. A proper size and balance of the configuration groups is typically essential in order to keep 
        the computations feasible.
    + It is assumed in this approach, that each physical level (uniquely) belongs to just 'one' group so that *double counting* of 
        rates, etc. canotn occur. If configurations are first grouped manually, this has to be made sure by the user.
    + A full (average-level) MCDF approach can be considered by this approximation if all levels of every ionization stage are just grouped
        alltogether; in practice, however, such a computation is likely unfeasible and of little interest.
    + Analogue to the SCA approach, energies, amplitudes, etc. are calculated for all level pairs within a given configuration group 
        or between two or more of such groups; these data are first kept in a CacscadeData structure before any evaluation of the data is made.
    + A MCA approach might be feasible until 1s-excited argon and, perhaps even, krypton but, likely, not much beyond.
    + **Properties** that can be considered in this MCA approach: 1) relative ion yields; ...

    ### `Multiple-configuration-shake approach (shakeMCA)`  (not yet implemented)

    + The shakeMCA approach is similar to the MSC approach but also includes shake-down and shake-up configurations that do not arise within
        a typically 'configuration picture'. These shake configurations can be treated as indidivual groups of configurations or together with
        configurations of the decay cascade.
    + A systematic approach to the incorporation of shake processes is supported by selecting a (we--defined) number of shake-displacements, 
        either for the whole cascade or for some individual step in the (auto-) ionization procedure.
    + Standard and conjugated shake-processes can be distinguished and selected independently of each other. Once a configuration list with
        *shaked configuations included* has been generated, all other computational steps and control is the same as for the MSC approach.
    + While there ist still one set of orbitals for each group of configurations, the electron relaxation is planned to be included by
        a biorthogonal transformation in order to allow proper non-vanishing shake amplitudes


    ### `Comments:`

    + For each of these four (cascade) approaches above, the mean computational effort should increase by (at least) one order of 
        magnitude, perhaps also by several orders.
    + Although the overall computational effort is difficult to estimate and formalize for arbitrary atomic systems, these four approaches 
        should be clearly discernible from each other and the *costs* of all simpler approaches should be typically negligble, when compared 
        with the **last and most sophisticated approach** that is still feasibleto be carried out. 
    + In practice, the averaged SCA should be always realistic for all atoms and ions, and independent of how well it describes the underlying 
        physics of the autoionization cascade.
    + The MSC and shakeMCA approaches enable one to include electron-electron correlations to a limited extent by choosing carefully
        different *groups* of configurations but by neglecting all admixtures from other groups. Nevertheless, all the presently-designed
        approaches neglect almost all correlations within the continuum.

    + **Comments for the implememtation:** A list of all generated electron configurations, whose levels are taken into accout within a
        given approach, should be made available, ordered by their averaged energies with regard to some specified *reference* energy.
    + In such a graphical representation of a cascade, the ionization stage of the atom should form the *columns* of either a table of 
        graph; for large cascade, the whole table/graph should be easily 'reducible' to some energy interval and/or a selected range of
        ionization stages.  Typically, such a list/graph will provide the user with a first overview about te cascade, enables one 
        to *zoom in* and it will likely form the basis for more sophisticated graphs  (which also allow access to *transitions* between 
        fine-structure levels of these configuration).
    + An automatically generated list of electron configurations forms the *framework* for the averageSCA and SCA approaches, and is the
        starting point to 'group' configurations together in the MCA and shakeMCA approaches. In these latter approaches, these configurations
        are often manipulated manually.
    + It need to be easily possible to add or remove electron configurations *by hand* in order to reduce the size of the configuration 
        basis by physical insight and practical needs; one should be able to mark those configurations, that are to be removed from a given 
        list and perhaps with some information about the actual size of the simplified computation before this *removal* is made explicit.
    + Since we have a list of levels (with more or less full physical specification with regard to coupling and quantum numbers), it should 
        be not difficult to make use of this list to *set* or *re-define* some of the rates or transition amplitudes. Such a setting 'by hand' 
        will allow to incorporate some selected *decay steps* into the cascade, such a direct DA, two-photon decay, etc, which are not supported
        by JAC itself. 
    + We shall support a *propagation* of the occupation probability to other levels due to the calculated Auger and radiative transition
        amplitudes; if all the probabilities become constant after a finite mumber of steps, the final-level distribution is reached. A display 
        of such a *probability flux* could be controlled by some resolution parameter in order to decide wether nearby-in-energy pathes are 
        to be summed up or are shown only as *summed flux*. For example, a resolution of 1 eV would mean that all channels with energy separation 
        <1 eV are shown by a single arrow together.
    + A Monte-Carlo path generator and distributor of the initial probability might be desirable in order to determine final-level distributions
        and ion yields for complex cascades. This would enable one to determine the probability flux through the cascade by simply counting 
        how often a transition rate is 'used' from some database to form a individual step in some path.
    """
    function decayCascades()  end




    """
    ### `Design principles and limitations of the JAC program`**

    With JAC, we aim to develop and implement tools that are based on a (descriptive) high-level language and which help:

    + perform simple as well as (more) sophisticated atomic comptutations, both at an ab-inito and sometimes also at a semi-empirical basis;
    + provide wave functions, level and transition amplitudes and that supports the simulation of various atomic processes;
    + take advantage of existing atomic codes from the literature or developed world-wide, whenever reasonable;
    + provide a simple access to graphical interfaces and representations of various data and, hence, to enlarge the range of
        possible atomic computations,
    + provide guidance and a (well-designed) framework for implementing future code for modelling (even more) complex processes.


    Despite of its rather wide range of applications, however, **JAC is restricted** to:

    - The computation of atomic properties and models for which the (many-electron) level structures and wave functions 
        can be modelled reasonably well by means of a **radial-spherical representation of all atomic orbitals.**
    - To (atomic) systems and processes for which the atomic states can be properly expanded in a **(finite) basis of 
        time-independent CSF.**
    - This includes the computation of level properties, such a hyperfine-related or induced properties, as well as of atoms in weak and even
        time-dependent fields, as long the level structure of the atoms an ions remain intact. This also enables one to deal with time-dependent
        atomic density matrices and several other atomic properties and processes.
    - **These restrictions exclude** however the treatment of molecules, atoms in (very) strong or short-pulse fields, solid-state effects, etc.


    Three important `design principles of JAC` are:

    ### `Design of a high-level atomic language`

    + **Precise and highly descriptive language:**  A precise language of the underlying physics comes always first. Moreover, a proper 
        decomposition of some physics task into (coarse-grained) computational steps should typically support already a pseudo-code description 
        of the problem and should be reasonably close to its realization within JAC.
    + **ciently general for a wide range of atomic structures, processes and applications:**  Apart from level and transition properties,
        we wish to support the computations of excitation and/or decay cascades, typically within 2-3 different models (average, fine-structure 
        resolved, with shake transitions) in order to cover a good range of applications. 
    + **Non-cryptic and well understandable from a general physics viewpoint:** A language that highlights the underlying physics and avoids 
        (unnecessary) technical slang. The purpose (and use of) all methods should be understandable to a typical atomic physicist without much
        additional training.
    + A **moderate number of (generic) commands:** that highlight the particular action/task, i.e. `actively` describe the task that will be 
        performed and where further details are provided by keystrings and special structs (for instance, to detail some special settings).
    + **Covers the basics functionality of other established codes:** cf. Section Comparison with other (established) atomic structure codes
        above.
    + **Tutorials and examples at different level of complexity:** Cf. the method JAC.tutorials().
    + **Supports and (re-) definition of physical constants, frequent settings and frameworks:** See the help pages ?JAC.Defaults.setDefaults 
        and ?Defaults.convertUnits.
    + **Open-code project:**  We wish to encourage other user to make suggestions, request and improvements to the code.


    ### `Program design and definition of data struct's`

    + **Common source code design:**  The JAC program is designed with a rather 'flat' hierarchy; in the source code, all commands should 
        almost always be called by their full name, e.g. Basics.compute(), JAC.Basics.display(), even if these methods are 'exportet' to the 
        user's level.
    + All **data struct's should be well-adapted to the underlying physics, general assumptions and frameworks:** Indeed, these structs should 
        clearly reflect the 'building blocks' of atomic struture theory, and good care has been undertaken to find a proper balance between 
        types and subtypes. Therefore, a small number of arguments often enables one to perform rather complex tasks.
    + **Close resemblance between pseudo- and working code:** The examples found in JAC.tutorials() demonstrate that there is only a minor 
        difference between pseudo-code and the actual realization of some task with the JAC program.
    """
    function design()  end




    """
    **`Glossary: Further notations from atomic physics`**

    Here, we list and explain (in alphabetic order) a few further terms from atomic structure and collisions theory that might be helpful to
    better understand the structure and notations within the JAC program.


    + `Amplitude`: Transition amplitudes are the many-electron matrix elements upon which atomic structure and collision theory is usually built on.

    + `Basis`: In JAC, a basis refers to a many-electron basis that is specified in terms of a CSF list and the radial orbitals of all (equivalent)
        electrons. Each CSF is specified uniquely by a proper set of quantum numbers; here, the (so-called) seniority scheme is applied for the
        unique classification of the CSF and the evaluation of all matrix elements.

    + `(Non-relativistic electron) configuration`: Describes the occupation of shells within the atomic shell model, for instance, 
        1s^2 2s^2 2p^6 3s^2 .... In the JAC program, closed-shell configurations can also be abbreviated, for instance, by [Ne] 3s^2; cf. 
        Configuration, Atomic.Computation, and at several places elsewhere.

    + `(Relativistic electron) configuration`: Describes the occupation of subshells within the atomic shell model. In relativistic theory, 
        each shell splits into two shell with j = l +- 1/2, apart from the ns_1/2 shells; cf. ConfigurationR.

    + `CSF`: A so-called configuration state function, i.e. a geometrical well-defined state vector of the many-electron Hilbert space, 
        Apart from its geometrical part (due to the its angular structure), it requires the specification of all readial orbitals; cf. Basis.

    + `Grid`: In JAC, all radial orbital functions are always represented on a grid. Only their grid representation is applied in the evaluation
        of all (single- and many-electron) matrix elements. Predefined grids refer to a 'exponential grid', suitable for bound-state computations, 
        as well as a log-lin grid, which increases exponentially in the inner part and linearly in the outer part. The latter grid is suitable for
        collision processes and for dealing with the electron continuum;  cf. Radial.Grid, Radial.Potential, Radial.Primitives.

    + `Level`: Most naturally, the quantum state of free atoms and ions is characterized in terms of (electron) configurations, levels and (magnetic) 
        substates. Since the substates are degenerate, (only) levels are energetically distinguishable for free atoms and are often classified either
        by their energy, total angular momentum and parity |alpha JP> or in LS notation by their term ^2S+1 L_J and total angular momentum. In JAC,
        levels (cf. Level) are the building blocks to describs, level energies, multiplets and atomic processes. Apart from the usual bound-state
        levels of atoms, the same level classification is applied to distinguish between quantum states embedded into the continuum or if they occur
        as 'unperturbed states' in second- and higher-order perturbation theory. 

    + `Light pulse`: a pulse of the em field that is characterized by its shape, (maximum) intensity, duratio, polarization. Typically, such pulses
        are first specified in a form as specified in the literature or by experiment, and are only later adopted to an internal representation 
        that is suitable for a proper time evolution of the density matrix.

    + `Line`: atomic transition that is characterized in terms of a well-defined initial and final level as well as it occurs for the computation
        of various properties, such as cross sections or rates, angular distribution parameters. Typically, a line contains various channels
        (sublines),for instance, due to occurance of multiploles or partial waves in the decomposition of the many-electron matrix elements; 
        cf. Auger.Line, PhotoEmission.Line, PhotoExcitation.Line, PhotoRecombination.Line and elsewhere.

    + `Multiplet`: Atomic levels are naturally 'grouped together' into (so-called) mutiplets; most often, this term just stands for all fine-
        structure levels of a given configurations. More generally, multiplets may refer to any group of levels, for instance, to groups of the same
        same J and/or parity, fine-structure levels of closely related configurations, etc.

    + `Orbital`: In atomic physics, an orbital typically refers to a one-electron function in an radial-spherical representation. Often, only the
        radial function(s) are meant. In the relativistic theory, of course, one needs to distinghuish between the large and small components of
        the orbital, following Dirac's theory. Orbitals [cf. the struct Orbital] are central in JAC to define the atomic bases, multiplets and 
        to built-up proper (many-electron) continua.

    + `Primitives`:  Though the radial parts of all orbital functions are always represented on a radial grid, the self-consistent field as well as
        solutions of the free electron can often be obtained more conviniently within a basis set. Primitives are elements of a (radial) basis and
        may be defined either locally (B-splines, ...) or globally (Slater-type orbitals, Gaussian-type orbitals, ...). In JAC, we make use
        of B-splines ...

    + `Pathway`: In contrast to an (atomic) line, that is characterized by an initial- and a final-level (and the corresponding multiplets), a
        pathway describe a sequence of three or more levels, and which correspond to different atomic processes. These levels are usually referred 
        to as initial, (one or several) intermediate and final level. Pathways occur naturally in dielectronic recombination as well as in various 
        excitation-ionization of excitation-autoioization processes.

    + `Shell/Subshell`: Shells and subshells are the building blocks for the atomic shell model, in which electrons occupy a particular shell or 
        make transitions (quantum jumps) between such shells. In the relativistic theory, each non-relativistic nl-shell (apart from the ns-shells) 
        splits into two relativistic subshells due to j = l +- 1/2. In JAC, the shell and subshell notations are consequently used in order to 
        denote the electron configurations, CSF and the orbitals for equivalent electrons. There are special data structs available to easily deal 
        and communicate the (sub-) shell specification of levels and wave functions.
    """
    function glossary()  end




    """
    **`Data types, structs and name conventions of the JAC module`**

    The use of a proper terminology and data structures has been found essential for developing the JAC module. Below, we list and briefly explain
    these data types and how they appear in atomic theory. Although we presently support just a (small) number of frequently requested *tasks* in 
    atomic structure and collision theory, we tried to define data types that are flexible enough to further extend these tools in the future.
    Following the Julia's standard conventions, all types (struct) are named in CamelCase notation.

    ### `Basic data types`

    + AngularJ64                   ... (positive, half-integer) angular momentum, j = 0, 1/2, 1, 3/2, ... .
    + AngularM64                   ... (half-integer) projection of ang. momentum, m = -1/2, 0, 1/2, ... can be initialized also
                                        w.r.t AngularJ64().
    + ContinuumNormalization       ... method for dealing with the normalization of continuum orbitals.
    + ContinuumPhase               ... method for determining the phase of continuum orbitals.
    + ContinuumSolutions           ... method for solving continuum orbitals.
    + Eigen                        ... represents eigenvalues and eigenvectors if different diagonalization procedures are used.
    + EmMultipole                  ... a multipole of the em field.
    + EmGauge                      ... an allowed gauge form for the em field, for instance, Coulomb, Babushkin, Magnetic, ...
    + EmProperty                   ... a given property in Coulomb (velocity) as well as Babushkin (length) gauge.
    + EmStokes                     ... (computed) Stokes parameter for the polarization of emitted radiation.
    + ExpStokes                    ... (experimentally) given Stokes parameter for the polarization of incoming radiation.
    + Guint                        ... specifier for dealing with graphical user interfaces (GUI).
    + LevelSymmetry                ... total level symmetry (J, parity).
    + Model                        ... to keep the all nuclear parameters.
    + Parity                       ... standard parity values
    + Shell                        ... a non-relativistic shell.
    + SolidAngle`                  ... defines a type for a solid angle Omega = (theta, phi).
    + Subshell                     ... a relativistic subshell.
    + SubshellStateR               ... a relativistic antisymmetric subshell state within the seniority scheme.
    + TensorComp                   ... component of the statistical tensor as associated with an atomic level. 
    + UseGauge                     ... an allowed gauge form requested for explicit computations: UseCoulomb or UseBabushkin.
    + Warnings                     ... for dealing with warnings that are made during a run or REPL session.

    ### `Data types from many-electron theory`

    + AsfSettings                  ... settings for SCF and CI computations.
    + Atomic.Computation           ... atomic computation of a multiplet, including the SCF, CI and transition properties.
    + Basis                        ... (relativistic) atomic basis, including the configuration space and radial orbitals.
    + Bspline                      ... set of B-splines.
    + Configuration                ... (non-relativistic) electron configuration as specified by its shells and their occupation.
    + ConfigurationR               ... (relativistic) electron configuration as specified by its subshells and their occupation.
    + Level                        ... atomic level in terms of its quantum number, energy and a (possible) representation.
    + Multiplet                    ... an ordered list of atomic levels with a name.
    + Orbital                      ... (relativistic) radial orbital function that appears as 'building block' to define many-electron 
                                        states; more often than not, it just occurs as radial orbital on a given (radial) grid while the 
                                        angular dependence is given by the subshell label.
    + Radial.Grid                  ... radial grid to represent the (radial) orbitals.
    + Radial.Potential             ... radial potential function.
    + Radial.Primitives            ... a list of radial functions, that may serve as a set of primitives in SCF computations, together 
                                        with several parameters for its definition.
    + Radial.SingleSymOrbitals     ... a list of radial orbitals with large and small component but of the same symmetry (kappa); such a 
                                        list may serves as (complete) single-electron basis to deal with second- and higher-order processes.
        

    ### `Data types calculating level properties`

    + AtomicLevelProperty          ... an atomic level property that is supported by the JAC module, such as HFS, IsotopeShift, ....
    + Einstein.Settings            ... settings for Einstein A and B coefficients, calculated within a single given Multiplet.
    + Einstein.Outcome             ... (results of the) Einstein A and B coefficients for a single line.
    + Hfs.Settings                 ... settings for HFS A and B coefficients.
    + Hfs.Outcome                  ... (results of the) HFS A and B coefficients for a single level.
    + IsotopeShift.Outcome         ... (results of the) M and F isotope-shift parameters for a single level.
    + IsotopeShift.Settings        ... settings for the M and F isotope-shift parameters.
    + LandeZeeman.sublevelJ        ... specifies a magnetic sublevel with well-defined J.
    + LandeZeeman.sublevelF        ... specifies a magnetic hyperfine sublevel with well-defined F, M_f.
    + LandeZeeman.Outcome          ... (results of the) Lande factors and Zeeman splittings for a single level.
    + LandeZeeman.Settings         ... settings for the Lande factors and Zeeman splitting in an external magnetic field.
    + PlasmaShift.PlasmaModel      ... Model for dealing with the plasma environment.
    + PlasmaShift.Outcome          ... (results of the) plasma shifts and energies for a single level.
    + PlasmaShift.Settings         ... settings for including plasma interactions into the CI matrix.


    ### `Data types for calculating (time-independent) atomic processes`

    + AtomicProcess                     ... an atomic process that is supported by the JAC module, such as Auger, photo, ....
    + AlphaVariation.Outcome            ... outcome of a alpha-variation computation, such as the K enhancement.
    + AlphaVariation.Settings           ... seetings for computing alpha variation parameters.
    + Auger.Channel                     ... Auger channel of well-defined energy and partial outgoing wave.
    + Auger.Line                        ... Auger line between (two) specified initial- and final-state levels and with (possible) subchannels.
    + Auger.Settings                    ... settings for computing Auger lines.
    + CoulombExcitation.Channel         ... Coulomb excitation channel of well-defined energy and partial wave.
    + CoulombExcitation.Line            ... Coulomb excitation line with (possible) subchannels.
    + CoulombExcitation.Settings        ... settings for computing Coulomb excitation  lines.
    + CoulombIonization.Channel         ... Coulomb ionization channel of well-defined energy and partial wave.
    + CoulombIonization.Line            ... Coulomb ionization line with (possible) subchannels.
    + CoulombIonization.Settings        ... settings for computing Coulomb ionization  lines.
    + Dielectronic.Channel              ... dielectronic-recombination channel of well-defined multipolarity and gauge as well as energy and 
                                            partial incoming wave.
    + Dielectronic.Line                 ... dielectronic recombination line between (three) specified initial-, intermediate and final-state 
                                            levels and with (possible) subchannels.
    + Dielectronic.Resonance            ... single dielectronic resonance that summarizes all Dielectronic.Line's for some fixed intermediate
                                            level within the continuum. 
    + Dielectronic.Settings             ... settings for computing dielectronic recombination lines.
    + DecayYield.Outcome                ... outcome of a decay yield computation.
    + DecayYield.Settings               ... settings for computing decay yields lines.
    + DoubleAuger.Channel               ... DoubleAuger channel of two partial outgoing waves with well-defined energy.
    + DoubleAuger.Line                  ... DoubleAuger line between (two) specified initial- and final-state levels and with (possible) 
                                            subchannels.
    + DoubleAuger.Settings              ... settings for computing DoubleAuger lines.
    + ImpactExcitation.Channel          ... electron-impact excitation channel of well-defined energies. partial waves and phases of the 
                                            incoming and outgoing electrons.
    + ImpactExcitation.Line             ... electron-impact excitation line between (two) specified initial- and final-state levels and with 
                                            (possible) subchannels.
    + ImpactExcitation.Settings         ... settings for computing electron-impact excitation lines.
    + ImpactExcitationAutoion.Channel   ... electron-impact excitation channel of well-defined energies, partial waves and phases of the 
                                            incoming and outgoing electrons.
    + ImpactExcitationAutoion.Pathway   ... electron-impact excitation line between (two) specified initial- and final-state levels and with 
                                            (possible) subchannels.
    + ImpactExcitationAutoion.Settings  ... settings for computing electron-impact excitation lines.
    + ImpactIonization.Channel          ... electron-impact ionization channel of well-defined energies, partial waves and phases of the 
                                            incoming and outgoing electrons.
    + ImpactIonization.Line             ... electron-impact ionization line between (two) specified initial- and final-state levels and with 
                                            (possible) subchannels.
    + ImpactIonization.Settings         ... settings for computing electron-impact ionization lines.
    + MultiPhotonDeExcitation.Channel   ... multi-photon excitation or decay channel with well-defined multipolarities and gauge.
    + MultiPhotonDeExcitation.Line      ... multi-photon excitation or decay line between (two) specified initial- and final-state levels 
                                            and with (possible) subchannels.
    + MultiPhotonDeExcitation.Settings  ... settings for computing multi-photon excitation or decay lines.
    + MultiPhotonIonization.Channel     ... multi-photon ionization channel with well-defined multipolarities, gauge as well as energy and 
                                            partial wave of the outgoing electron.
    + MultiPhotonIonization.Line        ... multi-photon ionization line between (two) specified initial- and final-state levels and with 
                                            (possible) subchannels.
    + MultiPhotonIonization.Settings    ... settings for computing multi-photon ionization lines.
    + MultiPhotonDoubleIon.Channel      ... multi-photon double ionization channel with well-defined multipolarities, gauge as well as energy
                                            and partial waves of the (two) outgoing electrons.
    + MultiPhotonDoubleIon.Line         ... multi-photon double ionization line between (two) specified initial- and final-state levels and 
                                            with (possible) subchannels.
    + MultiPhotonDoubleIon.Settings     ... settings for computing multi-photon double ionization lines.
    + PairAnnihilation1Photon.Channel   ... positron-bound-electron pair annihilation (PEPA) with single-photon emission channel of 
                                            well-defined multipolarity, gauge as well as energy and partial incoming (positron) wave.
    + PairAnnihilation1Photon.Line      ... PEPA with single-photon emission line between (two) specified initial- and final-state levels and 
                                            with (possible) subchannels.
    + PairAnnihilation1Photon.Settings  ... settings for computing PEPA with single-photon emission lines.
    + PairAnnihilation2Photon.Channel   ... positron-bound-electron pair annihilation (PEPA) with two-photon emission channel of well-defined 
                                            multipolarities, gauge as well as energy and partial incoming (positron) wave.
    + PairAnnihilation2Photon.Line      ... PEPA with two-photon emission line between (two) specified initial- and final-state levels and 
                                            with (possible) subchannels.
    + PairAnnihilation2Photon.Settings  ... settings for computing PEPA with two-photon emission lines.
    + PairProduction.Channel            ... positron-bound-electron pair production (PEPP) by single-photon absorption channel of well-defined 
                                            multipolarity, gauge as well as energy and partial outgoing (positron) wave.
    + PairProduction.Line               ... PEPP by single-photon absorption line between (two) specified initial- and final-state levels and 
                                            with (possible) subchannels.
    + PairProduction.Settings           ... settings for computing PEPP lines.
    + PhotoExcitation.Line              ... photoexcitation line between (two) specified initial- and final- state levels and with (possible
                                            JAC.PhotoEmission.Channel) subchannels.
    + PhotoExcitation.Settings          ... settings for computing photoexcitation lines.
    + PhotoExcitationAutoion.Channel    ... photo-excitation autoionization channel of well-defined energies of the incoming photon as well as 
                                            the partial wave and phase of the outgoing electron.
    + PhotoExcitationAutoion.Pathway    ... photo-excitation autoionization pathways between (three) specified initial-, intermediate and 
                                            final-state levels and with (possible) subchannels.
    + PhotoExcitationAutoion.Settings   ... settings for computing photo-excitation autoionization pathways.
    + PhotoIonization.Channel           ... photoionization channel of well-defined multipolarity, gauge as well as energy and partial 
                                            outgoing wave.
    + PhotoIonization.Line              ... photoionization line between (two) specified initial- and final-state levels and with (possible) 
                                            subchannels.
    + PhotoIonization.Settings          ... settings for computing photoionization lines.
    + PhotoRecombination.Channel        ... Rec channel of well-defined multipolarity and gauge as well as energy and partial incoming wave.
    + PhotoRecombination.Line           ... radiative electron capture line between (two) specified initial- and final-state levels and with 
                                            (possible) subchannels.
    + PhotoRecombination.Settings       ... settings for computing radiative electron capture lines.
    + PhotoEmission.Channel                 ... radiative channel of well-defined multipolarity and gauge.
    + PhotoEmission.Line                    ... radiative line between (two) specified initial- and final-state levels and with (possible) sublines.
    + PhotoEmission.Settings                ... settings for computing radiative lines.
    + RadiativeAuger.Channel            ... RadiativeAuger channel of a partial outgoing waves and one photon with well-defined energy.
    + RadiativeAuger.Line               ... RadiativeAuger line between (two) specified initial- and final-state levels and with (possible) 
                                            subchannels.
    + RadiativeAuger.Settings           ... settings for computing RadiativeAuger lines.
    + Radiative.Settings                ... settings for computing radiative lines.
    + RayleighCompton.Channel           ... RayleighCompton channel of an incoming and outgoing photon with well-defined energy.
    + RayleighCompton.Line              ... RayleighCompton line between (two) specified initial- and final-state levels and with (possible) 
                                            subchannels.
    + RayleighCompton.Settings          ... settings for computing RayleighCompton lines.
    + REDA.Channel                      ... resonant electron-excitation (sequential) double-autoionization (REDA) channel of well-defined 
                                            energies, partial waves and phases of the incoming and outgoing electrons.
    + REDA.Pathway                      ... resonant electron-excitation (sequential) double-autoionization (REDA) pathways.between (four) 
                                            specified initial-, (two) intrmediate and final-state levels and with (possible) subchannels.
    + REDA.Settings                     ... settings for computing resonant electron-excitation (sequential) double-autoionization (REDA) 
                                            pathways.

    ### `Data types for calculating (time-dependent) atomic processes`

    + Pulse.Envelope                       ... defines a type for the envelope (function) of an em pulse with well-defined time delay, 
                                                amplitude and (normalized) shape function.
    + Pulse.ExperimentalCharacterization   ... to characterized an experimental or physically described em pulse in terms of its 
                                                propagation direction, frequency, maximum intensity, pulse length or No. of cycles, time-delay, 
                                                polarization, etc., i.e. of what is easily accesssible by an experiment.
    + Pulse.Gaussian                       ... a Gaussian light pulse that is used for evaluating time-dependent statistical tensors.
    + Pulse.Polarization                   ... defines the polarization of an em pulse in terms of its linear and circular degrees, the 
                                                direction of the polarization vector or some generalized polarization coefficients.
    + Pulse.PolarizationType               ... defines the polarization of an experimentally described light pulse as linear, left-circular, ..
    + Pulse.Shape                          ... defines a shape of a general em pulse as Gaussian, SineSquared, etc.
    + Pulse.SineSquared                    ... a SinSquared light pulse that is used for evaluating time-dependent statistical tensors.

    ### `Data types for dealing wiht (time-dependent) statistical tensors`

    + Statistical.ResonanceR      ... a resonance state in the continuum with a well-defined bound-ionic core, one or several electrons 
                                        in the continuum, a widths as well as a loss rate due to *additional* decay processes that cannot be 
                                        accounted for explicitly.
    + Statistical.Tensor          ... represents a statistical tensor of given rank k, projection q and which generally depends upon two 
                                        resonances.

    ### `Data types for advanced computations`

    + Atomic.CasComputation          ... an individual or a series of systematically enlarged SCF computations.
    + Atomic.CasStep                 ... single-step in an (systematically enlarged) SCF calculation.
    + Atomic.CasSettings             ... settings for CAS computation.

    + Cascade.Approach               ... a particular (computational) approach in which a cascade is considered.
    + Cascade.Block                  ... a block of configurations that are treatet together within a given cascade. 
    + Cascade.Data                   ... all transition data of a cascade as given by a list of lines (of different type).
    + Cascade.Computation            ... definition of an atomic exciation/decay cascade from which the actual computations can be derived.
    + Cascade.Level (mutable)        ... defines a level specification for dealing with cascade transitions.
    + Cascade.LineIndex              ... defines a line index with regard to the various lineLists of data::Cascade.Data.
    + Cascade.Step                   ... an individual step of a Cascade.Computation that generally combines two ionization states of ions.
    + Cascade.Simulation             ... simulation of cascade data.
    + Cascade.SimulationSettings     ... defines settings for performing the simulation of some cascade (data).
    + Cascade.Settings               ... settings for cascade computations (not yet).
    """
    function datatypes()  end




    """
    **`Interactive use of JAC procedures`**

    Various functions (methods) are available in JAC.Basics which support the interactive generation and manipulation of data; 
    for example, use ? Basics.compute for all further details:

    +      compute()            ... to compute angular coefficients; CI matrix; radial orbitals; radial potentials; ...
    +  Defaults.convertUnits()  ... to convert units for energy, cross sections, rates, times, etc.; wave numbers; ...
    +  Basics.display()         ... to display constants and settings; ...
    +      estimate()           ... to estimate ionization potentials; ...
    +      generate()           ... to generate condensed multiplets; relativistic and non-relativistic configuration lists;
                                    order shell and subshell lists; single-electron spectra; ...
    +      getDefaults()        ... to give various physical constants; order shell and subshell lists; ...
    +      integrate()          ... to integrate a function on a grid; ...
    +      merge()              ... to merge (two or more) atomic basis; multiplete; ...
    +      setDefaults()        ... to define the framework; units for energy, cross sections, rates, etc.; standard grid; 
                                    subshell list; methods for solving continuum orbitals and their normalization; ...
    +      sort()               ... to sort the levels of (hyperfine) mutliplets by energy; ...
    +      tabulate()           ... to tabulate (hyperfine) mutliplets by energy; ...

    """
    function interactive()  end



    """
    **`Atomic processes`**

    In addition, the computation of the following excitation, ionization and decay processes is supported for atoms and ions with N electrons:

    + Auger              ......... **Auger transitions**, i.e. the autoionization or emission of a (free) electron into 
                                    the continuum; |i(N)> --> |f(N-1)> + e_A. 
    + Compton            ......... **Rayleigh and Compton** elastic and inelastic photon scattering cross sections; 
                                    |i(N)> + omega --> |f(N)> + omega'.
    + Conversion         ......... **Internal conversion**; |i(N)> + [nucleus^*]] --> |f(N-1)> + e_c.
    + CoulExc            ......... **Coulomb excitation of target or projectile ions** by fast, heavy ion projectiles; 
                                    |i(N)> + Ion --> |f(N)> + Ion.
    + Dierec             ......... **Di-electronic recombination**, i.e. the dielectronic capture of a free electron and 
                                    the subsequent emission of photons; |i(N)> + e --> |m(N+1)> --> |f(N+1)> + omega.
    + Eimex              ......... **electron-impact excitation** cross sections and collision strengths; 
                                    |i(N)> + e --> |f(N)> + e' (not yet).
    + EimexAuto          ......... **electron-impact excitation and subsequent autoionization**; cross sections; 
                                    |i(N)> + e --> |m(N)> + e' --> |f(N-1)> + e' + e_A (not yet).
    + EimexIon           ......... **electron-impact ionization** cross sections; |i(N)> + e --> |f(N-1)> + e' + e_i (not yet).
    + MultiPhoton        ......... **multi-photon excitation and decay**; amplitudes and cross sections; 
                                    |i(N)> + n*omega --> |f(N)>  or  |i(N)> --> |f(N)> + n*omega.
    + MultiIon           ......... **multi-photon ionization**; |i(N)> + n*omega --> |f(N-1)> + e_p (not yet).
    + MultiDoubleIon     ......... **multi-photon double ionization**; |i(N)> + n*omega --> |f(N-1)> + e_p + e_p' (not yet).
    + PhotoExc           ......... **photon-impact excitation** cross sections; |i(N)> + omega --> |f(N)>.
    + PhotoExcFluor      ......... **photon-impact excitation and fluorescene** cross sections; 
                                    |i(N)> + omega --> |m(N)> --> |f(N)>+ omega.
    + PhotoExcAuto       ......... **photon-impact excitation and autoionization** cross sections; 
                                    |i(N)> + omega --> |n(N)> --> |f(N-1)> + e_A.
    + PhotoIon           ......... **Photoionization processes**, i.e. the emission of a single free electron into the 
                                    continuum due to an external light field.; |i(N)> + omega --> |f(N-1)> + e_p.
    + PhotoIonFluor      ......... **Photoionization processes with subsequent fluorescence**; 
                                    |i(N)> + omega --> |m(N-1)> + e_p --> |f(N-1)>  + e_p + omega'.
    + PhotoIonAuto       ......... **Photoionization processes with subsequent autoionization**; 
                                    |i(N)> + omega --> |m(N-1)> + e_p --> |f(N-1)>  + e_p + e_a  (not yet).
    + Radiative         ......... **Radiative (multipole) transitions** between bound-state levels of the same charge state;
                                    transition probabilities; |i(N)> --> |f(N)>  + omega.
    + RadAuger           ......... **Radiative Auger**; simultaneous photon emission and autoionization cross sections; 
                                    |i(N)> --> |f(N-1)> + e_A + omega.
    + Rec                ......... **Radiative electron capture**, i.e. the capture of a free electron under the simultaneous 
                                    emission of a photon; |i(N)> + e --> |f(N+1)> + omega.

    """
    function processes()  end



    """
    **`Atomic properties`**

    Apart from approximate level energies and eigenvectors, JAC (will) support the computation of the following level properties:

    + AlphaX        ... alpha variations; differential sensitivity parameters.
    + Einstein      ... Einstein A, B coefficients and oscillator strength; although these coefficients are not an 
                        original level property, the Einstein module treats these computations within a single 
                        basis/multiplet and, hence, cannot include relaxation effects, etc. The Einstein 
                        feature of JAC can be used, however, for a quick overview to transition probabilities 
                        or in order to simplify cascade computations.
    + FormF         ... Standard and modified atomic form factors.
    + Greens        ... Greens function of an atomic level.
    + HFS           ... Hyperfine A and B parameters.
    + Isotope       ... Isotope shift M and F parameters.
    + LandeJ        ... Lande g_J factors.
    + LandeF        ... Lande g_F factors.
    + Polarity      ... Static and dynamic polarizibilities of atomic levels.
    + Plasma        ... CI computations including interactions from various plasma models.
    + Yields        ... Fluoerescence and Auger yields of atomic levels.

    """
    function properties()  end



    """

    ### `Use of (em) light pulses in the time evolution of statistical tensors`**

    Much recent interest focused on the interaction of atoms and ions with intense light fields. If the incident fields are not to strong,
    i.e. for Keldysh parameters lambda << 1, the time-evolution of an atomic cloud can be described by a Liouville-type equation for the atomic 
    density matrix rho. This density matrix, or the closely related statistical tensors, then evolve within a pre-defined atomic (many-electron) 
    basis  {Phi_n == (alpha J M_J), n = 1,...,n_max},  where J and M_J are as usual the total angular momentum and its projection, and where
    alpha all further quantum numbers that are needed to specify the atomic state uniquely. The time-evolution of the density matrix proceeds 
    due to the interaction of the atoms with the electric and magnetic multipole components (operators) of the radiation field. To properly 
    describe such an evolution, of course, care has to be taken about the characterization of the light pulses, their shape and intensity as well as 
    polarization:


    ### `(Experimental) Characterization of light pulses`

    + Like in experiments, a pulse is characterized by its frequency, direction of propagation, (normalized) shape function, maximum intensity,
        time delay and the polarization of the pulse.
    + To specify the propagation direction and polarization, all angles are defined with respect to a fixed laboratory frame Sigma = (x,y,z), 
        and where all coordinate systems are assumed to be `right-handed'.
    + The wave vector of the `first pulse' typically defines the z-axis and its linear polarization (with angle phi_0) the x-axis of the 
        fixed frame, i.e. the x-z scattering plane. 
    + For a circularly polarized pulse, in contrast, the x-axis can either be specified arbitrarely of indirectly due to the the definition
        of other pulses.
    + For all other pulses (n > 1), the program requests to specify the  angles  phi_k and theta_k  to determine the z'-axis as well as 
        bar{phi}_k, bar{theta}_k  to determine the x'-axis; then, the angle  chi_k  specifies again the linear polarization in the primed system 
        Sigma'. It is checked internally that the z' and x'-axes are specified orthogonal to each other.
    + The maximum intensity of the pulses is given in W/cm^2; here, we need a convention how this maximum intensity is 'converted' into a proper
        amplitude of the vector potential. In practice, the relationship between intensity and electric field slightly depends also on 
        polarization, and if a linearly polarized radiation beam is assumed. In addition, this conversion factor also depends of the pulse 
        shape. Since, however, the experimentalists never know the intensity very accurately, a simplified normalization procedure is typically
        sufficient to deal with most situations.
    + The atomic unit of electric field is  m_e^2 e^5/ hbar^4 = 5.1422 x 10^9 V/cm.  From this, we find 
        E_0  [a.u.] = 5.34 x 10^{-9} sqrt{I [W/cm^2]}. This conversion factor has been often used in various TDSE codes in order to obtain
        the electric field envelope amplitude, when experimentalists say that `the (peak) laser intensity in the pulse is I [W/cm^2]`.

    ### `Pulse shapes`

    + The present implementation of the code deals with several (pre-) defined symmetric shapes of em pulses, such as a **Gaussian** pulse, a
        **SineSquared** pulse or a **SineSquaredPlateau**.
    + These symmetric pulses are declared in terms of their *central* frequency omega, the delay time T_d, i.e. the delay of the central 
        time of the pulse shape with regard to the reference time t_0 = 0. 
    + The pulse is described by a (normalized) shape function f_s (t) = f(t)/f_o = 1  and, hence, the intensity of the em field must be 
        captured into the constant f_o.

    ### `Pulse polarization`

    + Four frequently occuring types of polarization are predefined: linear, left-circular, right-circular and elliptical.
    + In a first implementation, we shall restrict ourselves to pure polarization states of the incoming photons, i.e. to cases in
        which P_1^2 + P_2^2 +P_3^2 = 1; this is equivalent to  a pure polarization density matrix  rho^2 = rho. In this case, an arbitrary
        polarization state can be described as linear combination of left- and right-circular components.
    + While the pulses are defined from the viewpoint of experiment, they need to be *expressed* in tensor form to be appropriate for a pulse
        propagation. To this end, the (complex-valued) polarization parameters g_{+1} and g_{-1} need to be determined.

    ### `Time-integration of the statistical tensors`

    + Various methods can be applied to solve the coupled first-order odg's for the time-evolution of the statistical tensors. Apart from 
        methods, which are able to solve the equations with predefined accuracy but typically require an evaluation of the rhs at any times t,
        we are mainly interested here in so-called `shooting' or predictor-corrector methods. These latter methods are purely based on the 
        knowledge of rho and d rho / dt == rhs on some equidistant time grid.
    + Two particularly simple methods with a equidistant grid are Euler's method and the Adams-Bashford methods.


    ### `Observables from the time-dependent statistical tensors`

    + Level population during and after the interaction with one or several given em pulses.
    + Polarization of a selected atomic level alpha J during and after the interaction with one or several given em pulses.
    + Spontaneous photon emission from the atom perpendicular to the pulse propagation
    + Angular distribution of emitted electrons, either in total or associated with some particular final-ionic state
    """
    function pulses()  end




    """
    **`Why Julia ?`**

    Here, we recall a few remarks from the literature as well as some own experience why Julia is helpful for developing the JAC program. 
    Many of these arguments have been adapted from the paper by Bezanson et al. (2017) and by Post and Kendall (2004):

    +  Julia stands for the combination of productivity and performance through a careful language design and carefully chosen technologies; 
        it never forces the user to resort to C or Fortran for fast computations. -- Julia's design allows for gradual learning of modern concepts 
        in scientific computing; from a manner familiar to many users and towards well-structured and high-performance code. In Julia, however,
        it is sometimes beliefed that class-based methods are not scientifically powerful enough to express full abstraction in scientific computing.

    + `High-level languages`: Most traditional high-level languages are hampered by the overhead from the interpretor and which typically results into
                                more run-time processing that are strictly necessary. One of these hindrances is (missing) type information, and which
                                then results in the request for supporting vectorization. Julia is a 'verb'-based language in contrast to most 
                                object-oriented 'noun'-based language, in which the generic functions play a more important role than the datatypes.

    + `Code selection`: Julia name space allows the use of the same vocabulary in different circumstances, and which makes programs easier to read.
                        In particular, it uses the same mechanism of code selection at the lowest and highest levels and is therefore able to select
                        the right method, either already at compile time or later at run time.

    + `Multiple dispatch`: refers to the dynamically selected implementation and to the concept of running the right code at the right time.
                            This is achieved by overloading by multiple-argument function, a very powerful abstraction. Multiple dispatch makes 
                            it easier to structure the programs close to the underlying science. It also reduces the needs for argument checking
                            at the begin of a function. The overloading of functions by multiple dispatch is also called ad-hoc polymophism.
                            Instead of encapsulating methods inside classes, Julia's multiple dispatch is a paradigm in wich methods are defined on
                            combinations of data types (classes). Julia shows that this is remarkably well-suited for numerical computing.

    + `Code re-use`: In good language design, on should be able to write a general implementation and, as long as the necessary operations are
                        available, the code should just work.

    + `Type system`: Julia's expressive type system that allows opional type annotations; this type system supports an agressive code specializiation
                        against run-time types. Over a large extent, however, Julia code can be used without any mentioning of types (in contrast to
                        C and Fortran); this is achieved by data-flow interference.-- User's own types are also first class in Julia, that is there 
                        is not meaningful distinction between built-in and user-defined types. There are mutable and (default: immutable) 
                        composite types.

    + `Parallelization`: One of the central motivation to built Julia was the design of a parallel computing language. Therefore, Julia provides
                            different facilities for parallelism. Two important concepts are 'remote calls'  and 'remote references'; cf.
                            @parallel. In contrast, vectorization is in Julia NOT considered as a pre-requisite for performance.

    + `Performance`: There are helpful macros, such as @timing function_call(parameters) or @benchmark function_call(parameters) to analyze the
                        performance of the program and to find (and resolve) bottlenecks.

    + `Code developers`: In numerical and scientific computing, people with special skills are sometimes called library or package writers, while 
                        many others just use these libraries and packages.

    + `LAPACK`: All of LAPACK is available in Julia, not just the most common functions. LAPACK wrappers are fully implemented by 'ccall' and can
                be called directly from the Julia promt.

    + `Macros`: A macro is a function that runs at parse time. It takes symbolic expressions in and returns transformed expressions out, which are
                inserted into the code for later compilation. The output of macros is often 'inlined' into the code.

    + `Physical models`: However, better physics is more important than better computer science. It is recommended to use modern but well-proven
                            computer-science techniques, and a 'physics code' should not be a computer-science research project. Instead, one should
                            use best engineering practices to improve quality rather tha processes. Emphasis should be given to improvements of the 
                            physics capabilities. Do not use the latest computer-science features; let the new ideas mature first. Better physics is the
                            most important product of the code.
    
    + The scale of code-development can be truly immense; a good overview/quantitative database about (previously) successful software projects is
        required for good estimation for resources and schedules. It is easy to loose motivation on a project that last years and which has few
        incremental deliveries.

    + `Success criteria`: One of the important success criteria is the costumer focus. What do the user really need.

    + `Visualization`: Visualization of large data sets is essential for debugging, problem generation and analysis of results.

    + `Code evolution`: Continues replacement of code modules is recommended as better tools and techniques are developed. Every code development
                        typically proceeds in steps: First develop a core capability (with a small team) and let this small core be tested by users 
                        and, if successful, add further capabilities (so-called incremental delivery).

    + `Team work`: Good teams are more important than good processes; they are characterized by the ability to share ideas and to work informally
                    together. The chracks need to evolve into senior mentors and code architects who communicate their vision and expertise. --- Poor 
                    performers not only do not complete their work, they also tend to impede and discourage others.

    + `Code specification`: Some flexibility in the requirement specification phase is essential because it is difficult to predict when (or if) a 
                            new algorithm/approach will be available. There is a need to pursue multiple approaches for algorithms and modules near 
                            to the critical path. If one approach is not feasible, another one can be used.

    + `Quality tests`: If a code is not verified and validated by proper examples, the users will not believe and see its connection to reality. 
                        Therefore, develop and execute a verification and validation program.
    """
    function whyJulia()  end



end # module

    
