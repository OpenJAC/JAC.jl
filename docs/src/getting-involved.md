
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
writing some new *demos* (Pluto notebooks), in addition to those provided above, in order to navigate others to the task of
a new user. Also, new graphical user interface and plotting features on different outcomes of atomic computations will be very 
helpful for the community. A few further suggestions for extending and improving JAC can be found in section 1.7 in the 
[User Guide, Compendium & Theoretical Background to JAC](https://github.com/AlokaSahoo/JAC.jl/blob/master/docs/UserGuide-Jac.pdf).

