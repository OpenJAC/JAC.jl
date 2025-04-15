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

# Installation

The development version of JAC can installed installed directly in the interactive Julia REPL (Read-Eval-Print-Loop)
Press `]` to enter the Julia package mode and type `add https://github.com/OpenJAC/JAC.jl` to install JAC.

```
julia> ]
(@v1.10) pkg> add https://github.com/OpenJAC/JAC.jl
```

# Development

As a **Scientific package** users are encouraged to explore the capabilities of JAC, and modify the package as per their requirement. JAC comes with an MIT 'Expat' [License](@ref jac-license).

For development purpose, user can install JAC from GitHub to their **desired direcory** using `git clone` in the terminal. This will download a copy of JAC from the development branch of GitHub repository.
```
git clone https://github.com/OpenJAC/JAC.jl.git
```

Next `cd` to the `JAC.jl` directory and start a new Julia session

```
shell> cd JAC.jl
shell> julia
```

Then activate the `JAC.jl` environment and install the dependancy packages of `JAC.jl`
```
julia> ]                   # Change to the Julia package mode
(@v1.10) pkg> activate .
(JAC) pkg> instantiate     # Installs the dependancy packages of JAC
(JAC) pkg> resolve         # Required **Only if** there is a cinflict in the dependancies
(JAC) pkg> develop .       # Adds the current directory path to Julia local Registry
```