# Getting started with Julia (in REPL)

!!! info 
    Link to the Pluto jl and direct start with  example...

Here, we shall **not introduce** Julia's syntax and concepts for which many tutorials are available on the web. Instead, we just wish to remind and highlight some simple (syntax) features that help to go easier around with JAC, and especially for occasional users from experiment or teaching. This reminder aims to lower the initial *threshold* for users that have been trained on other languages in the past. Here, we shall *pick up* some issues whose physics background is explained only later in other tutorial. Obviously, however, Julia is a very rich and powerful language with many features that go well beyond of what is (and will ever) needed for JAC.

In brief, JAC provides tools for performing atomic (structure) calculations of different kind and complexity, and for which further details are given in the tutorials below. To see anything from JAC, we shall first invoke the tools by:

```@repl
using JAC
```

a line that will appear at the beginning of all subsequent tutorials. -- A first powerful and frequently needed feature refers to Julia's help pages or just the "?". By typing, for instance, atom or computation

```
? atom
```

```
search: atomic Atomic AtomicState AtomicCompass AtomicStructure @atomic @atomicswap @atomicreplace

Couldn't find atom
Perhaps you meant atomic, atan, acot, acos, htol, hton, ltoh, ntoh, Atomic, @atomic, pathof, atand, atanh, cat, match, catch, stat, acotd, acoth, add, ans, abs, abs2, acosd, acosh, acsc, all, all!, any, any!, asec, asin, axes, as, tan or Beam
No documentation found.

Binding atom does not exist.
```

we see, that `atom` itself is not a well-defined term in the JAC toolbox but that there exists a number of related terms, such as `Atomic`, `AtomicState` (two modules of JAC) and others. We shall not enter here the modular structure of the JAC toolbox but start much simpler with: 

```
? Orbital
```

```@docs    ; canonical=false
Orbital
```

which, apart from its formal meaning, is a particular data structure (`struct`) of JAC and which represents a relativistic orbital (function) including additional information that appears helpful in the given implementation. There are very many (say, more than 300) of such data struct's specified in the JAC toolbox, and thus quite obvious that nobody will remember the details of all these definitions. Indeed, the "?" is the right and a powerful means to remind yourself and make use of these data structures whenever necessary. Special care has been taken that all data structures and functions/methods comes with a reasonable explanation (docstring) in order to work efficiently with JAC.

For instance, we might ask of what can be *added* to each other in JAC:
```
? add
```

```@docs    ; canonical=false
add
```

Apart from a short explanation, these docstring always tell the user (i) in which module the method is defined; (ii) which arguments it takes, including Julia's *multiple dispatch* feature as well as (iii) the type of the return value. All these information are typically relevant to the user, especially if some input or output does not behave as it should. Indeed, the complexity can grow quite rapidly, for instance, if we ask for help of what we can `generate`:

```
? generate
```

```@docs    ; canonical=false
generate
```

Well, this is quite a lot, and we shall explain some of these methods below; a similar or even larger output, you can generate by `? perform` as well as few other terms that are central to the implementation of JAC.

**Constructors & program control:$\quad$** Another frequent use of the (help) "?" concerns the data flow and control of almost all computations. In JAC, we often make use of (so-called) `Settings` that enable the user to overwrite default values or to *control* the computation to the extent, he or she wishes to have control. These `Settings` are context dependent and are different for each atomic property or process that can be computed by the JAC toolbox. They are defined in the various modules and need to be specified accordingly. For instance, to control the computation of transition probabilities for the (fine-structure) levels between given initial- and final-state configuration, one has to overwrite the (defaults) settings:

```
? PhotoEmission.Settings
```

```@docs    ; canonical=false
PhotoEmission.Settings
```

We shall meet these and (many) other settings quite often in the tutorials below. --- Beside of Julia's help features (?), however, it is sometimes difficult to remember the right term or function name. In this case, it easy to make a <double-tab> after the dot (notation) or to make use of the (Unix/Linux) `grep` command within the `JAC/src` directly. Similar line-search commands will exist also at other platforms. In particular, for those of you who wishes to support and extend the JAC toolbox, the dot expansion and the `grep` command will be found very helpful, perhaps more than other search tools.


!!! @info "Info"
    Users can also refer to the API Reference section that provides a comprehensive list of JAC declared `Type`s and `Function`s for selected atomic processes.

**Use of constructors:$\quad$** Another Julia feature, that is frequently applied in JAC, is the successive definition of constructors in order to set-up complex data structures. This features is applied, for instance, in order to define an `Atomic.Computation` or a `Cascade.Computation` as a whole. We shall explain these rather complex data types below in different tutorials. The same issue appears however already at a much simpler level. For example, if we wish to select (specify) a number of levels from a multiplet prior to some particular -- configuration interaction -- computation, we can make use of a
```
? LevelSelection
```

```@docs ; canonical=false
LevelSelection
```
Apart from the logical flag `active`, such a level selection requires to either specify a list of level numbers (indices) or *level symmetries*

```@repl
using JAC   # hide
LevelSelection(true, indices= [i for i in 1:11])
LevelSelection(true, symmetries= [LevelSymmetry(1//2, Basics.plus)])
```

Here, we made use of a LevelSymmetry to specify the overall rotational $J^P$ symmetry of atomic levels.

```
? LevelSymmetry
```

```@docs ; canonical=false
LevelSymmetry
```

As seen from this definition, the level symmetry just comprises the total angular momentum (of type `AngularJ64`) and the parity of the level (of type `Parity`). Therefore, the specification of a list of level symmetries in `LevelSelection` already requires to nest four constructors in order make the level selection explicit: (i) For the angular momentum, (ii) the parity, (iii) the level symmetry and (iv) to create a list (array) of such level symmetries. All the constructors can be specified and built together also in subsequent steps, such as:

```@julia
using JAC   # hide
J1    = AngularJ64(1//2);           J2 = AngularJ64(5//2)  
pl    = Basics.plus;                mn = Basics.minus
lsym1 = LevelSymmetry(J1, pl);      lsym2 = LevelSymmetry(J2, mn)
levelsyms = [lsym1, lsym2]
```

or simply by *nesting* all the information within a single step
```@repl
using JAC   # hide
levelsyms = [LevelSymmetry(AngularJ64(1//2), Basics.plus), LevelSymmetry(AngularJ64(5//2), Basics.minus)]
```

Both way have their pros and cons, and often some *mixture* is applied where complex constructors are first assigned to some variables, and which are later utilized to built up constructors of higher complexity. --- To finally specify aǹ instance of a `LevelSelection`, we use (onc more) its second constructor above:

```@repl
LevelSelection(true, indices=[1,2,3], symmetries=levelsyms)
```

and which will tell the JAC program to compute the lowest three levels (1, 2, 3) as well as all levels with 1/2+ and 5/2- symmetry. Apart from the selection of individual levels, it is often helpful for the computation of atomic processes to make a prior `LineSelection` and in some cases even a `PathwaySelection` as, for instance, for dielectronic recombination processes. Some of these features will be explained below in subsequent tutorials of JAC:

**Functions & methods:$\quad$** Like most other languages, Julia is based on the successive work through functions and methods; a **function** is first of all specified by its name and it maps a tuple of argument values upon a return value. For instance, the function

```@repl
function addSomething(a, b)
    c = a + b
end
addSomething(3.1, 5//2)
```

can be used to *add* two numbers, rather independent of their particular type, and which are *infered* here automatically. However, additional type declarations might help to *specialize* a function and to ensure **type stability**:

```@repl
function addSomething(a::Int64, b::Int64, c::Int64)
    d = a + b + c
end
```

While the function name is the same in both of these examples above, Julia carefully distinguishes between these two **methods** of the function `addSomething` that may differ by the type *and/or* the number of arguments. This multiple *use* **(dispatch)** of function name enables the user to write highly specialized code. Although a proper (and specialized) definition of functions is often very important for the performance of the program, we shall not discuss such technical issues here. Let us just mention, that a function/method may also return `nothing`:

```@repl
typeof(nothing)
```

In JAC, the value `nothing` is usually returned by all *display* functions that print some selected data or tabulation to screen or elsewhere but does not return a value otherwise.

**Code failures:$\quad$** Beside of its large flexibility and user-friendliness, JAC might terminate from time to time for *non-obvious* reasons. Since JAC is first of all a *physics code*, no attempt has been made that all possible errors are fully captured and recovered by the program. Wrong input parameters or an inappropriate use of contructors will often lead to errors that cannot be resolved by the program. While some of this input can be readily recognized as wrong, and then lead to a proper error message, other wrong data may appear dynamically and cannot be captured with a reasonable overhead of the code.  In JAC, therefore, several conditional `if ... elseif ... else ... end` blocks include an additional clause `error("stop a")` or similar; these are clauses, which due to a first design of a function should never be entered, but this appears not to be true in all cases. The use of these (fully) *non-instructive* error message have still a great advantage due to Julia: If not switched-off explicitly, Julia always reports for all program failures the hierarchy of call's, that are made before the error occurs, and lists these calls together with the file and line number of source code. For this reason, an `error("stop a")` readily shows the position where something unexpected occurs. A short inspection of the corresponding (line of the) source code often help to understand of what went wrong internally.

**Julia macros:$\quad$**  What can one do, if the (source) code itself does not tell so much about the problem ? --- In this case, it is often useful to include some additional **printouts**  near to the line in question into the code and to re-run it again. There are different ways (`@show`, `print()`, `println()`) to place printout in the code; cf. https://julialang.org/learning/  A particular quick and useful way makes use of the Julia macro:

```@repl
using JAC   # hide
@show levelsyms
```

which simply repeats the *names* of the variables together with their values. Of course, the values of several such variables can be shown within the same call:

```@repl
wa = 5;   wb = [2.0, pi];   wc = ones(3)
@show wa, wb, wc
```

Indeed, this `@show` macro makes printout very easy. There are many macros (all starting with `@`) in Julia which need not to be considered here. We just mention that `@time` in front of a Julia command (block) will take and display the CPU time that is necessary to run this line(s):

```@repl
@time rand(50000)   ;
```
