### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ f18a9488-13b8-11f0-342b-f99afaf3ac19
begin
	using Pkg
	Pkg.develop(path="/home/fritzsch/fri/JAC.jl")
end

# ╔═╡ fd4ef417-311f-43aa-a601-bf8a9d66b890
using JAC

# ╔═╡ 79374db5-1fbd-4499-b8d9-c7b840bdee52
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(150px, 10%);
    	padding-right: max(150px, 10%);
	}
</style>
"""

# ╔═╡ fcec6737-f8e8-4a8f-91b4-34b8d8a71869
md"""
# Using JAC for the first time
"""

# ╔═╡ 75061d53-4fb7-45a6-b8bb-c199d52479da
md"""
Let us first invoke JAC, either from the own source-code basis or simply by **using JAC**
"""

# ╔═╡ 956d7bd6-d59f-4daa-88b6-9b8879918d79
md"""
We have seen already how the *LiveDocs* in the right-lower corner can be used to access Julia's help features by looking up, for instance, an `atom` or `Orbital`. While an atom itself is not a well-defined term within the JAC toolbox, there exists a number of related terms, such as Atomic, AtomicState (two modules of JAC) as well as several others. If you search for `Orbital`, you will find a `struct  Radial.Orbital`  that defines a tsingle-electron radial orbital function with a large and small component, and which can refer to either the standard or some explicitly given grid due to the logical flag useStandardGrid. More generally, this data type represents a relativistic orbital (function) including additional information that facilitates its use. There are very many (say, more than 300) of such data types (`struct`'s) specified in the JAC toolbox. Since nobody will remember the details of all these definitions, the *LiveDocs* is the right place to remind yourself and make use of these data structures whenever necessary. Special care has been taken that all data structures and functions/methods comes with a reasonable explanation (docstring) in order to work efficiently with JAC. 
"""

# ╔═╡ b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
md"""
Moreover, we can ask the *LiveDocs* of what we can `add`, `extract` or `generate`. The search for `add` show several methods, such as `Basics.add(ma::AngularM64, mb::AngularM64)` to add the projections of two angular momenta ma + mb or `Basics.add(pota::Radial.Potential, potb::Radial.Potential)` to add two radial potentials together, if defined on the same grid. --- Apart from a short explanation, these docstring for the methods always tell the user (i) in which module the given method is defined; (ii) which arguments it takes, including Julia's multiple dispatch feature as well as (iii) the type of the return value. All these information are typically relevant to the user, especially if some input or output does not behave as it should. Indeed, the complexity can grow quite rapidly, for instance, if we search for `generate`, where quite a long list of methods occur.
"""

# ╔═╡ 4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
md"""
**Constructors & program control:** Another frequent use of the *LiveDocs* or Julia's (help) "?" features concerns the data flow and control of almost all computations. In JAC, we often make use of (so-called) `Settings` that enable the user to overwrite default values or to control the computation to the extent, he or she wishes to have control. These `Settings` are context dependent and are different for each atomic property or process that can be computed by the JAC toolbox. They are defined in the various modules and need to be specified accordingly. For instance, to control the computation of transition probabilities for the (fine-structure) levels between given initial- and final-state configuration, one has to overwrite the (defaults) settings: `PhotoEmission.Settings`. Although the photo emission and excitation (absorption) processes are, of course, closely related to each other, the `Settings` will be different since the excitation processes are generally affected also by the properties of the incident radiation; see `PhotoExcitation.Settings`. We shall meet these and (many) other settings quite often in the notebooks below. --- 
Finally, it is sometimes difficult to remember the right term or function name. In this case, it is easy also to make use of the (Unix/Linux) grep command within the JAC/src directly. For those of you who wishes to support and extend the JAC toolbox, the dot expansion and the grep command will be found very helpful, perhaps more than other search tools.
"""

# ╔═╡ 759da271-0adf-4ae7-9642-8620c4c34690
md"""
**Use of constructors:** Another Julia feature, that is frequently applied in JAC, is the successive definition of constructors in order to set-up complex data structures. This features is applied, for instance, in order to define an `Atomic.Computation` or a `Cascade.Computation` as a whole. We shall apply these rather complex data types below in different notebooks. The same issue appears however already for much simpler input. For example, if we wish to select (specify) a number of levels from a multiplet prior to some particular -- configuration interaction -- computation, we can make use of a `LevelSelection`. Apart from the logical flag active, such a level selection requires to either specify a list of level numbers (indices) or level symmetries
"""

# ╔═╡ 644e31bc-d769-4afb-8f81-829a1985a3f3
begin
	LevelSelection(true, indices= [i for i in 1:11])
	LevelSelection(true, symmetries= [LevelSymmetry(1//2, Basics.plus)])	
end

# ╔═╡ 0f7bb544-cd3b-4c5b-8275-18fb09ce035c
md"""
Here, we made use of (the data type) `LevelSymmetry` to specify the overall rotational symmetry of atomic levels. As seen from its definition in the *LiveDocs*, the level symmetry just comprises the total angular momentum (of type `AngularJ64`) and the parity of the level (of type `Parity`). Therefore, the specification of a list of level symmetries in LevelSelection already requires to nest four constructors in order make the level selection explicit: (i) For the angular momentum, (ii) the parity, (iii) the level symmetry and (iv) to create a list (array) of such level symmetries. All the constructors can be specified and built together also in subsequent steps or simply by nesting all the information within a single step, such as:
"""

# ╔═╡ aeaac74e-2087-436c-8bc1-32cec2931f0b
begin
	J1    = AngularJ64(1//2);           J2 = AngularJ64(5//2)  
	pl    = Basics.plus;                mn = Basics.minus
	lsym1 = LevelSymmetry(J1, pl);      lsym2 = LevelSymmetry(J2, mn)
	levelsyms = [lsym1, lsym2]
	levelsyms = [LevelSymmetry(AngularJ64(1//2), Basics.plus), LevelSymmetry(AngularJ64(5//2), Basics.minus)]
end 

# ╔═╡ 610f9dfd-e8df-4706-b750-434ee8f67fb2
md"""
Both ways have their pros and cons, and often some mixture is applied where complex constructors are first assigned to some variables, and which are later utilized to built up constructors of higher complexity. --- To finally specify aǹ instance of a LevelSelection, we use (onc more) its second constructor above, and which will tell the JAC program to compute the lowest three levels (1, 2, 3) as well as all levels with 1/2+ and 5/2- symmetry. Apart from the selection of individual levels, it is often helpful for the computation of atomic processes to make a prior LineSelection and in some cases even a PathwaySelection as, for instance, for dielectronic recombination processes.
"""

# ╔═╡ 993a97ab-ec0a-4bdc-8e61-0bd1692c2524
begin
	LevelSelection(true, indices=[1,2,3], symmetries=levelsyms)
end

# ╔═╡ 4ff08e31-a10f-4aa0-b9a1-aaa0d2996df9
md"""
**Code failures:** Beside of its large flexibility and user-friendliness, JAC might terminate from time to time for non-obvious reasons. Since JAC is first of all a *physics code*, no attempt has been made that all possible errors are fully captured and recovered by the program. Wrong input parameters or an inappropriate use of contructors will often lead to errors that cannot be resolved by the program. While some of this input can be readily recognized as wrong, and then lead to a proper error message, other wrong data may appear dynamically and cannot be captured with a reasonable overhead of the code. In JAC, therefore, several conditional if ... elseif ... else ... end blocks include an additional clause error("stop a") or similar; these are clauses, which due to a first design of a function should never be entered, but this appears not to be true in all cases. The use of these *highly non-instructive* error message have still a great advantage due to Julia: If not switched-off explicitly, Julia always reports for all program failures the hierarchy of call's, that are made before the error occurs, and lists these calls together with the file and line number of source code. For this reason, an error("stop a") readily shows the position where something unexpected occurs. A short inspection of the corresponding (line of the) source code often help to understand of what went wrong internally.
"""

# ╔═╡ 370dd5de-df2d-4f71-827f-90117d1d7e75
md"""
**Julia macros:** What can one do, if the (source) code itself does not tell so much about the problem ? --- In this case, it is often useful to include some additional printouts near to the line in question into the code and to re-run it again. There are different ways (`@show`, `print()`, `println()`) to place printout in the code; cf. https://julialang.org/learning/ A particular quick and useful way makes use of the Julia macro: `@show wa, wb, wc` which simply repeats the names of the variables together with their values. Of course, the values of several such variables can be shown within the same call.  Indeed, this `@show` macro makes printout very easy. There are many macros (all starting with @) in Julia which need not to be considered here. We just mention that `@time` in front of a Julia command (block) will take and display the CPU time that is necessary to run this line(s).
"""

# ╔═╡ cead1e2a-7753-4caa-ab49-ed4c9ccbdd7d
md"""
**Graphics & plots:** One of the powerful features of Julia, finally, is (or could be) its fast and efficient use of graphical representations. In JAC, we wish to make explicit use of such graphical elements, although not too many high-level plot methods are yet prepared. Here, we just wish to briefly remind the user that by `using Plots` or a similar package, one can readily draw (colorful) plots of different dimension. While this is a topic by its own, we make partly use of `Plots` in JAC with the aim to prepare specialized methods that help display the data from complex data structs, such as orbitals, (radial) potentials or synthetic spectra. However, since these plot packages are rather large (and slow package in loading them), they are not loaded automatically but need to called explictly by the user.
"""

# ╔═╡ Cell order:
# ╟─79374db5-1fbd-4499-b8d9-c7b840bdee52
# ╟─fcec6737-f8e8-4a8f-91b4-34b8d8a71869
# ╟─75061d53-4fb7-45a6-b8bb-c199d52479da
# ╟─f18a9488-13b8-11f0-342b-f99afaf3ac19
# ╠═fd4ef417-311f-43aa-a601-bf8a9d66b890
# ╟─956d7bd6-d59f-4daa-88b6-9b8879918d79
# ╟─b0f8586c-db7e-4872-b4e3-ff6beafb4c6c
# ╟─4c2b4d2a-4b34-4f6a-9690-6efc9f2d3fea
# ╟─759da271-0adf-4ae7-9642-8620c4c34690
# ╠═644e31bc-d769-4afb-8f81-829a1985a3f3
# ╟─0f7bb544-cd3b-4c5b-8275-18fb09ce035c
# ╠═aeaac74e-2087-436c-8bc1-32cec2931f0b
# ╟─610f9dfd-e8df-4706-b750-434ee8f67fb2
# ╠═993a97ab-ec0a-4bdc-8e61-0bd1692c2524
# ╟─4ff08e31-a10f-4aa0-b9a1-aaa0d2996df9
# ╟─370dd5de-df2d-4f71-827f-90117d1d7e75
# ╟─cead1e2a-7753-4caa-ab49-ed4c9ccbdd7d
