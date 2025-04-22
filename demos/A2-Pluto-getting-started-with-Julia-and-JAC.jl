### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 9d77f36f-7b62-4f3a-9094-d7c341655220
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(200px, 10%);
    	padding-right: max(200px, 10%);
	}
</style>
"""

# ╔═╡ 9cd6811e-0c7e-4527-8c3f-17dd1da68638
md"""
# Julia Programming Language
!!! quoto "A high-level, dynamic language built for speed and simplicity"
Julia is a high-level, high-performance programming language designed for numerical and scientific computing. The initial design of Julia was started about 2009, and a first beta version was launched around 2012 (well before I got aware of it). Since then, the Julia community has grown very rapidly. While Julia has been developed as a general-purpose language for many different applications, its central features are well suited for high-performance numerical analysis and computational science, and as especially required for modelling quantum many-particle systems. Indeed, Julia comprises a very careful design and a good number of modern technologies that enable the user to gradually learn modern concepts in scientific computing.

Julia offers various features that are otherwise just known from (so-called) productivity languages, including rapid development cycles; exploratory programming without the need to type declarations and memory management; language extensibility via multiple dispatch as well as various features for meta-programming and parallelization. In particular, Julia provides the user with dynamic typing, automatic garbage collection as well as a type-specializing just-in-time compilation of code, and which enables one to port code to other platforms with only moderate adaptions. Julia comes with comprehensive and well-maintained documentation, which can be easily assecible online at <https://docs.julialang.org/en/v1.10/>. With the rapid growth of the Julia user community, there has been developed a large number of tutorials <https://julialang.org/learning/> which help the user to get started. 

Here, we will not cover the detailed instructions for installing Julia. Instead, users are recommended to visit the Julia download page <https://julialang.org/downloads/> for Juliaup at <https://github.com/JuliaLang/juliaup> for more platform-specific instructions. 
!!! info "Installation Tip!!!"
	Users are encouraged to install Julia with **Juliaup**. Detailed instructions to install Julia in different Operating Systems are available at (<https://github.com/JuliaLang/juliaup>). Juliaup gives the flexibility to choose between different versions of Julia which helps to select the correct version of Julia as per the project's compatibility. This [YouTube video](https://www.youtube.com/watch?v=14zfdbzq5BM) explains more on Juliaup installation and its usage.
"""

# ╔═╡ 3c91b372-9ed4-453a-8d19-31ba5e25bbad
md"""
## Jena Atomic Calculator (JAC)
"""

# ╔═╡ 791c73a4-ec96-421c-8233-7b7376ceac9a
md"""
The modern features of Julia are quite in contrast to many codes that are presently (still) applied in computational atomic physics or, more generally, in atomic, molecular and optical (AMO) physics. Mainly for historical reason (and owing to more or less powerful existing codes), most atomic computations make use of Fortran or C. Since atomic structure (and collision) theory has been developed quite in parallel with Fortran, a large (or even huge) amount of Fortran code exists and is applied by the community, a situation with sometimes undesired consequences: Newcomers to computional atomic or AMO either need to learn basic Fortran programming or, vice versa, might belief that (the use of) Fortran and atomic physics are both "old-fashioned". From a physics viewpoint, however, the opposite is true: Until today, atomic pyhsics has been proven -- again and again -- as great playground for developing new ideas, concepts or theoretical methods in quantum many-particle physics. With the design of the **Jena Atomic Calculator (JAC)**, we wish to show that atomic computations can be made simple and applied towards modern and emerging fields of physics.

Here, we shall not introduce Julia's syntax and concepts for which various tutorials are available on the web. Instead, we just wish to remind and highlight some simple (syntax) features that help to go easier around with JAC, and especially for occasional users from experiment or teaching. This reminder aims for lowering the initial threshold for those users who have been trained on other languages in the past. Here, we shall pick up some issues whose physics background is explained later in other tutorial. Obviously, however, Julia is a very rich and powerful language with many features that go well beyond of what is (and will ever be) needed for JAC.
"""

# ╔═╡ 686bebce-9dba-43ab-b02d-f56a9d931a3d
md"""
#### Installation of JAC in Julia REPL / terminal (not in Pluto notebook!!!)

!!! info "With Julia package manger mode"
	- press **]** to enter Julia package manager mode in the Julia terminal (REPL mode)
	- type **add JAC**

!!! info "With 'Pkg.jl' ine Julia REPL"
	using Pkg\
	Pkg.add("JAC")
!!! warning "Not required in Pluto.jl"
	- Installation of JAC is **not explicitly required** in Pluto notebooks which can be simply called with **using JAC**
	- Pluto package manager automatically installs any package that is called/required inside of the notebook by the **using** or **import** keywords; for example, JAC can be called (and also will be installed for the first time) just by **using JAC**
	- More details on Pluto.jl package manager at <https://plutojl.org/en/docs/packages/>
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.9"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╟─9d77f36f-7b62-4f3a-9094-d7c341655220
# ╟─9cd6811e-0c7e-4527-8c3f-17dd1da68638
# ╟─3c91b372-9ed4-453a-8d19-31ba5e25bbad
# ╟─791c73a4-ec96-421c-8233-7b7376ceac9a
# ╠═686bebce-9dba-43ab-b02d-f56a9d931a3d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
