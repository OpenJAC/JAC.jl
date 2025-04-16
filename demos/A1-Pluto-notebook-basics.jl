### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 76fd22f7-01b0-452a-a785-2d385da856a8
md"""
# Introduction to Pluto.jl notebook
"""

# ╔═╡ 55cb5bf8-b26c-4410-9e1f-0acbb6d1cf84
md"""
!!! info "Attention!!!"
	- This particular Pluto notebook provides some details how Pluto can be utilized in order to explore the basic capablities of Jena Atomic Calculator (JAC).
	- A more comprehensive documentation for Pluto.jl is available at **<https://plutojl.org/en/docs/>**
	- Users are highly recommended to go through the official documentation for further details on the installation and usage of Pluto.jl for Julia
"""

# ╔═╡ beb799bd-d87e-4641-bbe8-130f458456fb
md"""
#### Installation of Pluto.jl and Usages
- Pluto.jl can be installed analogue to all other Julia packages
- <https://plutojl.org/#install> provides further details on the installation of Pluto.jl

"""

# ╔═╡ 74b89f11-822c-48ea-91cd-edf6ba26c096
md"""
#### Accessing the documentation in the Pluto notebook
- The documentation or doc strings about any (Julia) function can be accessed inside of a Pluto notebook by using the **Live Docs** at the **bottom-right** corner of the windows; this feature replaces the "?" at Julia's terminal/REPL.
#### Giving input in a cell
- Single-line inputs can be given in a (Pluto) cell as usually given at the Julia REPL.
- Multi-line inputs in a single cell need to be wrapped with **"begin ... end"**
#### Executing a cell
- a selected cell can executed by either clicking the **run cell** button at the right bottom of the cell or just by (Shift + Enter)
- after the execution has been completed, the output of a cell is shown **on top of the cell**, while all printed output during execution (the "logs") become visible below of the cell. For example:
"""

# ╔═╡ d34976ee-e73c-497e-a36f-215e97f59594
begin
	a = 4 + 3
	println("The sum of 4 and 3 is: ", a )
	a
end

# ╔═╡ 6ade9c7f-8703-47f6-94aa-7600c3743fe7
md"""
#### Calling a package inside of a Pluto notebook
The package manager in Pluto works a bit differently when compared to Julia's package manager. If the package called is not installed in Julia, then the Pluto package manager installs the called package automatically. More details on Pluto's package manager can be found at <https://plutojl.org/en/docs/packages/>
"""

# ╔═╡ 57a590b0-c27c-4eca-b9ab-7c0212cdd074
md"""
!!! info "Better Visualization in Pluto notebook!!!"
	By default, Pluto notebooks open with rather a narrow width, which is often insufficient to visualize the output or to write longer inputs.

	This can be easily handled with a workaround to increase the width of each cell by including the following **html** input in any cell of the notebook (preferably, however, as the first cell).
"""

# ╔═╡ 489d1501-b32e-4b4f-a8d7-decfa97c4bc9
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

# ╔═╡ Cell order:
# ╟─76fd22f7-01b0-452a-a785-2d385da856a8
# ╟─55cb5bf8-b26c-4410-9e1f-0acbb6d1cf84
# ╟─beb799bd-d87e-4641-bbe8-130f458456fb
# ╟─74b89f11-822c-48ea-91cd-edf6ba26c096
# ╠═d34976ee-e73c-497e-a36f-215e97f59594
# ╟─6ade9c7f-8703-47f6-94aa-7600c3743fe7
# ╟─57a590b0-c27c-4eca-b9ab-7c0212cdd074
# ╠═489d1501-b32e-4b4f-a8d7-decfa97c4bc9
