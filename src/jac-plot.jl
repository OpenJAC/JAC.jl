
# using Plots
# pyplot()

using RecipesBase

"""
`Basics.plot()`  ... plots various quantities, often in a new window.

  + `("radial potentials", potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N::Int64 = 0)`  
    ... to plot one or more radial potentials, and where N::Int64 describes the number of grid points to be considered.
        call:  using Plots; pyplot()    ... to access this method by plot(...)
"""
function Basics.plot(sa::String, potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N::Int64 = 0)   
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
end


@recipe function f(sa::String, potentials::Array{Radial.Potential,1}, grid::Radial.Grid; N = 0)
    !(sa == "radial potentials")   &&   error("Unsupported keystring = $sa")
    wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(potentials)
    for  pot in potentials
        wb = wc
        nx = min(length(pot.Zr), N);   wb[1:nx] = pot.Zr[1:nx]
        append!(wa, wb)
        push!(labels, pot.name)
    end
    x = grid.r[1:N];     y = reshape(wa, (N, np))
    label --> reshape(labels, (1,np))
    x, y
end



"""
  + `("radial orbitals: large", orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`  
    ... to plot the large component of one or more radial orbitals, and where N::Int64 describes the number of grid 
        points to be considered.
  + `("radial orbitals: small", orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`  
    ... to plot the small component of one or more radial orbitals.
  + `("radial orbitals: both",  orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N::Int64 = 0)`  
    ... to plot the large and small component of one or more radial orbitals.

    call:  using Plots; pyplot()    ... to access this method by plot(...)
"""
function Basics.plot(sa::String, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...)")
end


@recipe function f(sa::String, orbitals::Array{Radial.Orbital,1}, grid::Radial.Grid; N = 0)
    wa = Float64[];   wc = [NaN for i=1:N];   labels = String[];   np = length(orbitals)

    if       sa == "radial orbitals: large"
        for  orb in orbitals
            wb = wc
            nx = min(length(orb.P), N);   wb[1:nx] = orb.P[1:nx]
            append!(wa, wb)
            push!(labels, "$(orb.subshell):large")
        end
        ##x print("labels = $labels ")
        x = grid.r[1:N];     y = reshape(wa, (N, np))

    elseif   sa == "radial orbitals: small"
        for  orb in orbitals
            wb = wc
            nx = min(length(orb.Q), N);   wb[1:nx] = orb.Q[1:nx]
            append!(wa, wb)
            push!(labels, "$(orb.subshell):small")
        end
        ##x print("labels = $labels ")
        x = grid.r[1:N];     y = reshape(wa, (N, np))

    elseif   sa == "radial orbitals: both"
        for  orb in orbitals
            wb = wc;    nx = min(length(orb.P), N);   wb[1:nx] = orb.P[1:nx];    append!(wa, wb)
            wb = wc;    nx = min(length(orb.Q), N);   wb[1:nx] = orb.Q[1:nx];    append!(wa, wb)
            push!(labels, "$(orb.subshell):large");   push!(labels, "$(orb.subshell):small")
        end
        ##x print("labels = $labels ")
        x = grid.r[1:N];     y = reshape(wa, (N, 2np))

    else   error("Unsupported keystring = $sa") 
    end
    label --> permutedims(labels)
    x, y
end



"""
  + `("spectrum: transition rates over energy", lines::Array{PhotoEmission.Line,1})`  
    ... to plot the transition rates of all lines as function of their transition energies. The plot is shown in a new 
        window but nothing is returned otherwise. **Not yet implemented !**

  + `("spectrum: oscillator strength over energy, emission", lines::Array{PhotoEmission.Line,1})` or
    `("spectrum: oscillator strength over energy, absorption", lines::Array{PhotoEmission.Line,1})` 
    ... to plot the absorption oscillator strength of all lines as function of their transition energies. Again, a new 
        window is opened but nothing returned by this method. **Not yet implemented !**
"""
function Basics.plot(sa::String, lines::Array{PhotoEmission.Line,1})
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...) ... not yet implemented !")
end



"""
  + `("spectrum: transition rates over energy, Gaussian", lines::Array{PhotoEmission.Line,1}; widths=value::Float64)` or
    `("spectrum: transition rates over energy, Lorentzian", lines::Array{PhotoEmission.Line,1}; widths=value::Float64)`
    ... to plot the transition rates of all lines as function of their transition energies but with a Gaussian or Lorentzian 
        distribution. Again, a new window is opened but nothing returned by this method. It still need to be decided how 
        the widths (and, perhaps, other parameters) will be communicated to the method. **Not yet implemented !**
"""
function Basics.plot(sa::String, lines::Array{PhotoEmission.Line,1}, widths::Float64)
    error("call instead:  using Plots; pyplot()    ... to access this method simply by plot(...) ... not yet implemented !")
end

