using MacroTools: @capture, postwalk

using ..Basics: ⊕, projections

# If no upper or lower bound for a summation index can be determined, this cutoff is used,
# i.e., the range of summation indices is considered a subset of -CUTOFF:1/2:CUTOFF.
const CUTOFF = 30

"""
    @racahsum expression var1[=range1] var2[=range2] …

Evaluate the sum over `expression` with the summation indices `var1`, `var2`, etc. The
summand `expression` must contain at least one call to one of the functions
[`clebschgordan`](@ref), [`wigner3j`](@ref), [`wigner6j`](@ref), and [`wigner9j`](@ref). The
bounds for each of the specified summation indices are determined from the selection rules
for the arguments of these functions.

The summation index may be restricted to a range of values (which may be reduced further by
applying the selection rules). The type of such a range should be a subtype of
`OrdinalRange{<:Union{HalfInt,Integer}}`, e.g. `UnitRange{Int}` or
`StepRange{HalfInt,HalfInt}`.

# Examples

A typical call to this macro might look like this:
```jldoctest
julia> @racahsum (2j₃+1)*wigner6j(1,3/2,j₃,0,5/2,3/2)*wigner6j(1,3/2,j₃,0,5/2,3/2) j₃
0.24999999999999986
```
Here, we have verified the orthogonality relation for the Wigner 6-j symbol by summing over
the index `j₃`, which according to the selection rules can only take the value `5/2`, so the
sum consists of only one summand.

To verify the orthogonality relation of Wigner 9-j symbols, two summation indices are needed:
```jldoctest
julia> @racahsum (2j₇+1)*(2j₈+1)*wigner9j(1/2,1/2,1,3,3,2,j₇,j₈,1)*wigner9j(1/2,1/2,1,3,3,2,j₇,j₈,1) j₇ j₈
0.06666666666666662
```

If no upper bound can be placed on a summation index, it is set to $CUTOFF (unless a range is
specified explicitly):
```jldoctest
julia> @racahsum wigner6j(j,0,j,0,j,0) j       # j ∈ 0:1/2:$CUTOFF
0.6686547347577827

julia> @racahsum wigner6j(j,0,j,0,j,0) j=1:10  # j ∈ 0:10
3.2630912526125706
```
"""
macro racahsum(ex::Expr, sumindices...)
    # The following arrays hold all information about the summation indices:
    # * sumvars: the names of the variables
    # * sumtypes: the types of the variables, determined from an explicit type annotation or
    #   from their position in the Wigner symbols or Clebsch–Gordan coefficients. Possible
    #   values are :AngularL, :AngularJ, :AngularM, and :undef (if no type has been
    #   determined yet). The types are only used to restrict the possible range of each
    #   summation index (i.e., an AngularM may be an arbitrary half-integer, an AngularJ
    #   must be a non-negative half-integer, and an AngularL must be a non-negative
    #   integer), the actual type of a summation index is always HalfInt. If a symbol
    #   appears in places for different types (e.g., in an AngularM field of a 3-j symbol
    #   and in a 6-j symbol), the most restrictive type is chosen (AngularL ⊂ AngularJ ⊂
    #   AngularM).
    # * sumranges: expressions that evaluate to ranges which restrict the summation indices,
    #   e.g. `:(projections(j))` or `:(j₁ ⊕ j₂)`.
    # * sumvalues: expressions that evaluate to values that the summation index must be
    #   equal to, e.g. `:(-m₁-m₂)` for `m₃` in a Wigner 3-j symbol.
    # * sumlbounds: expressions that evaluate lower bounds on summation indices, i.e.,
    #   `:(abs(m₁))` implies that the summation index may only take values abs(m₁):1:∞.
    # Any expression that restricts a summation index (in sumranges, sumvalues, or
    # sumlbounds) may only use the values of outermore summation indices (i.e., the ones
    # that were declared earlier in the list).
    sumvars  = Vector{Symbol}(undef, length(sumindices))
    sumtypes = Vector{Symbol}(undef, length(sumindices))
    sumranges  = [[] for _ = sumindices]
    sumvalues  = [[] for _ = sumindices]
    sumlbounds = [[] for _ = sumindices]
    # Parse list of summation indices (in reverse order, so innermost index comes first)
    # with optional type annotation and range
    for (i,sumvar) in enumerate(Iterators.reverse(sumindices))
        var, typ, range = _parse_sumvar(sumvar)
        var in sumvars[1:i-1] && error("summation index $var declared multiple times.")
        sumvars[i] = var
        sumtypes[i] = something(typ, :undef)
        if range !== nothing
            _in_sumvars(range, sumvars[1:i]) &&
                error("summation range contains innermore summation index: $sumvar.")
            push!(sumranges[i], esc(range))
        end
    end
    # Analyse “Racah expression” to determine bounds for each summation index
    postwalk(ex) do x
        if @capture(x, clebschgordan(args__))
            _clebschgordanranges!(sumranges, sumtypes, sumvalues, sumlbounds, sumvars, args)
        elseif @capture(x, wigner3j(args__))
            _wigner3jranges!(sumranges, sumtypes, sumvalues, sumlbounds, sumvars, args)
        elseif @capture(x, wigner6j(args__))
            _wigner6jranges!(sumranges, sumtypes, sumvars, args)
        elseif @capture(x, wigner9j(args__))
            _wigner9jranges!(sumranges, sumtypes, sumvars, args)
        end
        x
    end
    # Assemble loops and return expression. The summation index that was declared last
    # (i.e., the *first* variable in sumvars) defines the innermost for loop or let block.
    loops = :(sum += $(esc(ex)))
    for (var,typ,ranges,values,lbounds) in zip(sumvars,sumtypes,sumranges,sumvalues,sumlbounds)
        typ === :undef && error("$var does not appear in any Wigner symbol.")
        loops = _addloop(loops, var, typ, ranges, values, lbounds)
    end
    :(local sum = 0.0; $loops; sum)
end

# Return a for loop with the loop body `body`, looping over the variable `var` which is
# restricted by the type `typ` (:AngularJ, :AngularL, or :AngularM), the ranges `ranges` and
# the lower bounds `lbounds`. If `values` is not empty, a let block is created instead of a
# loop, setting `var` to the first value in `values` and executing `body` only if `var` is
# equal to all other elements of `values` and is contained in each range specified by
# `ranges` and `lbounds`.
function _addloop(body, var, typ, ranges, values, lbounds)
    allranges = copy(ranges)
    _append_lbounds!(allranges, lbounds, isempty(ranges) & isempty(values))
    if isempty(values)
        _append_typerange!(allranges, typ, !isempty(lbounds) || _hasoplus(ranges))
        if isempty(allranges)
            error("failed to determine summation bounds for $var.")
        elseif length(allranges) == 1
            ex = :(for $(esc(var)) = StepRange{HalfInt,HalfInt}($(allranges[1])); $body; end)
        else
            ex = :(for $(esc(var)) = StepRange{HalfInt,HalfInt}(∩($(allranges...))); $body; end)
        end
    else
        conds = [:($(esc(var)) == $val) for val = values[2:end]]
        if isempty(lbounds) && !_hasoplus(ranges) && typ ∈ (:AngularJ,:AngularL)
            push!(conds, :($(esc(var)) ≥ 0))
        end
        if typ === :AngularL
            push!(conds, :(isinteger($(esc(var)))))
        end
        append!(conds, [:($(esc(var)) ∈ $range) for range = allranges])
        condition = :(true)
        for cond in Iterators.reverse(conds)
            condition = :($cond && $condition)
        end
        ex = :(let $(esc(var)) = HalfInt($(values[1])); if $condition; $body; end; end)
    end
    ex
end

# Append lbounds to ranges, using lbound:CUTOFF or lbound:typemax(HalfInt)
# depending on the value of usecutoff
function _append_lbounds!(ranges, lbounds, usecutoff)
    for lbound in lbounds
        if usecutoff
            push!(ranges, :($lbound:$CUTOFF))
        else
            push!(ranges, :($lbound:$(typemax(HalfInt))))
        end
    end
end

# Append boundaries based on the type (:AngularJ, :AngularL, or :AngularM) to ranges.
# If `isnonnegative` is true, the ranges are assumed to already exclude negative numbers.
function _append_typerange!(ranges, typ, isnonnegative)
    if typ === :AngularJ
        if isempty(ranges)
            push!(ranges, 0:HalfInt(1,2):CUTOFF)
        elseif !isnonnegative
            push!(ranges, 0:HalfInt(1,2):typemax(HalfInt))
        end
    elseif typ === :AngularL
        isempty(ranges) ? push!(ranges, 0:CUTOFF) : push!(ranges, 0:typemax(HalfInt))
    elseif typ === :AngularM
        isempty(ranges) && push!(ranges, -CUTOFF:onehalf(HalfInt):CUTOFF)
    else
        error("invalid type: $typ.")
    end
end

# Check whether any of the expressions in `ranges` is an angular momentum addition `a ⊕ b`,
# where `a` and `b` may be any other object.
_hasoplus(ranges) = any(x -> @capture(x, a_ ⊕ b_), ranges)

# For an expression `var::Type=value` where `::Type` and `=value` are optional, return the
# tuple (var, Type, value). If `Type` and/or `value` are not given, `nothing` is returned in
# their place.
function _parse_sumvar(expr)
    if @capture(expr, var_=range_)
        _extract_sumvartype(var)..., range
    else
        _extract_sumvartype(expr)..., nothing
    end
end
_parse_sumvar(x::Symbol) = x, nothing, nothing

# For an expression `var::Type` where `::Type` is optional, return the tuple (var, Type). If
# `Type` is not given, `nothing` is returned in its place.
function _extract_sumvartype(expr)
    if @capture(expr, var_::T_) && var isa Symbol && T isa Symbol
        T in (:AngularL, :AngularJ, :AngularM) ||
            error("invalid type annotation in summation index: $expr.")
        return var, T
    end
    error("invalid summation index expression: $expr.")
end
_extract_sumvartype(var::Symbol) = var, nothing

# Return true iff the expression `a` contains any of the symbols in `otherjs`
_in_sumvars(a, otherjs) = a ∈ otherjs
_in_sumvars(a::Expr, otherjs) =
    _in_sumvars(a.head, otherjs) || any(x -> _in_sumvars(x, otherjs), a.args)

# For a clebschgordan call with the arguments `args`, extract the restrictions that can be
# placed on each variable in `sumvars` due to selection rules and append the resulting
# ranges/types/values/lbounds to the appropriate arrays. Only restrictions that do not
# depend on innermore summation indices (i.e., symbols that appear earlier in `sumvars`) or
# the summation index itself are considered.
function _clebschgordanranges!(sumranges, sumtypes, sumvalues, sumlbounds, sumvars, args)
    _check_argcount(:clebschgordan, 6, args)
    for ind = 1:2:5
        j = args[ind]
        m = args[ind+1]
        jaindex, jbindex = mod1(ind+2, 6), mod1(ind+4, 6)
        jindex = findfirst(isequal(j), sumvars)
        if jindex !== nothing
            _setvartype!(sumtypes, jindex, :AngularJ)
            _addjrange!(sumranges[jindex], sumvars[1:jindex], args[jaindex], args[jbindex])
            _addlbound!(sumlbounds[jindex], sumvars[1:jindex], m)
        end
    end
    j₁, m₁, j₂, m₂, J, M = args
    if @capture(m₁, -mm₁_)
        if (mindex=findfirst(isequal(mm₁),sumvars); mindex !== nothing)
            _setvartype!(sumtypes, mindex, :AngularM)
            _addmrange!(sumranges[mindex], sumvars[1:mindex], j₁)
            _addvalue_pm!(sumvalues[mindex], sumvars[1:mindex], m₂, M)
        end
    elseif (mindex=findfirst(isequal(m₁),sumvars); mindex !== nothing)
        _setvartype!(sumtypes, mindex, :AngularM)
        _addmrange!(sumranges[mindex], sumvars[1:mindex], j₁)
        _addvalue_pm!(sumvalues[mindex], sumvars[1:mindex], M, m₂)
    end
    if @capture(m₂, -mm₂_)
        if (mindex=findfirst(isequal(mm₂),sumvars); mindex !== nothing)
            _setvartype!(sumtypes, mindex, :AngularM)
            _addmrange!(sumranges[mindex], sumvars[1:mindex], j₂)
            _addvalue_pm!(sumvalues[mindex], sumvars[1:mindex], m₁, M)
        end
    elseif (mindex=findfirst(isequal(m₂),sumvars); mindex !== nothing)
        _setvartype!(sumtypes, mindex, :AngularM)
        _addmrange!(sumranges[mindex], sumvars[1:mindex], j₂)
        _addvalue_pm!(sumvalues[mindex], sumvars[1:mindex], M, m₁)
    end
    if @capture(M, -mM_)
        if (mindex=findfirst(isequal(mM),sumvars); mindex !== nothing)
            _setvartype!(sumtypes, mindex, :AngularM)
            _addmrange!(sumranges[mindex], sumvars[1:mindex], J)
            _addvalue_mm!(sumvalues[mindex], sumvars[1:mindex], m₁, m₂)
        end
    elseif (mindex=findfirst(isequal(M),sumvars); mindex !== nothing)
        _setvartype!(sumtypes, mindex, :AngularM)
        _addmrange!(sumranges[mindex], sumvars[1:mindex], J)
        _addvalue_pp!(sumvalues[mindex], sumvars[1:mindex], m₁, m₂)
    end
end

# Like _clebschgordanranges!, but for wigner3j
function _wigner3jranges!(sumranges, sumtypes, sumvalues, sumlbounds, sumvars, args)
    _check_argcount(:wigner3j, 6, args)
    for ind = 1:3
        j = args[ind]
        m = args[ind+3]
        jaindex, jbindex = mod1(ind+1, 3), mod1(ind+2, 3)
        maindex, mbindex = jaindex+3, jbindex+3
        if (jindex=findfirst(isequal(j),sumvars); jindex !== nothing)
            _setvartype!(sumtypes, jindex, :AngularJ)
            _addjrange!(sumranges[jindex], sumvars[1:jindex], args[jaindex], args[jbindex])
            _addlbound!(sumlbounds[jindex], sumvars[1:jindex], m)
        end
        if @capture(m, -mm_)
            if (mindex=findfirst(isequal(mm),sumvars); mindex !== nothing)
                _setvartype!(sumtypes, mindex, :AngularM)
                _addmrange!(sumranges[mindex], sumvars[1:mindex], j)
                _addvalue_pp!(sumvalues[mindex], sumvars[1:mindex], args[maindex], args[mbindex])
            end
        elseif (mindex=findfirst(isequal(m),sumvars); mindex !== nothing)
            _setvartype!(sumtypes, mindex, :AngularM)
            _addmrange!(sumranges[mindex], sumvars[1:mindex], j)
            _addvalue_mm!(sumvalues[mindex], sumvars[1:mindex], args[maindex], args[mbindex])
        end
    end
end

# For a wigner6j call with the arguments `args`, extract the restrictions that can be placed
# on each variable in `sumvars` due to selection rules and append the resulting ranges/types
# to the appropriate arrays. Only restrictions that do not depend on innermore summation
# indices (i.e., symbols that appear earlier in `sumvars`) or the summation index itself are
# considered.
function _wigner6jranges!(sumranges, sumtypes, sumvars, js)
    _check_argcount(:wigner6j, 6, js)
    for (ind,j) = enumerate(js)
        jindex = findfirst(isequal(j), sumvars)
        jindex === nothing && continue
        _setvartype!(sumtypes, jindex, :AngularJ)
        div, rem = divrem(ind-1, 3)
        aindex, bindex = mod1(rem+2, 3), mod1(rem, 3)
        cindex, dindex = aindex+3, bindex+3
        if iszero(div) # j is in the top row of the 6-j symbol
            _addjrange!(sumranges[jindex], sumvars[1:jindex], js[aindex], js[bindex])
            _addjrange!(sumranges[jindex], sumvars[1:jindex], js[cindex], js[dindex])
        else # j is in the bottom row of the 6-j symbol
            _addjrange!(sumranges[jindex], sumvars[1:jindex], js[aindex], js[dindex])
            _addjrange!(sumranges[jindex], sumvars[1:jindex], js[cindex], js[bindex])
        end
    end
end

# Same as _wigner6jranges! but for wigner9j
function _wigner9jranges!(sumranges, sumtypes, sumvars, js)
    _check_argcount(:wigner9j, 9, js)
    for (ind,j) = enumerate(js)
        jindex = findfirst(isequal(j), sumvars)
        jindex === nothing && continue
        _setvartype!(sumtypes, jindex, :AngularJ)
        div, rem = divrem(ind-1, 3)
        aindex = 3div + mod1(rem+2, 3)
        bindex = 3div + mod1(rem, 3)
        cindex = 3mod(div+1, 3) + rem + 1
        dindex = 3mod(div+2, 3) + rem + 1
        _addjrange!(sumranges[jindex], sumvars[1:jindex], js[aindex], js[bindex])
        _addjrange!(sumranges[jindex], sumvars[1:jindex], js[cindex], js[dindex])
    end
end

function _check_argcount(func::Symbol, expected, args)
    length(args) == expected || error("wrong number of arguments for $func: " *
                                      "expected $expected, got $(length(args)).")
end

# Set `types[index] = newtype`, but only if `newtype` is more restrictive than the type that
# is currently stored in `types[index]`.
function _setvartype!(types, index, newtype)
    order = (:undef, :AngularM, :AngularJ, :AngularL)
    oldindex = findfirst(isequal(types[index]), order)
    newindex = findfirst(isequal(newtype), order)
    oldindex === nothing && error("invalid type annotation: $(types[index]).")
    newindex === nothing && error("invalid type annotation: $newtype.")
    types[index] = order[max(oldindex,newindex)]
end

# Add the range `j₁ ⊕ j₂` to `jranges`, unless it is already in there or contains
# symbols from `otherjs`.
function _addjrange!(jranges, otherjs, j₁, j₂)
    if !_in_sumvars(j₁, otherjs) && !_in_sumvars(j₂, otherjs) &&
        (:($(esc(j₁)) ⊕ $(esc(j₂))) ∉ jranges) && (:($(esc(j₂)) ⊕ $(esc(j₁))) ∉ jranges)
        push!(jranges, :($(esc(j₁)) ⊕ $(esc(j₂))))
    end
end

# Add the lower bound `abs(m)` to `jlbounds`, unless it is already in there or contains
# symbols from `otherjs`.
function _addlbound!(jlbounds, otherjs, m)
    if !_in_sumvars(m, otherjs) && (:(abs(HalfInt($(esc(m))))) ∉ jlbounds)
        push!(jlbounds, :(abs(HalfInt($(esc(m))))))
    end
end

# Add the range `projections(j)` (i.e., `-j:j`) to `mranges`, unless it is already in there
# or contains symbols from `otherjs`.
function _addmrange!(mranges, otherjs, j)
    if !_in_sumvars(j, otherjs) && (:(projections($(esc(j)))) ∉ mranges)
        push!(mranges, :(projections($(esc(j)))))
    end
end

# Add the value `-m₁-m₂` to `mvalues`, unless it is already in there or contains symbols
# from `otherjs`.
function _addvalue_mm!(mvalues, otherjs, m₁, m₂)
    if !_in_sumvars(m₁, otherjs) && !_in_sumvars(m₂, otherjs) &&
        (:(-$(esc(m₁))-$(esc(m₂))) ∉ mvalues) && (:(-$(esc(m₂))-$(esc(m₁))) ∉ mvalues)
        push!(mvalues, :(-$(esc(m₁))-$(esc(m₂))))
    end
end

# Add the value `M-m` to `mvalues`, unless it is already in there or contains symbols
# from `otherjs`.
function _addvalue_pm!(mvalues, otherjs, M, m)
    if !_in_sumvars(M, otherjs) && !_in_sumvars(m, otherjs) && (:($(esc(M))-$(esc(m))) ∉ mvalues)
        push!(mvalues, :($(esc(M))-$(esc(m))))
    end
end

# Add the value `m₁+m₂` to `mvalues`, unless it is already in there or contains symbols
# from `otherjs`.
function _addvalue_pp!(mvalues, otherjs, m₁, m₂)
    if !_in_sumvars(m₁, otherjs) && !_in_sumvars(m₂, otherjs) &&
        (:($(esc(m₁))+$(esc(m₂))) ∉ mvalues) && (:($(esc(m₂))+$(esc(m₁))) ∉ mvalues)
        push!(mvalues, :($(esc(m₁))+$(esc(m₂))))
    end
end
