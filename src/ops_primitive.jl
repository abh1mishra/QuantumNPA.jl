# Abstract type representing primitive operators.
# Specific ones are defined in the file
abstract type Operator end

# Default multiplication rule for operators, which can be specialised.
#
# Assumption for the moment: multiplication returns a pair
#   (c, [ops...])
# consisting of a coefficient and a possibly empty list of operators.
#
# By default we return c = 1 and a list of the same two operators given as
# inputs, i.e., we just concatenate the operators.
Base.:*(x::Operator, y::Operator) = (1, [x, y])

"""
    join_ops(opsx::Vector{Operator}, opsy::Vector{Operator})
This controls how lists of operators are multiplied.
It is not very general at the moment.
Assumption: inputs opsx and opsy both contain at least one element.
Join two lists of operators `opsx` and `opsy`. The function multiplies the last operator
of `opsx` with the first operator of `opsy`, and continues until a non-empty list of operators
is produced or one of the lists is exhausted.

# Arguments
- `opsx::Vector{Operator}`: The first list of operators.
- `opsy::Vector{Operator}`: The second list of operators.

# Returns
- A tuple `(c, ops)` where `c` is the coefficient and `ops` is the resulting list of operators.
"""
function join_ops(opsx::Vector{Operator}, opsy::Vector{Operator})
    j = length(opsx)
    k = 1
    K = 1 + length(opsy)
    c = 1

    while true
        opx, opy = opsx[j], opsy[k]
        (c1, op) = opx * opy

        if c1 == 0
            return (0, Operator[])
        end

        c = rmul(c, c1)
        j -= 1
        k += 1

        if (op != []) || (j == 0) || (k == K)
            ops = Vector{Operator}(vcat(opsx[1:j], op, opsy[k:end]))
            return (c, ops)
        end
    end
end

# Default equality test. This should be specialised for specific types of
# operators (e.g., projectors), so if this default one is called it means the
# arguments are not the same type and are therefore not equal.
Base.:(==)(x::Operator, y::Operator) = false

"""
    Base.isless(x::Operator, y::Operator)
Default order. This again should be specialised for operators of the same
type so here we just see if the type names are ordered lexicographically.
Determine if operator `x` is less than operator `y` based on their type names.

# Arguments
- `x::Operator`: The first operator.
- `y::Operator`: The second operator.

# Returns
- `Bool`: `true` if `x` is less than `y`, `false` otherwise.
"""
function Base.isless(x::Operator, y::Operator)
    return isless(nameof(typeof(x)), nameof(typeof(y)))
end

# This needs to be redefined for operators that aren't Hermitian.
Base.conj(o::Operator) = o

Base.adjoint(o::Operator) = conj(o)

Base.show(io::IO, o::Operator) = print(io, string(o))



"""
    intrace_reduce(ops::Vector{Operator})

Simplify (if possible) the trace of a list of operators `ops`.

# Arguments
- `ops::Vector{Operator}`: The list of operators.
# Details
1. Initialize the coefficient `c` to 1.
2. While the length of `ops` is greater than 2:
   - Multiply the last operator with the first operator.
   - Update the coefficient and the list of operators accordingly.
3. Return the final coefficient `c` and the simplified list of operators `ops`.

# Returns
- A tuple `(c, ops)` where `c` is the coefficient and `ops` is the simplified list of operators.
"""
function intrace_reduce(ops::Vector{Operator})
    c = 1

    while length(ops) > 2
        x = ops[end]
        y = ops[1]
        zs = ops[2:end-1]

        (c1, xy) = x*y

        if length(xy) >= 2
            break
        end

        c = rmul(c, c1)
        ops = Vector{Operator}(vcat(xy, zs))
    end

    return (c, ops)
end

function opcycles(ops::Vector{Operator})
    return (Vector{Operator}(vcat(ops[j:end], ops[1:j-1]))
            for j in 1:length(ops))
end

# function opcycles(ops,subSys::Bool)
#     return reorderMonomial.(Monomial.([vcat(ops[j:end],ops[1:j-1]) for j in 1:length(ops)]))
# end

"""
    opcycles(ops, trace::Bool)

Generate all cyclic permutations of a list of operators `ops` and reorder them.

# Arguments
- `ops`: The list of operators.
- `trace::Bool`: A boolean flag (currently unused).

# Returns
- An array of reordered cyclic permutations of `ops`.
"""
function opcycles(ops, trace::Bool)
    cycleMonArr=[]
    N=length(ops)
    if N==0
        return Id
    end
    for j in 1:N
        ops=vcat(ops[2:end],ops[1])
        cycledOps=reorderMonomial(Monomial(ops))
        if cycledOps == 0
            return 0
        end
        push!(cycleMonArr,cycledOps)
    end
    return cycleMonArr
end


function simplifyOps(ops)
    monend=join_monomials(Monomial([ops[end-1]]),Monomial([ops[end]]))
    if monend isa Tuple
        return nothing

    else
        # might need to flatten the returned monomial
        return (monend==0) ? 0 : vcat(ops[1:end-2],monend.word)
    end

end


"""
    trace(ops::Vector{Operator})

Compute the trace of a list of operators `ops`.

# Arguments
- `ops::Vector{Operator}`: The list of operators.

# Returns
- A tuple `(c, ops)` where `c` is the coefficient and `ops` is the minimum cyclic permutation of the operators.
"""
function trace(ops::Vector{Operator})
    (c, ops) = intrace_reduce(ops)
    return (c, minimum(opcycles(ops)))
end


# Definition of the @operator macro to define operators.

IndexRange = Union{UnitRange{<:Integer},Array{<:Integer}}



"""
    getfields(expr)

Get the fields from an expression `expr`.

# Arguments
- `expr`: The expression.

# Returns
- A list of fields.
"""
getfields(expr) = (expr.head == :tuple) ? expr.args : [expr]

"""
    getfieldnames(fields)

Get the names of the fields.

# Arguments
- `fields`: The list of fields.

# Returns
- A list of field names.
"""
getfieldnames(fields) = map(argname, fields)

"""
    argname(arg)

Get the name of an argument `arg`.

# Arguments
- `arg`: The argument.

# Returns
- The name of the argument.
"""
argname(arg::Symbol) = arg
argname(arg::Expr) = arg.args[1]

"""
    instance_fields(instance, names)

Get the instance fields for a given instance and field names.

# Arguments
- `instance`: The instance.
- `names`: The list of field names.

# Returns
- A list of instance fields.
"""
instance_fields(instance, names) = [:($instance.$f) for f in names]

"""
    fmt_remove(fmt, s::Symbol)

Remove a symbol `s` from a format expression `fmt`.

# Arguments
- `fmt`: The format expression.
- `s::Symbol`: The symbol to remove.

# Returns
- The modified format expression.
"""
function fmt_remove(fmt, s::Symbol)
    return Expr(fmt.head, filter(!isequal(s), fmt.args)...)
end

"""
    parse_fmt(fmt)

Parse a format expression `fmt`.

# Arguments
- `fmt`: The format expression.

# Returns
- A tuple of two format expressions.
"""
function parse_fmt(fmt)
    if fmt.head === :string
        return (fmt, fmt_remove(fmt, :party))
    else
        fmts = fmt.args
        return (fmts[1], fmts[2])
    end
end

"""
    stringfdef(name, fmt)

Generate a string function definition for a given name and format.

# Arguments
- `name`: The name of the operator.
- `fmt`: The format expression.

# Returns
- The string function definition.
"""
function stringfdef(name, fmt)
    fieldnames = [f for f in fmt.args if f isa Symbol]

    args = if (:party in fieldnames)
               (:(x::$name), :(party::Array{Int64,1}))
           else
               (:(x::$name),)
           end

    function fix_special(f)
        if f === :conj
            return :((x.$f ? "*" : ""))
        elseif f === :party
            return :(party_str(party))
        else
            return :(x.$f)
        end
    end

    bindings = (:($f = $(fix_special(f))) for f in fieldnames)

    return :(function Base.string($(args...))
                 $(bindings...)
                 return $fmt
             end)
end

"""
    string_fdefs(name, fmt)

Generate string function definitions for a given name and format.

# Arguments
- `name`: The name of the operator.
- `fmt`: The format expression.

# Returns
- A tuple of two string function definitions.
"""
function string_fdefs(name, fmt)
    (fmt_party, fmt_noparty) = parse_fmt(fmt)
    return (stringfdef(name, fmt_party), stringfdef(name, fmt_noparty))
end

"""
    conjfalse(field)

Set the `conj` field to false if it exists.

# Arguments
- `field`: The field.

# Returns
- The modified field.
"""
function conjfalse(field)
    return ((argname(field) === :conj) ? Expr(:kw, field, false) : field)
end

"""
    conj_def(name, fieldnames)

Generate a `conj` function definition for a given name and field names.

# Arguments
- `name`: The name of the operator.
- `fieldnames`: The list of field names.

# Returns
- The `conj` function definition.
"""
function conj_def(name, fieldnames)
    if :conj in fieldnames
        cxfields = replace(instance_fields(:x, fieldnames),
                           :(x.conj) => :(!x.conj))
        conjf = :( Base.conj(x::$name) = $name($(cxfields...)) )
    else
        conjf = nothing
    end

    return conjf
end

"""
    parse_ctoropt(s::Bool)

Parse a constructor option.

# Arguments
- `s::Bool`: The constructor option.

# Returns
- A tuple of the constructor option.
"""
parse_ctoropt(s::Bool) = (s, s)

"""
    parse_ctoropt(x::Expr)

Parse a constructor option.

# Arguments
- `x::Expr`: The constructor option.

# Returns
- A tuple of the constructor option arguments.
"""
parse_ctoropt(x::Expr) = x.args

"""
    int_args(fields)

Get the integer arguments from a list of fields.

# Arguments
- `fields`: The list of fields.

# Returns
- A list of integer arguments.
"""
function int_args(fields)
    return filter(fields) do f
        !isa(f, Symbol) && (f.args[2] === :Integer)
    end
end

"""
    replace_seq(collection, replacements)

Replace elements in a collection with replacements.

# Arguments
- `collection`: The collection.
- `replacements`: The replacements.

# Returns
- The modified collection.
"""
function replace_seq(collection, replacements)
    collection = copy(collection)

    for r in replacements
        replace!(collection, r)
    end

    return collection
end

"""
    chtype(arg::Expr, type::Symbol)

Change the type of an argument.

# Arguments
- `arg::Expr`: The argument.
- `type::Symbol`: The new type.

# Returns
- The modified argument.
"""
chtype(arg::Expr, type::Symbol) = Expr(arg.head, arg.args[1], type)

"""
    chtype(args::Array{Expr,1}, type::Symbol)

Change the type of a list of arguments.

# Arguments
- `args::Array{Expr,1}`: The list of arguments.
- `type::Symbol`: The new type.

# Returns
- The modified list of arguments.
"""
chtype(args::Array{Expr,1}, type::Symbol) = [chtype(a, type) for a in args]

"""
    chtypes(fields, args::Array{Expr,1}, type::Symbol)

Change the types of fields based on a list of arguments.

# Arguments
- `fields`: The list of fields.
- `args::Array{Expr,1}`: The list of arguments.
- `type::Symbol`: The new type.

# Returns
- The modified list of fields.
"""
function chtypes(fields, args::Array{Expr,1}, type::Symbol)
    return replace_seq(fields, (a => chtype(a, type) for a in args))
end

"""
    ctor_ranges(lcname, fields, fieldnames)

Generate constructor ranges for a given name, fields, and field names.

# Arguments
- `lcname`: The lowercase name of the operator.
- `fields`: The list of fields.
- `fieldnames`: The list of field names.

# Returns
- A list of constructor range functions.
"""
function ctor_ranges(lcname, fields, fieldnames)
    ifields = int_args(fields)
    to_sub = drop(powerset(ifields), 1)

    function mkrange(sub)
        fnames = getfieldnames(sub)
        nfields = chtypes(fields, sub, :IndexRange)
        lfcall = :($lcname(party, $(fieldnames...)))
        lassgms = [:($a = $a) for a in fnames]
        comp = Expr(:comprehension, Expr(:generator, lfcall, lassgms...))

        return :(function $lcname(party, $(nfields...))
                     return $comp
                 end)
    end

    return map(mkrange, to_sub)
end

"""
    constructor_defs(make_constructors, name, fields, fieldnames)

Generate constructor definitions for a given name, fields, and field names.

# Arguments
- `make_constructors`: A boolean or expression indicating whether to make constructors.
- `name`: The name of the operator.
- `fields`: The list of fields.
- `fieldnames`: The list of field names.

# Returns
- A tuple of the main constructor and constructor range functions.
"""
function constructor_defs(make_constructors, name, fields, fieldnames)
    lcname = Symbol(lowercase(string(name)))
    (mk_ctor, mk_crange) = parse_ctoropt(make_constructors)

    cargs = Any[conjfalse(f) for f in fields]

    if mk_ctor
        mctor = :(function $(esc(lcname))(party, $(cargs...))
                      return Monomial(party, $name($(fieldnames...)))
                  end)
    else
        mctor = nothing
    end

    if mk_crange
        mcrange = ctor_ranges(esc(lcname), cargs, fieldnames)
    else
        mcrange = ()
    end

    return (mctor, mcrange)
end

"""
    macro operator(ctor::Expr, fmt, make_constructors=true, order=nothing)

Define a new type of operator with a given name (e.g., Projector), fields,
and format string. In addition to generating the struct definition this also
generates method definitions for the following generic functions:

  * Base.hash,
  * Base.:(==),
  * Base.isless,
  * Base.string,

as well as a constructor method with the name in lowercase (e.g., projector)
that creates a monomial containing a single operator associated to a given
party.

If one of the fields is named conj a Base.conj method is also generated.

# Arguments
- `ctor::Expr`: The constructor expression.
- `fmt`: The format string.
- `make_constructors`: A boolean indicating whether to make constructors (default: true).
- `order`: The order of the fields (default: nothing).

# Returns
- The generated operator type and associated methods.
"""
macro operator(ctor::Expr, fmt, make_constructors=true, order=nothing)
    name = ctor.args[1]
    fields = ctor.args[2:end]

    fieldnames = getfieldnames(fields)

    order = (isnothing(order) ? reverse(fieldnames) : getfields(order))
    xfields = :($(instance_fields(:x, order)...),)
    yfields = :($(instance_fields(:y, order)...),)

    (strf_party, strf_noparty) = string_fdefs(name, fmt)

    conjf = conj_def(name, fieldnames)
    (mctor, mcrange) = constructor_defs(make_constructors,
                                        name,
                                        fields,
                                        fieldnames)

    methods = [strf_party, strf_noparty, conjf, mctor, mcrange...]
    methods = filter(!isnothing, methods)

    return quote
        struct $name <: Operator
            $(fields...)
        end
        Base.hash(x::$name, h::UInt) = hash($xfields, h)
        Base.:(==)(x::$name, y::$name) = ($xfields == $yfields)
        Base.isless(x::$name, y::$name) = ($xfields < $yfields)
        $(methods...)
    end
end