struct PPolynomial
    terms::Dict{PMonomial,Number}
end

PPolynomial() = PPolynomial(Dict{PMonomial,Number}())

PPolynomial(x::Number) = PPolynomial((x != 0) ? Dict(Id => demote(x)) : Dict())

PPolynomial(x::PMonomial) = PPolynomial(Dict(x => 1))

function PPolynomial(x::Number, y::PMonomial)
    return (x != 0) ? PPolynomial(Dict(y => demote(x))) : PPolynomial()
end

PPolynomial(x::PPolynomial) = x

function PPolynomial(x::Base.Generator)
    return PPolynomial(Dict((m, demote(c)) for (c, m) in x))
end

new_ppolynomial(x) = PPolynomial(x)
new_ppolynomial(p::PPolynomial) = copy(p)



# Iteration over Polynomials.

function swap_mc(result::Tuple{Pair{PMonomial, Number}, Int64})
    if isnothing(result)
        return nothing
    else
        ((m, c), state) = result
        return (Pair{Number,PMonomial}(c, m), state)
    end
end

Base.iterate(x::PPolynomial) = swap_mc(iterate(x.terms))

Base.iterate(x::PPolynomial, state) = swap_mc(iterate(x.terms, state))



Base.hash(p::PPolynomial, h::UInt) = hash(p.terms, h)

Base.getindex(x::PPolynomial, y::PMonomial) = get(x.terms, y, 0)

function Base.setindex!(x::PPolynomial, y::Number, z::PMonomial)
    if y == 0
        delete!(x.terms, z)
    else
        x.terms[z] = demote(y)
    end
end

function Base.copy(x::PPolynomial)
    return PPolynomial(copy(x.terms))
end



function Base.show(io::IO, p::PPolynomial)
    if isempty(p)
        print(io, "0")
    else
        c2s = firstcoeff2string

        for (c, m) in (p)
            print(io,c2s(c))
            show(io, m)
            c2s = coeff2string
        end
    end
end

"""
Degree of a PPolynomial.
degree(0) returns negative infinity. With this definition the following rules
hold even if one or both of P or Q are zero:
  degree(P + Q) == max(degree(P), degree(Q))
  degree(P - Q) <= max(degree(P), degree(Q))
  degree(P * Q) == degree(P) + degree(Q)
"""
function degree(p::PPolynomial)
    return !iszero(p) ? maximum(degree.(monomials(p))) : -Inf
end

degree_less(n::Integer) = (x -> (degree(x) < n))



"Add y to PPolynomial x, modifying x."
function add!(x::PPolynomial, y::Number)
    x[Id] += y
    return x
end

function add!(x::PPolynomial, y::PMonomial)
    x[y] += 1
    return x
end

function add!(x::PPolynomial, y::PPolynomial)
    for (c, m) in y
        x[m] += c
    end

    return x
end



"Add y*z to the PPolynomial x, modifying x. y has to be a number."
function addmul!(x::PPolynomial, y::Number, z::Number)
    x[Id] += y*z
    return x
end

function addmul!(x::PPolynomial, y::Number, z::PMonomial)
    x[z] += y
    return x
end

function addmul!(x::PPolynomial, y::Number, z::PPolynomial)
    for (c, m) in z
        x[m] += c*y
    end

    return x
end



"Subtract y from PPolynomial x."
function sub!(x::PPolynomial, y::Number)
    x[Id] -= y
    return x
end

function sub!(x::PPolynomial, y::PMonomial)
    x[y] -= 1
    return x
end

function sub!(x::PPolynomial, y::PPolynomial)
    for (c, m) in y
        x[m] -= c
    end

    return x
end



function mul!(x::PPolynomial, y::Number)
    for (c, m) in x
        x[m] = c*y
    end

    return x
end



"Return a PPolynomial consisting of the sum of items in s."
function psum(s)
    z = PPolynomial()

    for x in s
        add!(z, x)
    end

    return z
end

psum(m::PMonomial) = PPolynomial(m)



Base.:+(x::Number, y::PMonomial) = add!(PPolynomial(y), x)
Base.:+(x::PMonomial, y::Number) = add!(PPolynomial(x), y)

function Base.:+(x::PMonomial, y::PMonomial)
    return PPolynomial((x != y) ? Dict(x => 1, y => 1) : Dict(x => 2))
end

Base.:+(x::Number, y::PPolynomial) = add!(copy(y), x)
Base.:+(x::PPolynomial, y::Number) = add!(copy(x), y)

Base.:+(x::PMonomial, y::PPolynomial) = add!(copy(y), x)
Base.:+(x::PPolynomial, y::PMonomial) = add!(copy(x), y)

Base.:+(x::PPolynomial, y::PPolynomial) = add!(copy(x), y)



Base.:-(x::PMonomial) = PPolynomial(Dict(x => -1))
Base.:-(x::PPolynomial) = PPolynomial((-c, m) for (c, m) in x)



Base.:-(x::Number, y::PMonomial) = add!(-y, x)
Base.:-(x::PMonomial, y::Number) = sub!(PPolynomial(x), y)

function Base.:-(x::PMonomial, y::PMonomial)
    return PPolynomial((x != y) ? Dict(x => 1, y => -1) : Dict())
end

Base.:-(x::Number, y::PPolynomial) = add!(-y, x)
Base.:-(x::PPolynomial, y::Number) = sub!(copy(x), y)

Base.:-(x::PMonomial, y::PPolynomial) = sub!(PPolynomial(x), y)
Base.:-(x::PPolynomial, y::PMonomial) = sub!(copy(x), y)

Base.:-(x::PPolynomial, y::PPolynomial) = sub!(copy(x), y)



Base.:*(x::Number, y::PMonomial) = PPolynomial(x, y)
Base.:*(x::PMonomial, y::Number) = PPolynomial(y, x)

function Base.:*(x::PMonomial, y::PMonomial)
    product = join_monomials(x, y)

    if product isa Tuple
        (c, m) = product
        return PPolynomial(c, m)
    else
        return product
    end
end

function Base.:*(x::Number, y::PPolynomial)
    return (x != 0) ? PPolynomial((x*c, m) for (c, m) in y) : 0
end

function Base.:*(x::PPolynomial, y::Number)
    return (y != 0) ? PPolynomial((y*c, m) for (c, m) in x) : 0
end

function Base.:*(x::PMonomial, y::PPolynomial)
    z = PPolynomial()

    for (c, m) in y
        addmul!(z, c, x*m)
    end

    return z
end

function Base.:*(x::PPolynomial, y::PMonomial)
    z = PPolynomial()

    for (c, m) in x
        addmul!(z, c, m*y)
    end

    return z
end

function Base.:*(x::PPolynomial, y::PPolynomial)
    z = PPolynomial()

    for (cx, mx) in x
        for (cy, my) in y
            addmul!(z, cx*cy, mx*my)
        end
    end

    return z
end



Base.:/(x::PMonomial, y::Number) = PPolynomial(rdiv(1, y), x)

function Base.:/(x::PPolynomial, y::Number)
    divs = ((rdiv(c, y), m) for (c, m) in x)
    return PPolynomial((c, m) for (c, m) in divs if c != 0)
end



function Base.:^(x::Union{PMonomial,PPolynomial}, p::Integer)
    @assert p >= 0

    result = ((p % 2 == 1) ? x : Id)
    p >>= 1

    while p > 0
        x *= x

        if p % 2 == 1
            result *= x
        end

        p >>= 1
    end

    return result
end



comm(x::Number, y::Number) = 0
comm(x::Number, y::PMonomial) = 0
comm(x::PMonomial, y::Number) = 0
comm(x, y) = x*y - y*x

acomm(x::Number, y::Number) = 2*rmul(x, y)
acomm(x::Number, y::PMonomial) = PPolynomial(2*x, y)
acomm(x::PMonomial, y::Number) = Polynomail(2*y, x)
acomm(x, y) = x*y + y*x



Base.:(==)(x::Number, y::PPolynomial) = isempty(y - x)
Base.:(==)(x::PPolynomial, y::Number) = isempty(x - y)

Base.:(==)(x::PMonomial, y::PPolynomial) = isempty(y - x)
Base.:(==)(x::PPolynomial, y::PMonomial) = isempty(x - y)

Base.:(==)(x::PPolynomial, y::PPolynomial) = isempty(x - y)



function Base.conj(x::PPolynomial)
    return PPolynomial((conj(c), conj(m)) for (c, m) in x)
end

function Base.adjoint(x::PPolynomial)
    return PPolynomial((adjoint(m), adjoint(c)) for (c, m) in x)
end

Base.zero(::PPolynomial) = PPolynomial()

Base.length(p::PPolynomial) = length(p.terms)



function conj_min(p::PPolynomial,cyclic::Bool)
    return psum(conj_min(c) * conj_min(m,true) for (c, m) in p)
end

# function cyclic_conj_min(p::PPolynomial)
#     return psum(cyclic_conj_min(c) * cyclic_conj_min(m) for (c, m) in p)
# end


function trace(m::PMonomial)
    result = ctrace(m)

    if result isa Tuple
        (c, tm) = result
        return PPolynomial(c, tm)
    else
        return result
    end
end

trace(p::PPolynomial) = psum(c*trace(m) for (c, m) in p)

function P2PP(p::Polynomial)
    pp=PPolynomial(Dict())
    for (key,value) in p.terms
        if !haskey(pp.terms,M2PM(key))
            pp.terms[M2PM(key)]=value
        else
            pp.terms[M2PM(key)]+=value
        end
    end
    return pp
end