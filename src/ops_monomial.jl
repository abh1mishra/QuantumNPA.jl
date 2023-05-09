struct Monomial
    word::Array{Tuple{Array{Int64,1},Array{Operator,1}},1}
end
struct PMonomial
    pword::Dict{Int64,Array{Tuple{Array{Int64,1},Array{Operator,1}},1} }
end
function Monomial(party::Array{Int64,1}, operator::Operator)
    @assert all(>(0),party)
    return Monomial([(party, [operator])])
end

Monomial(party, operator::Operator) = Monomial(party_num(party), operator)

Id = Monomial([])

isidentity(m::Monomial) = isempty(m)
isidentity(m::PMonomial) = isempty(m.pword)
Base.iterate(m::Monomial) = iterate(m.word)
Base.iterate(m::Monomial, state) = iterate(m.word, state)

Base.length(m::Monomial) = length(m.word)

Base.hash(m::Monomial, h::UInt) = hash(m.word, h)
Base.hash(pm::PMonomial, h::UInt) = hash(length(pm.pword), h)


function Base.show(io::IO, m::Monomial)
    if isidentity(m)
        print(io, "Id")
    else
        sep = ""

        for (party, ops) in m
            for o in ops
                print(io, sep)
                print(io, string(o, party))
                sep = " "
            end
        end
    end
end

function Base.show(io::IO, x::PMonomial)
    if isidentity(x)
        print(io, "Id")
    else
        for (key,value) in x.pword
            print(io,key,"-->")
            for i in value
                print(io,Monomial(i[1],i[2]...)," ")
            end
        end
    end
end



degree(x::Number) = !iszero(x) ? 0 : -Inf

function degree(m::Monomial)
    result = 0

    for (_, ops) in m.word
        result += length(ops)
    end

    return result
end



Base.:(==)(x::Number, y::Monomial) = (x == 1) && isempty(y)

Base.:(==)(x::Monomial, y::Number) = (y == 1) && isempty(x)

Base.:(==)(x::Monomial, y::Monomial) = (x.word == y.word)

function Base.isless(x::Monomial, y::Monomial)
    ox, oy = degree(x), degree(y)

    if ox != oy
        return ox < oy
    end

    for ((p1, ops1), (p2, ops2)) in zip(x, y)
        if p1 != p2
            return p1 < p2
        end

        l1, l2 = length(ops1), length(ops2)

        if l1 != l2
            return l1 > l2
        end

        for (o1, o2) in zip(ops1, ops2)
            if o1 != o2
                return o1 < o2
            end
        end
    end

    return false
end



# function Base.conj(m::Monomial)
#     return
# end

# function Base.conj(m::Monomial,cyclic::Bool)
#     if cyclic
#         # can simplify conj for subsystems by doing
#         return reorderMonomial(Monomial([(party, reverse!([conj(op) for op in ops]))
#                      for (party, ops) in reverse(m.word)]))
#     else
#         return Monomial([(party, reverse!([conj(op) for op in ops]))
#                          for (party, ops) in m])
#     end
# end

function Base.conj(m::Monomial,cyclic::Bool)
    if cyclic
        # can simplify conj for subsystems by doing
        return reorderMonomial(Monomial([(party, reverse!([conj(op) for op in ops]))
                     for (party, ops) in reverse(m.word)]))
    else
        return Monomial( reverse!([(party, reverse!([conj(op) for op in ops]))
                         for (party, ops) in m]))
    end
end

function reorderMonomial(m::Monomial)
    monArr=[Monomial([o]) for o in m.word]
    if length(m)<=1
        return m
    end
    return *(monArr...)

end

function Base.adjoint(m::Monomial)
    return Monomial([(party, reverse!([adjoint(op) for op in ops]))
                     for (party, ops) in m])
end

Base.zero(m::Monomial) = 0



conj_min(x::Number) = real(x)

function conj_min(m::Monomial,cyclic::Bool)
    if m==Id
        return m
    end
    if !cyclic
        return min(m, conj(m))
    else
        monCycles=opcycles(flatMonomial(m).word,true)
        #=


        Need to put conjugate of operators in conjMonCycles, now not necessary as dealing with projectors


        =#
        conjMonCycles=opcycles(reverse!(flatMonomial(m).word),true)
        # println(monCycles)
        # println(conjMonCycles)
        return ((monCycles==0) | (conjMonCycles==0)) ? 0 : min(vcat(monCycles,conjMonCycles)...)
    end
end

# function cyclic_conj_min(m::Monomial)
#
# end



"""
Concatenate two monomials. This is used later to decide what the result
of multiplying two monomials is.
"""

function join_monomials(x::Monomial, y::Monomial)
    coeff = 1

    if (M = length(x)) == 0
        return y
    end

    if (N = length(y)) == 0
        return x
    end
    word=deepcopy(x.word)
    for (py,opsy) in y.word
        for j in length(word):-1:1
            (px,opsx)=word[j]
            if (intersect(py,px)==Int64[]) && ( sort(py) > sort(px) )
                insert!(word,j+1,(py,opsy))
                break
            end

            if intersect(py,px)!=Int64[]
                if py!=px
                   insert!(word,j+1,(py,opsy))
                   break
                end
                if py==px
                    (c, ops) = join_ops(opsx, opsy)
                    if c == 0
                       return 0
                    end
                    coeff *= c
                    if !isempty(ops)
                       deleteat!(word,j)
                       insert!(word,j, (px, ops))
                       break
                    else
                        deleteat!(word,j)
                        break
                    end
                end

            end
            if j==1
               insert!(word,1,(py,opsy))
            end

        end
    end
    if isempty(word)
        return Id
    end
    m = Monomial(word)
    return (coeff == 1) ? m : (coeff, m)

end



#             # if px < py
#             #     push!(word, x.word[j])
#             #     j += 1
#             # elseif py < px
#             #     push!(word, y.word[k])
#             #     k += 1
#             else
#                 (c, ops) = join_ops(opsx, opsy)
#
#                 if c == 0
#                     return 0
#                 end
#
#                 coeff *= c
#
#                 if !isempty(ops)
#                     push!(word, (px, ops))
#                 end
#
#                 j += 1
#                 k += 1
#             end
#
#     end
#
#     append!(word, x.word[j:end])
#     append!(word, y.word[k:end])
#
#     m = Monomial(word)
#
#     return (coeff == 1) ? m : (coeff, m)
# end



function ctrace(m::Monomial)
    coeff = 1

    pcops = [(p, trace(ops)) for (p, ops) in m.word]

    for (_, (c, _)) in pcops
        coeff = rmul(coeff, c)
    end

    m = Monomial([(p, ops) for (p, (_, ops)) in pcops])

    return (coeff == 1) ? m : (coeff, m)
end

function flatMonomial(m::Monomial)
    return Monomial([(s,[ops]) for (s,opsArr) in m for ops in opsArr ])
end

M2PM(x::Int64)=x
function M2PM(m::Monomial)
    a=PMonomial(Dict())
    m=flatMonomial(m)
    for i in m.word
           for j in i[1]
               if haskey(a.pword,j)
                   push!(a.pword[j],i)
               else
                   a.pword[j]=[i]
               end
           end
     end

     for (key,value) in a.pword
        monArr=[Monomial([i]) for i in value]
        push!(monArr,Id)
         a.pword[key]=flatMonomial(*(monArr...)).word
         if (Monomial([a.pword[key][1]]) == Monomial([a.pword[key][length(a.pword[key])]]) && length(a.pword[key])>1) & (typeof(a.pword[key][1][2][1])==Projector)
             pop!(a.pword[key])
         end
     end
     return a
end
 


function Base.:(==)(x::PMonomial, y::PMonomial)
    if keys(x.pword)==keys(y.pword)
        for (key,value) in x.pword
            if length(value)==length(y.pword[key])
                Str=string(Monomial(y.pword[key]))
                conjStr=string(conj(Monomial(y.pword[key]),false))
                if !( occursin(string(Monomial(value)), Str*" "*Str ) || occursin(string(Monomial(value)), conjStr*" "*conjStr ) )
                    return false
                end
            else
                return false
            end
        end
        return true
    else
        return false
    end
end
