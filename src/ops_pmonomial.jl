struct PMonomial
    pword::Dict{Int64,Array{Tuple{Array{Int64,1},Array{Operator,1}},1} }
end

function PMonomial(party::Array{Int64,1}, operator::Operator)
    @assert all(>(0),party)
    return PMonomial(Dict(Pair(p,[(party,[operator])]) for p in party))
end

PMonomial(party, operator::Operator) = PMonomial(party_num(party), operator)

PId = PMonomial(Dict{Int64,Array{Tuple{Array{Int64,1},Array{Operator,1}},1}}())

isidentity(m::PMonomial) = isempty(m.pword)


Base.iterate(m::PMonomial) = iterate(m.word)
Base.iterate(m::PMonomial, state) = iterate(m.word, state)

Base.length(m::PMonomial) = length(m.word)

Base.hash(pm::PMonomial, h::UInt) = hash(pm.pword, h)


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
Base.conj(x::Int64,cyclic::Bool)    = x

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









function flatMonomial(m::Monomial)
    return Monomial([(s,[ops]) for (s,opsArr) in m for ops in opsArr ])
end


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
        
        while true
            if length(value)<2
                break
            end
            pfirst=value[1][1]
            plast=value[end][1]
            if pfirst!=plast
                break
            end
            ofirst=value[1][2][1]
            olast=value[end][2][end]
            prod=ofirst*olast
            if(prod[1]==0)
                return 0
            end
            if(isempty(prod[2]))
                pop!(value)
                popfirst!(value)
                continue
            end
            if(length(prod[2])==1)
                value[end]=(pfirst,prod[2])
                popfirst!(value)
                continue
            end
            break
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
