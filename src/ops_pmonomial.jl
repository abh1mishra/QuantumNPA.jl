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
Base.hash(pm::PMonomial, h::UInt) = hash(length(pm.pword), h)

function Base.show(io::IO, x::PMonomial)
    if isidentity(x)
        print(io, "Id")
    else
        for k in sort(collect(keys(x.pword)))
            print(io,"_(",Monomial(x.pword[k]),")_")
        end
    end
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
        
        while true
            if length(value)<2
                break
            end
            pfirst=value[1][1]
            plast=value[end][1]
            if pfirst!=plast
                break
            end
            if sort(pfirst)[1]< key
                break
            end
            p_rest=subsElems(pfirst,key)
            ofirst=value[1][2][1]
            olast=value[end][2][end]
            if !(all([ofirst==a.pword[i][1][2][1] for i in p_rest]) && all([olast==a.pword[i][end][1] for i in p_rest]))
                break
            end
            prod=ofirst*olast
            if(prod[1]==0)
                return 0
            end
            if(isempty(prod[2]))
                pop!(value)
                popfirst!(value)
                for i in p_rest
                    pop!(a.pword[i])
                    popfirst!(a.pword[i])
                end
                continue
            end
            if(length(prod[2])==1)
                value[end]=(pfirst,prod[2])
                popfirst!(value)
                for i in p_rest
                    a.pword[i][end]=(a.pword[i][1][1],prod[2])
                    popfirst!(a.pword[i])
                end
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
                # conjStr=string(conj(Monomial(y.pword[key]),false))
                # if !( occursin(string(Monomial(value)), Str*" "*Str ) || occursin(string(Monomial(value)), conjStr*" "*conjStr ) )
                if !occursin(string(Monomial(value)), Str*" "*Str )
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