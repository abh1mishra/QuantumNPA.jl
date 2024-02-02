include("groebner_utils.jl")
include("ops_utils.jl")
using BlockDiagonals
using Base.Iterators
using SparseArrays
using Mosek, MosekTools
using JuMP
Moments = Dict{slpalg{n_Q}}{VariableRef}

function sparse_sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    # if i != j
    #     matrix[j, i] += val
    # end
end

function sparse_Nt(I, J, V, M, N)
    return SparseMatrixCSC{Float64}{Int}(sparse(I, J, V, M, N))
end

function sparse_sym(N, i, j, val)
    
    if i == j
        return sparse_Nt([i], [i], [val], N, N)
    else
        # return sparse_Nt([i, j], [j, i], [val, val], N, N)
        return sparse_Nt([i], [j], [val], N, N)


    end
end


function npa_moments_block(operators,cPoly,g::grobnerNPA,model::Model)
    N = length(operators)
    iops=sort(collect(operators))
    iops = collect(enumerate(iops))
    block = Array{AffExpr}(undef,N,N)
    mons=Moments()
    G = g.G[1]
    for (i, x) in iops
        for (j, y) in iops
            p=adj(x,g)*cPoly*y
            # p=1//2*(p+adj(p,g))
            p = reduce(p,G)
            println("$i,$j",p)
            if !haskey(mons, p)
                var_m=@variable(model)
                block[i,j] = 1*var_m
                mons[p]=var_m
            else
                block[i,j] = 1*mons[p]
            end
            # p=slpalg2dict(p,G)
            # for (m, c) in p
            #     if !haskey(block, m)
            #         block[m] = sparse_sym(N, i, j, c)
            #     else
            #         sparse_sym_add!(block[m], i, j, c)
            #     end
            # end
        end
    end

    return block,mons
end

"""
Construct the NPA moment matrix.

The argument operators can in general be an array of arrays of operators
(blocks), e.g. [[Id], [A1 + B1], [A1 - B1]]. It can also be a simple array of
operators, in which case it is treated the same as an array containing a
single array of operators, e.g., [[Id, A1, A2]]). In either case the return
value is a dictionary with:

  * as keys: monomials obtained by multiplying operators in the same blocks
    together.

  * as values: block-diagonal sparse matrices with coefficients obtained
    from multiplying the input operators together.

When the second argument is not specified, the function will build a 
principal moment matrix. When one includes a polynomial in the second 
argument, the ouput will be the localizing moment matrix of the 
polynomial specified.

"""
function npa_moments(operators, cPoly,g::grobnerNPA)
    if isempty(operators)
        return moments
    end

    if first(operators) isa Union{Number,slpalg{n_Q}}
        operators = [operators]
    end
    nblocks = length(operators)
    bsizes = length.(operators)
    blocks = [npa_moments_block(i,cPoly,g) for i in operators]
    ms = union([Set(keys(block)) for block in blocks]...)
    moments = Moments()

    for m in ms
        blocks_m = [(haskey(block, m)
                     ? block[m]
                     : (n -> spzeros(n, n))(bsizes[b]))
                    for (b, block) in enumerate(blocks)]

        moments[m] = BlockDiagonal(blocks_m)
    end

    return moments
end



function SparseArrays.dropzeros!(matrix::BlockDiagonal)
    for blk in blocks(matrix)
        dropzeros!(blk)
    end

    return matrix
end

sp1x1(x) = (iszero(x)
            ? spzeros(Float64, 1, 1)
            : SparseMatrixCSC{Float64,Int}(sparse([1], [1], x)))

"""
Generate the NPA relaxation for a given quantum optimisation problem (an
operator expr whose expectation we want to maximise with the expectation
values of the operators constraints set to zero).
"""

function bspzeros(bsizes)
    return BlockDiagonal([spzeros(n, n) for n in bsizes])
end

function Base.zero(bm::BlockDiagonal)
    return bspzeros(first.(blocksizes(bm)))
end

function BlockDiagonals.blocksizes(moments::Moments)
    if isempty(moments)
        return []
    else
        return first.(blocksizes(first(moments)[2]))
    end
end




"""
Generate the SDP relaxation of a generic polynomial optimization problem
given equality and inequality constraints. The level can be specified
and it refers to the level of the localizing matrices. The words can 
be defined commuting or non commuting.
"""
function npa_general( obj,level ,g::grobnerNPA; 
                    op_ge = 0 ,
                    av_eq = 0,
                    av_ge = 0,
                    show_moments=false,
                    verbose=false)
    pRing=g.symArr[1].parent
    Id=one(pRing)
    G = std(Ideal(g.ring,g.consArr))
    push!(g.G,G)
    # print(G)
    objP=slpalg2dict(obj,G)
    op_geP=slpalg2dict(op_ge,G)
    av_geP=slpalg2dict(av_ge,G)
    av_eqP=slpalg2dict(av_eq,G)


    opsStr=Arr2Str([obj,av_eq, av_ge, op_ge, g.consArr],pRing)
    ops_principal = ops_at_level(opsStr, level)
    
    # ops_principal = unique([ops_add[o]*ops[p] 
    #                         for o in 1:length(ops_add) for p in 1:length(ops)])
    
    model = Model(Mosek.Optimizer)

    moments_p = npa_moments(ops_principal,Id,g)
    # moments_p = npa_moments([Id,g.symArr...,g.symArr[1]*g.symArr[2],g.symArr[3]*g.symArr[1]],Id,g)
    mons_p = keys(moments_p)
    @variable(model, Γ[mons_p])
    @constraint(model,
                sum(Γ[m].*moments_p[m] for m in mons_p) >= 0,
                PSDCone())
    @constraint(model, Γ[Id]==1)

    if av_eq!=0
        av_eq=[[reduce(av_eq[x][1],G),av_eq[x][2]] for x in 1:length(av_eq)]
        [@constraint(model, sum(Γ[m]*av_eq[x][1][m] for m in mons_p) == av_eq[x][2]) for x in 1:length(av_eq) ]
    end

    if av_ge!=0
        av_ge=[[reduce(av_ge[x][1],G),av_ge[x][2]] for x in 1:length(av_ge)]
        [@constraint(model, sum(Γ[m]*av_ge[x][1][m] for m in mons_p) >= av_ge[x][2]) for x in 1:length(av_ge) ]
    end


    if op_ge!=0
        pol = Arr2Str(op_ge,pRing)
        lDeg = Int(floor((2*level-degree(pol))/2))
        if lDeg<1
        @error "Not enough level"
        end
        ops_add = ops_at_level(opsStr, lDeg)
        moments_ge = [npa_moments(ops_add,op_ge[x],g) for x in 1:length(op_ge)]
        mons_ge = [keys(moments_ge[x]) for x in 1:length(op_ge)]
    
        [@constraint(model,
                    sum(Γ[m].*moments_ge[x][m] for m in mons_ge[x]) >= 0,
                    PSDCone()) for x in 1:length(op_ge)]
    end
    @objective(model, Min, sum(objP[m]*Γ[m] for m in keys(objP)))
    if !verbose
        set_silent(model)
    end
    optimize!(model)
    return model,moments_p,ops_principal
    # println(termination_status(model))
    if show_moments==false
        return objective_value(model)
    else
        if op_ge==0
            return objective_value(model), sum(value(Γ[m])*moments_p[m] for m in mons_p)
        else
            return objective_value(model), sum(value(Γ[m])*moments_p[m] for m in mons_p), 
                [sum(value(Γ[m])*moments_ge[x][m] for m in mons_ge[x]) for x in 1:length(op_ge)]
        end
    end
end

