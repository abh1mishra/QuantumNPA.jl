PMoments = Dict{PMonomial}{BlockDiagonal}
Moments = Dict{Monomial}{BlockDiagonal}

function sparse_sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end
end

function sparse_Nt(I, J, V, M, N)
    return SparseMatrixCSC{Float64}{Int}(sparse(I, J, V, M, N))
end

function sparse_sym(N, i, j, val)
    if i == j
        return sparse_Nt([i], [i], [val], N, N)
    else
        return sparse_Nt([i, j], [j, i], [val, val], N, N)
    end
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


function npa_moments_block(operators,cPoly)
    N = length(operators)
    iops = collect(enumerate(operators))
    block = Dict{Monomial,SparseMatrixCSC}()

    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(conj(x,false)*cPoly*y)
            for (c, m) in p
                if !haskey(block, m)
                    block[m] = sparse_sym(N, i, j, c)
                else
                    sparse_sym_add!(block[m], i, j, c)
                end
            end
        end
    end

    return block
end

function cyclic_npa_moments_block(operators,cPoly)
    N = length(operators)
    iops = collect(enumerate(operators))
    block = Dict{PMonomial,SparseMatrixCSC}()
    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(conj(x,false)*cPoly*y)
            for (c, m) in p
                pm=M2PM(m)
                pmk=haskey(block, pm)
                if !pmk
                    block[pm] = sparse_sym(N, i, j, c)
                else
                    sparse_sym_add!(block[pm], i, j, c)
                end
            end
        end
    end

    return block
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
function npa_moments(operators,cPoly=Id)
    if isempty(operators)
        return moments
    end

    if first(operators) isa Union{Number,PMonomial,Polynomial,Monomial}
        operators = [operators]
    end

    nblocks = length(operators)
    bsizes = length.(operators)
    blocks = [npa_moments_block(i,cPoly) for i in operators]
    ms = monomials(keys(block) for block in blocks)

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

"""
Construct the NPA moment matrix for tracial hierarchy.

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
function cyclic_npa_moments(operators, cPoly=Id)
    if isempty(operators)
        return moments
    end

    if first(operators) isa Union{Number,Monomial,Polynomial}
        operators = [operators]
    end

    nblocks = length(operators)
    bsizes = length.(operators)
    blocks = [cyclic_npa_moments_block(i,cPoly) for i in operators]
    ms = monomials(keys(block) for block in blocks)

    Pmoments = PMoments()

    for m in ms
        blocks_m = [(haskey(block, m)
                     ? block[m]
                     : (n -> spzeros(n, n))(bsizes[b]))
                    for (b, block) in enumerate(blocks)]

        Pmoments[m] = BlockDiagonal(blocks_m)
    end

    return Pmoments
end


# Creates moment matrix starting from a polynomial Op
function cyclic_pol_moments(Op::Polynomial, operators)
    s=size(operators)[1]
    res=Dict(P2PP(operators[i]*Op*operators[j]) => zeros(s,s)
                                    for i in 1:s for j in 1:s)
    for i in 1:s
        for j in 1:s
            res[P2PP(operators[i]*Op*operators[j])][i,j]=1
        end
    end
    return res
end




"""
Generate the SDP relaxation of a generic polynomial optimization problem
given equality and inequality constraints. The level can be specified
and it refers to the level of the localizing matrices. The words can 
be defined commuting or non commuting.
"""

function get_monomials(obj, level; 
                        op_eq = 0, 
                        op_ge = 0,
                        tr_eq = 0,
                        tr_ge = 0)
    # for a given level of the localising moment matrices, it returns the operators at that level and the ones of the principal moment matix (which will be larger)
    ops = ops_at_level(Polynomial.(vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])), level)
    pol = 1+sum(op_ge)+sum(op_eq)
    deg = Int(ceil(degree(pol)/2))
    ops_add = ops_at_level([op_ge,op_eq], deg)
    ops_principal = unique([ops_add[o]*ops[p] 
                            for o in 1:length(ops_add) for p in 1:length(ops)])
    return ops, ops_principal
end


"""
    npa_general(obj, level; op_eq=0, op_ge=0, tr_eq=0, tr_ge=0, show_moments=false, verbose=false, termination=false, level_principal=0)

This function performs a Noncommutative Polynomial Algebra (NPA) hierarchy optimization using the Mosek optimizer. It constructs and solves a semidefinite programming (SDP) problem based on the provided polynomial constraints and objectives.

# Arguments
- `obj`: The objective polynomial to be minimized.
- `level`: The level of the NPA hierarchy.
- `op_eq`: A list of operator equality constraints. Default is 0.
- `op_ge`: A list of operator inequality constraints (greater than or equal to). Default is 0.
- `tr_eq`: A list of trace equality constraints. Default is 0.
- `tr_ge`: A list of trace inequality constraints (greater than or equal to). Default is 0.
- `show_moments`: A boolean flag indicating whether to return the moments. Default is false.
- `verbose`: A boolean flag indicating whether to print detailed solver output. Default is false.
- `termination`: A boolean flag indicating whether to return the termination status. Default is false.
- `level_principal`: The level of the principal operators in the NPA hierarchy. Default is 0.

# Returns
- If `show_moments` is false:
    - If `termination` is true: Returns a tuple `(objective_value, termination_status)`.
    - Otherwise: Returns `objective_value`.
- If `show_moments` is true:
    - If `op_ge` is 0:
        - If `termination` is true: Returns a tuple `(objective_value, moments_sum, termination_status)`.
        - Otherwise: Returns a tuple `(objective_value, moments_sum)`.
    - If `op_ge` is not 0:
        - If `termination` is true: Returns a tuple `(objective_value, moments_sum, moments_ge_sum, termination_status)`.
        - Otherwise: Returns a tuple `(objective_value, moments_sum, moments_ge_sum)`.
# Details
The function first determines the monomials and principal operators based on the level_principal parameter.
It then creates a model using the Mosek optimizer.
Moments and monomials for the principal operators are generated.
The function defines the variable Γ in the model.
A PSD constraint is added for the principal moments.
Constraints for tr_eq and tr_ge are added if they are not 0.
Constraints for op_eq and op_ge are added if they are not 0.
The objective function is defined.
The model is set to silent if verbose is false.
The model is optimized.
The termination status is printed.
The function returns the objective value and optionally the moments and termination status based on the flags show_moments and termination.
# Example
```julia
obj = Polynomial(...)
level = 2
op_eq = [...]
op_ge = [...]
tr_eq = [...]
tr_ge = [...]
result = npa_general(obj, level; op_eq=op_eq, op_ge=op_ge, tr_eq=tr_eq, tr_ge=tr_ge, show_moments=true, verbose=true, termination=true)
"""
function npa_general( obj, level;
                    op_eq = 0, 
                    op_ge = 0,
                    tr_eq = 0,
                    tr_ge = 0,
                    show_moments = false,
                    verbose = false,
                    termination = false,
                    level_principal=0)
                        
    if level_principal==0
            ops, ops_principal = get_monomials(obj,level; op_eq = op_eq, op_ge = op_ge, tr_eq = tr_eq, tr_ge = tr_ge)
        else
            ops_principal = ops_at_level(Polynomial.(vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])), level)
            ops = ops_at_level(Polynomial.(vcat([obj, op_ge..., op_eq...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])), level_principal)
    end

    model = Model(Mosek.Optimizer)
    moments_p = cyclic_npa_moments(ops_principal)
    mons_p = keys(moments_p)
    # return ops_principal,mons_p
    @variable(model, Γ[mons_p])
    @constraint(model,
                sum(Γ[m].*moments_p[m] for m in mons_p) >= 0,
                PSDCone())

    if tr_eq!=0
        tr_eq=[[P2PP(1*tr_eq[x][1]),tr_eq[x][2]] for x in 1:length(tr_eq)]
        [@constraint(model, sum(Γ[m]*tr_eq[x][1][m] for m in mons_p) == tr_eq[x][2]) for x in 1:length(tr_eq) ]
    end

    if tr_ge!=0
        tr_ge=[[P2PP(1*tr_ge[x][1]),tr_ge[x][2]] for x in 1:length(tr_ge)]
        [@constraint(model, sum(Γ[m]*tr_ge[x][1][m] for m in mons_p) >= tr_ge[x][2]) for x in 1:length(tr_ge) ]
    end

    if op_eq!=0
        moments_eq = [cyclic_npa_moments(ops,op_eq[x]) for x in 1:length(op_eq)]
        mons_eq = [keys(moments_eq[x]) for x in 1:length(op_eq)]

        [@constraint(model,
                sum(Γ[m].*moments_eq[x][m] for m in mons_eq[x]) >= 0,
                PSDCone()) for x in 1:length(op_eq)]

        [@constraint(model,
                sum(Γ[m].*moments_eq[x][m] for m in mons_eq[x]) <= 0,
                PSDCone()) for x in 1:length(op_eq)]
    end
    if op_ge!=0
        moments_ge = [cyclic_npa_moments(ops,op_ge[x]) for x in 1:length(op_ge)]
        mons_ge = [keys(moments_ge[x]) for x in 1:length(op_ge)]
        [@constraint(model,
                    sum(Γ[m].*moments_ge[x][m] for m in mons_ge[x]) >= 0,
                    PSDCone()) for x in 1:length(op_ge)]
    end
    obj=P2PP(obj)
    @objective(model, Min, sum(c*Γ[m] for (c,m) in obj))
    if !verbose
        set_silent(model)
    end
    optimize!(model)
    println(termination_status(model),objective_value(model),"\n")
    if show_moments==false
        if termination==true
            return objective_value(model), termination_status(model)
        else
            return objective_value(model)
        end
    else
        if op_ge==0
            if termination==true
                return objective_value(model), sum(value(Γ[m])*moments_p[m] for m in mons_p), termination_status(model)
            else
                return objective_value(model), sum(value(Γ[m])*moments_p[m] for m in mons_p)
            end
        else
            if termination==true
                return objective_value(model), sum(value(Γ[m])*moments_p[m] for m in mons_p), 
                    [sum(value(Γ[m])*moments_ge[x][m] for m in mons_ge[x]) for x in 1:length(op_ge)], termination_status(model)
            else
                return objective_value(model), sum(value(Γ[m])*moments_p[m] for m in mons_p), 
                    [sum(value(Γ[m])*moments_ge[x][m] for m in mons_ge[x]) for x in 1:length(op_ge)]
            end
        end
    end
end

function npa_model(level; obj=0,

    op_eq = 0, 
 
    op_ge = 0,
 
    tr_eq = 0,
 
    tr_ge = 0,
 
    add_ops=0)
 
 
    println("obj:",obj,"\nop_eq:",op_eq,"\nop_ge:",op_ge,"\ntr_eq:",tr_eq,"\ntr_ge:",tr_ge)
    print(vcat([obj, op_ge..., op_eq..., add_ops...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)]))
    ops = ops_at_level(Polynomial.(vcat([obj, op_ge..., op_eq..., add_ops...],[tr_ge[i][1] for i in 1:length(tr_ge)], [tr_eq[i][1] for i in 1:length(tr_eq)])), level)
 
    pol = 1+sum(op_ge)+sum(op_eq)
 
    deg = Int(ceil(degree(pol)/2))
 
    ops_add = ops_at_level([op_ge,op_eq], deg)
 
    ops_principal = unique([ops_add[o]*ops[p] 
 
                for o in 1:length(ops_add) for p in 1:length(ops)])
    
    println(length(ops_principal))
 
 
 
 
    model = Model(Mosek.Optimizer)
 
 
 
 
    moments_p = cyclic_npa_moments(ops_principal)
 
    mons_p = keys(moments_p)
    @variable(model, Γ[mons_p])
 
    @constraint(model,
 
    sum(Γ[m].*moments_p[m] for m in mons_p) >= 0,
 
    PSDCone())
    
 
 
    if tr_eq!=0
 
        tr_eq=[[P2PP(1*tr_eq[x][1]),tr_eq[x][2]] for x in 1:length(tr_eq)]
    
        [@constraint(model, sum(Γ[m]*tr_eq[x][1][m] for m in mons_p) == tr_eq[x][2]) for x in 1:length(tr_eq) ]
    
    end
 
 
 
 
    if tr_ge!=0
 
        tr_ge=[[P2PP(1*tr_ge[x][1]),tr_ge[x][2]] for x in 1:length(tr_ge)]
    
        [@constraint(model, sum(Γ[m]*tr_ge[x][1][m] for m in mons_p) >= tr_ge[x][2]) for x in 1:length(tr_ge) ]
 
    end
 
 
 
 
    if op_eq!=0
 
        moments_eq = [cyclic_npa_moments(ops,op_eq[x]) for x in 1:length(op_eq)]
    
        mons_eq = [keys(moments_eq[x]) for x in 1:length(op_eq)]
    
    
    
    
        [@constraint(model,
    
        sum(Γ[m].*moments_eq[x][m] for m in mons_eq[x]) .== 0,
    
        PSDCone()) for x in 1:length(op_eq)]
    
    
    end
 
    if op_ge!=0
 
        moments_ge = [cyclic_npa_moments(ops,op_ge[x]) for x in 1:length(op_ge)]
    
        mons_ge = [keys(moments_ge[x]) for x in 1:length(op_ge)]
    
    
    
    
        [@constraint(model,
    
            sum(Γ[m].*moments_ge[x][m] for m in mons_ge[x]) >= 0,
    
            PSDCone()) for x in 1:length(op_ge)]
 
    end
 
    if obj!=0
 
        obj=P2PP(obj)
 
        @objective(model, Min, sum(obj[m]*Γ[m] for m in mons_p))
 
    end
    
    return model, Γ, mons_p, moments_p
 
 end
function npa2sdp(pMatrix,lMatrixArray)

end