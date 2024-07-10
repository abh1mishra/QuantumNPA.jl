# QuantumNPA

Code to do cyclic NPA in Julia. In development - names of important functions or
even the entire project could change.

Prerequisites:
```julia
using Pkg; Pkg.add(["Combinatorics", "JuMP", "SCS", "BlockDiagonals", "Mosek", "MosekTools"])
```

Then to use or try out:
```julia
include("qnpa.jl");
```

# Example from arXiv:2007.16145 (Fig. 3 red and blue lines)

Let us start computing a point from the plot with pure states (red line)
```
G=0.8  # value of the guessing probability
σ = freeop([1],1)
ρ = projector([1], 1, 1:3)
PA = projector([1], 1:2, 4:5, full=true)
ge = [σ-1/3*ρ[i] for i in 1:3]
ge_constr = [[-σ,-G]]
eq_constr = [[ρ[x],1] for x in 1:3]
wit = -(2*PA[1,1]-1)*ρ[1]-(2*PA[1,2]-1)*ρ[1]-(2*PA[1,1]-1)*ρ[2]+(2*PA[1,2]-1)*ρ[2]+(2*PA[1,1]-1)*ρ[3]
level = 1 # this is the level of the localising matrices and it will be enought to retrive the values of the plot (they use level 3 of the principal moment matrix)
-npa_general(-wit, level; op_ge=ge, tr_eq=eq_constr, tr_ge=ge_constr)
```

Now we compute a point from the blue line (mixed states)
```
G=0.8  # value of the guessing probability
σ = freeop([1],1)
ρ = freeop([1],2:4) # the input numbers need to be different from the ones of σ or they will be recognised as the same operator
PA = projector([1], 1:2, 4:5, full=true)
ge = [σ-1/3*ρ[i] for i in 1:3]
ge_add = [ρ[i]-ρ[i]^2 for i in 1:3]
ge=vcat(ge,ge_add)
ge_constr = [[-σ,-G]]
eq_constr = [[ρ[x],1] for x in 1:3]
wit = -(2*PA[1,1]-1)*ρ[1]-(2*PA[1,2]-1)*ρ[1]-(2*PA[1,1]-1)*ρ[2]+(2*PA[1,2]-1)*ρ[2]+(2*PA[1,1]-1)*ρ[3]
level = 2 # level of the localising matrices. Here we need level 2 to retrive the values of the plot (and actually improve them because of our complete set of constraints)
-npa_general(-wit, level; op_ge=ge, tr_eq=eq_constr, tr_ge=ge_constr)
```
