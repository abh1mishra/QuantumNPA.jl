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

# Example from arXiv:2007.16145 (Fig. 3 blue line)

```
G=0.6
σ = freeop([1],1)
ρ = projector([1], 1, 1:3)
PA = projector([1], 1:2, 4:5)
ge = [σ-1/3*ρ[i] for i in 1:3]
eq = [PB[1,y]+PB[2,y]-Id for y in 1:2]
ge_constr = [[-σ,-G]]
eq_constr = [[ρ[x],1] for x in 1:3]
wit = -(2*PB[1,1]-1)*ρ[1]-(2*PB[1,2]-1)*ρ[1]-(2*PB[1,1]-1)*ρ[2]+(2*PB[1,2]-1)*ρ[2]+(2*PB[1,1]-1)*ρ[3]
level = 1
npa_general(-wit, level; op_eq=eq, op_ge=ge, tr_eq=eq_constr, tr_ge=ge_constr)
```
