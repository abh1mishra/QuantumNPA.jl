include("qnpa.jl")
G=0.8  # value of the guessing probability
σ = hermitian([1],1)
ρ = projector([1], 1, 1:3)
PA = projector([1], 1:2, 4:5, full=true)
ge = [σ-1/3*ρ[i] for i in 1:3]
ge_constr = [[-σ,-G]]
eq_constr = [[ρ[x],1] for x in 1:3]
wit = -(2*PA[1,1]-1)*ρ[1]-(2*PA[1,2]-1)*ρ[1]-(2*PA[1,1]-1)*ρ[2]+(2*PA[1,2]-1)*ρ[2]+(2*PA[1,1]-1)*ρ[3]
level = 1 # this is the level of the localising matrices and it will be enought to retrive the values of the plot (they use level 3 of the principal moment matrix)
-npa_general(-wit, level; op_ge=ge, tr_eq=eq_constr, tr_ge=ge_constr)
