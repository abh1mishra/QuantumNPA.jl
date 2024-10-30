include("qnpa.jl")
using Plots
points = collect(range(0.3, 1, step=0.05))
σ = hermitian([1],1)
ρ = hermitian([1],2:4) # the input numbers need to be different from the ones of σ or they will be recognised as the same operator
PA = projector([1], 1:2, 4:5, full=true)
ge = [σ-1/3*ρ[i] for i in 1:3]
ge_add = [ρ[i]-ρ[i]^2 for i in 1:3]
ge=vcat(ge,ge_add)
eq_constr = [[ρ[x],1] for x in 1:3]
wit = -(2*PA[1,1]-1)*ρ[1]-(2*PA[1,2]-1)*ρ[1]-(2*PA[1,1]-1)*ρ[2]+(2*PA[1,2]-1)*ρ[2]+(2*PA[1,1]-1)*ρ[3]
level = 2 # level of the localising matrices. Here we need level 2 to retrive the values of the plot (and actually improve them because of our complete set of constraints)
res=zeros(length(points))
for (i,G) in enumerate(points)
    ge_constr = [[-σ,-G]]
    res[i]=-npa_general(-wit, level; op_ge=ge, tr_eq=eq_constr, tr_ge=ge_constr)
end
plot(points,res)