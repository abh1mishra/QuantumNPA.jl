include("groebner_npa.jl")  

function mixer(x::slpalg{n_Q})
    letters = g.symArr
    ops = collect(exponent_words(x))
    dual=1
    ops=ops[1]
    for i in ops
        dual*=letters[(((i+7)%16)+1)]
    end
    return Γ[reduce(dual,g.G[1])]+Γ[reduce(x,g.G[1])]
end

g=grobnerNPA(16,0,4)

Aa1,Aa2,Ba1,Ba2,Ca1,Ca2,Da1,Da2,Ab1,Ab2,Bb1,Bb2,Cb1,Cb2,Db1,Db2=g.symArr

trichotor(g,g.symArr)
Aa = [Aa1, Aa2]
Ba = [Ba1, Ba2]
Ca = [Ca1, Ca2]
Da = [Da1, Da2]
a_ops=[Aa...,Ba...,Ca...,Da...]

Ab = [Ab1, Ab2]
Bb = [Bb1, Bb2]
Cb = [Cb1, Cb2]
Db = [Db1, Db2]
b_ops=[Ab...,Bb...,Cb...,Db...]

#Apply the commutation relations
for a in [Aa;Ca]
    for b in [Ba;Da]
        commute(g, [a, b])
    end
end

for a in [Ab;Cb]
    for b in [Bb;Db]
        commute(g, [a, b])
    end
end

commute(g,[Ca[1],Ca[2]])
commute(g,[Db[1],Db[2]])


pRing=g.symArr[1].parent
Id=one(pRing)
G = std(Ideal(g.ring,g.consArr))
println("Ideal generated")
push!(g.G,G)
mons_a=ops_at_level(Set([a_ops...]), 4)
println("Monomials_A created",length(mons_a))
mons_b=ops_at_level(Set([b_ops...]), 4)
model = Model(Mosek.Optimizer)
Γa,mons_a_var=npa_moments_block(mons_a,Id,g,model)
Γb,mons_b_var=npa_moments_block(mons_b,Id,g,model)
println("Principal matrix created")
# mons_p_a = keys(Γa)
# mons_p_b = keys(Γb)

# @variable(model, Γ_a[mons_p_a])
# @variable(model, Γ_b[mons_p_b])
Γ=Dict(union(mons_a_var,mons_b_var))
# Γ=Dict{slpalg{n_Q},VariableRef}(union([i => Γ_a[i] for i in mons_p_a if i!=Id],[i => Γ_b[i] for i in mons_p_b if i!=Id]))
println("vars created")
# function get_sum(Γ_,Γ, mons_p_)
#     sum_ = 0*Γ_[first(mons_p_)].*Γ[first(mons_p_)]
#     for m in mons_p_
#         sum_ += Γ_[m].*Γ[m]
#     end
#     return sum_
# end
# sum_a = get_sum(Γ_a,Γa,mons_p_a)
println("sum created")
# println("type of matrix",typeof(sum(Γ_a[m].*Γa[m] for m in mons_p_a)))
@constraint(model,
            Γa >= 0,
            PSDCone())
println("PSD constraints created for A")
@constraint(model,
            Γb >= 0,
            PSDCone())
println("PSD constraints created")

η_AS = 1
η_BS = 1
η_AL = 2/3
@variable(model, η_BL)
@constraint(model, η_BL<=1)
visibility = 1

function get_mixedmoment_eq(η_AS, η_AL, η_BS, η_BL, visibility)
    v = visibility
    for x in 1:2
        @constraint(model, mixer(Aa[x]) == 0)
        @constraint(model, mixer(Aa[x]^2) == η_AS)
        @constraint(model, mixer(Ba[x])== 0)
        @constraint(model, mixer(Ba[x]^2) == η_BS)
        @constraint(model, mixer(Ca[x])== 0)
        @constraint(model, mixer(Ca[x]^2)== η_AL)
        @constraint(model, mixer(Da[x])== 0)
        @constraint(model, mixer(Da[x]^2) == η_BL)
    end
    for x in 1:2
        for y in 1:2
            @constraint(model, mixer(Aa[x]*Ba[y]) - (v^2 * (-1)^((x-1)*(y-1)) * η_AS * η_BS * 1/√2) == 0)
            @constraint(model, mixer(Aa[x]*Da[y]) - (v^2 *(-1)^((x-1)*(y-1)) * η_AS * η_BL * 1/√2) == 0)
            @constraint(model, mixer(Ca[x]*Ba[y])- (v^2 *(-1)^((x-1)*(y-1)) * η_AL * η_BS * 1/√2) == 0)
            @constraint(model, mixer(Ca[x]*Da[y])- (v^2 *(-1)^((x-1)*(y-1)) * η_AL * η_BL * 1/√2) == 0)

            @constraint(model, mixer(Aa[x]^2*Ba[y]) == 0)
            @constraint(model, mixer(Aa[x]^2*Da[y]) == 0)
            @constraint(model, mixer(Ca[x]^2*Ba[y])== 0)
            @constraint(model, mixer(Ca[x]^2*Da[y]) == 0)

            @constraint(model, mixer(Aa[x]*Ba[y]^2)== 0)
            @constraint(model, mixer(Aa[x]*Da[y]^2)== 0)
            @constraint(model, mixer(Ca[x]*Ba[y]^2)== 0)
            @constraint(model, mixer(Ca[x]*Da[y]^2)== 0)

            @constraint(model, mixer(Aa[x]^2*Ba[y]^2) - η_AS * η_BS == 0)
            @constraint(model, mixer(Aa[x]^2*Da[y]^2) - η_AS * η_BL == 0)
            @constraint(model, mixer(Ca[x]^2*Ba[y]^2) - η_AL * η_BS == 0)
            @constraint(model, mixer(Ca[x]^2*Da[y]^2) - η_AL * η_BL == 0)

        end
    end
end

get_mixedmoment_eq(η_AS, η_AL, η_BS, η_BL, visibility)
println("addl. constraints created")
@objective(model, Max, η_BL)
@constraint(model, mons_a_var[Id]+mons_b_var[Id]==1)
# #@objective(model, Max, mixer(Aa[1]*Ba[1]) + mixer(Aa[1]*Ba[2]) + mixer(Aa[2]*Ba[1]) - mixer(Aa[2]*Ba[2]))
# #@objective(model, Max, mixer(Ca[1]*Da[1]) + mixer(Ca[1]*Da[2]) + mixer(Ca[2]*Da[1]) - mixer(Ca[2]*Da[2]))
optimize!(model)
println(termination_status(model),objective_value(model))
