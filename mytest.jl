include("qnpa.jl")
using LinearAlgebra
using DelimitedFiles
using Plots



# function key(x)
#   y= [1-log2(1+sqrt(2-x[i]*x[i]/4)) for i in 1:length(x)]      
#   return y
# end

N=2             

alpha=0
p=0.5*ones(2)
lvl="1+B E" # + B E + B B,E + B E B,E
#lvl="1+B E + B,E"

PB = projector([2], 1:2, 1:2, full=true)
PE = projector([5], 1:2, 1, full=true)
# if proj==true
#K = projector([2,5],1,1:4)
# else
K = hermitian([2,5],1:4)
ρ=[K[1] K[2];K[3] K[4]]
  #end
σ = hermitian([2,5],5)


G=(1+(N-1)*alpha)/N

B1=PB[1,1]-PB[2,1]
B2=PB[1,2]-PB[2,2]
R1=p[1]*ρ[1,1]-(1-p[1])*ρ[2,1]
R2=p[2]*ρ[1,2]-(1-p[2])*ρ[2,2]
S=(B1 + B2)*R1 + (B1 - B2)*R2


# if proj==true
#    op_ge = vcat([σ-(p[x]*ρ[1,x]+(1-p[x])*ρ[2,x])/N for x in 1:2],[G*Id-σ])
# else
#
op_ge = vcat([σ-(p[x]*ρ[1,x]+(1-p[x])*ρ[2,x])/N for x in 1:2],[G*Id-σ],[ ρ[a,x]- ρ[a,x]*ρ[a,x] for a in 1:2 for x in 1:2])
# end
#op_eq = vcat([  ρ[a,x]- ρ[a,x]*ρ[a,x] for a in 1:2 for x in 1:2], [p[1]*ρ[1,1]+(1-p[1])*ρ[2,1]-(p[2]*ρ[1,2]+(1-p[2])*ρ[2,2])])
#op_eq =vcat( [p[1]*ρ[1,1]+(1-p[1])*ρ[2,1]-(p[2]*ρ[1,2]+(1-p[2])*ρ[2,2])] )
#op_ge=vcat([  ρ[a,x]- ρ[a,x]*ρ[a,x] for a in 1:2 for x in 1:2])
tr_ge = [[-σ, -G]]

obj = (p[1]*ρ[1,1]*PE[1]+(1-p[1])*ρ[2,1]*PE[2])

# --------------------------------------------------------------------------
t1=time()
s_i=2.65                #Initial CHSH value
s_f=2*sqrt(2)-0.001     #Final CHSH value
tot=1                   #total no. of points
#s=zeros(tot)            #CHSH value
pg=zeros(tot)           #Guessing probability..later converted to min Entropy..Numerical key
ky=zeros(tot)          #Analytical key
num=1
s=2.6
#for num in 1:tot
  #s[num]=s_i*(tot-num)/(tot-1)+(num-1)*s_f/(tot-1)#(2*sqrt(2)*(1-alpha)+4*alpha)
  eq_constr = vcat( [[ ρ[a,x], 1] for a in 1:2 for x in 1:2], [[ S, s]])
  ob_j = npa_general(-obj, lvl; op_ge,tr_eq=eq_constr,tr_ge, show_moments=false)
  pg[num]=-ob_j 
  pg[num]=-log2(pg[num])
  ky[num]=1-log2(1+sqrt(2-s*s/4))
  println("Level: ", lvl, ", CHSH value: ", s,", Numerical key: ",pg[num],", Analytic key: ",ky[num])
#end

println("Elapsed time: ", time()-t1, " seconds")

#plot(s,[pg,key])

#pg[num]=-log2(pg[num])

#    return pg


#---------Datafile---------------
# fw = open("lvl2_a=0.txt", "a")#keyrate_alpha=0_p=0.5_lvl_1+B E
# writedlm(fw, [s pg ky]) 
# close(fw)
