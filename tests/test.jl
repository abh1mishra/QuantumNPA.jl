include("qnpa.jl")
PB=projector([2],1:2,1:2;full=true)
# PE=projector([5],1:2,1;full=true) # or uncomment next one if it doesn't work
PE=projector([5],1:2,1;full=true)

ops=ops_at_level([Id,PB[1,1:2],PE],2) # or just PE is you uncommented second command
moments= npa_moments(ops)
