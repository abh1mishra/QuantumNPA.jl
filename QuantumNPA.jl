module QuantumNPA

export npa_min, npa_max, npa2sdp, sdp2jump, npa2jump,
    set_solver!,
    Id, @dichotomic,
    dichotomic, fourier, unitary, projector, zbff,
    Monomial, Polynomial, PPolynomial, monomials, coefficients,
    cglmp,flatMonomial,PMonomial,M2PM,opcycles,freeop

include("qnpa.jl")

end
