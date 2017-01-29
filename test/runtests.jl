using PolynomialMatrices
using Polynomials
using Base.Test

# Include tests here
# include ("test1.jl")

# Checking ltriang
s = Poly([0, 1],:s)
P = PolyMatrix([s-1 2 0; s^2-1 2*s+2 3])
P2, Q = ltriang(P)
@test norm(coeffs(Q*P-P2)[0]) < 1e-15
@test norm(coeffs(Q*P-P2)[1]) < 1e-15
@test norm(coeffs(Q*P-P2)[2]) < 1e-15

# A low rank polynomial matrix
P = PolyMatrix([s^3 s^2 s; s^2 s 1; 0 1 2])
P2, Q = ltriang(P)
@test norm(coeffs(Q*P-P2)[0]) < 1e-15
@test norm(coeffs(Q*P-P2)[1]) < 1e-15
@test norm(coeffs(Q*P-P2)[2]) < 1e-15
