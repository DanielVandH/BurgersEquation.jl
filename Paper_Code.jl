using LinearAlgebra
using CairoMakie
using FastGaussQuadrature
using HypergeometricFunctions
using ForwardDiff
using SpecialFunctions
using ApproxFun
using Colors
using ComplexPortraits
using BurgersEquation

x = collect(-2:0.01:2) 
y = collect(-2:0.01:2)
Z = [1/(1+(x+im*y)^2) for x in x, y in y]
fig = Figure()
ax = Axis(fig[1, 1], zlabel = L"|1/(1+z^2)|")
surface!(ax, x, y, abs.(Z))
zlims!(ax, 0, 5)