using BurgersEquation
using CairoMakie
using LinearAlgebra
using FastGaussQuadrature
using Contour
using OrderedCollections
using Peaks
using DelimitedFiles
using ArbNumerics
using JLD2

include("constants.jl")
include("introduction.jl")
include("exact_solution.jl")
include("colour_wheel.jl")
include("parabolic_cylinder_function.jl")
include("steepest_descent_preliminaries.jl")
include("tracking_poles.jl")
include("balancing_advection_diffusion.jl")
include("large_time_solution.jl")
include("large_time_solution_roots.jl")
include("aaa.jl")
include("aaa_closest_poles.jl")
# include("aaa_tracking_poles.jl") < --- See tracking_pole_saddle_mu.jl
include("small_to_large_transition.jl")

## Some extra figures added after review
# 1/(1+x^2)^2
include("initial_condition_1o1px2a2.jl")

# 1/(1+x^2)^(1/2)
include("initial_condition_1osqrt1px2.jl")

# Inviscid solutions for 1/(1+x^2), 1/(1+x^2)^2, and 1/(1+x^2)^(1/2)
include("inviscid_solutions.jl")