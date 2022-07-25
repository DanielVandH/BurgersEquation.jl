using BurgersEquation
using CairoMakie
using LinearAlgebra
using FastGaussQuadrature
using Contour
using OrderedCollections
using Peaks
#using BaryRational
using DelimitedFiles

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