# The code below is one way to do this in Julia. It's not as accurate though... ApproxFun.jl would be nice to use. For now we paste some results from Chebfun.
#= 
t = 2.0
μ = [0.05, 0.1468, 0.5]
#z₁1 = 3.6149423868439254 + 0.45654346819226965im
#z₁2 = 3.605398629392598 + 1.3294133122075977im
#z₁3 = 3.7913588390007855 + 3.7531544895022875im
#z = [z₁1, z₁2, z₁3]
z = [2.0714 + 0.4821im, 2.1270 + 1.225im, 2.3123 + 2.7113im]
## Track the exact poles
N = 20000
Δt = t / N
t_vals = 0:Δt:t
pole_locs_exact = Vector{Vector{ComplexF64}}(undef, 3)
for j in 1:3
    _, pole_locs_exact[j] = tracking_poles_exact(z[j], t[1], μ[j])
end

## AAA 
t_vals_aaa = []
pole_locs_aaa = []
#t_max = 5.0
t_max = 2.0;
t_min = 1e-6
for i in 1:3
    taaa, paaa = tracking_poles_aaa(z[i], μ[i]; t_max, t_min, L=15)
    push!(t_vals_aaa, taaa)
    push!(pole_locs_aaa, paaa)
end

## Plot 
for (t, z) in zip(t_vals_aaa, pole_locs_aaa)
    scatter!(____fig.content[1], real.(z[1:25:end]), imag.(z[1:25:end]), color=:green, markersize=5)
    scatter!(____fig.content[2], t[1:25:end], imag.(z[1:25:end]), markersize=5, color=:green)
end
=#

# The constants below are obtained in the MATLAB script "burger_tracking_poles.m".
skip = 25
for j in 1:3
    if j > 2
        global ____a3
    end
    ____a3 = scatter!(____fig.content[1], real.(AAA_POLES[1:skip:end-1, j]), imag.(AAA_POLES[1:skip:end-1, j]),
        color=:green, markersize=8, label=L"$ $AAA", linewidth=3, linestyle=:dashdot)
    scatter!(____fig.content[2], T_AAA_POLES[1:skip:end-1], imag.(AAA_POLES[1:skip:end-1, j]),
        markersize=8, color=:green, linewidth=3, linestyle=:dashdot)
end

axislegend(ax2, [____a1, ____a2, ____a3], [L"$ $Exact", L"$ $Saddle Point", L"$ $AAA"], position=:lt)
resize_to_layout!(____fig)