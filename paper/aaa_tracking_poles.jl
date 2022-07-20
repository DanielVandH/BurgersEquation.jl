t = 5.0
μ = [0.05, 0.1468, 0.5]
z₁1 = 3.6149423868439254 + 0.45654346819226965im
z₁2 = 3.605398629392598 + 1.3294133122075977im
z₁3 = 3.7913588390007855 + 3.7531544895022875im
z = [z₁1, z₁2, z₁3]

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
t_max = 5.0
t_min = 1e-6
for i in 1:3
    taaa, paaa = tracking_poles_aaa(z[i], μ[i]; t_max, t_min, L=15)
    push!(t_vals_aaa, taaa)
    push!(pole_locs_aaa, paaa)
end

## Plot 
for (t, z) in zip(t_vals_aaa, pole_locs_aaa)
    scatter!(____fig.content[1], real.(z[1:25:end]), imag.(z[1:25:end]), color=:green, markersize=5)
    scatter!(____fig.content[2], t[1:25:end], imag.(z[1:25:end]), markersize = 5, color=:green)
end

legendentries = OrderedDict(L"$ $Exact" => LineElement(linestyle=nothing, linewidth=2.0, color=:blue),
    L"$ $Saddle Point" => LineElement(linestyle=:dash, linewidth=2.0, color=:red),
    L"$ $AAA" => LineElement(linestyle=:dash, linewidth=2.0, color=:green))
axislegend(ax2, [values(legendentries)...], [keys(legendentries)...], position=:lt)
resize_to_layout!(____fig)
