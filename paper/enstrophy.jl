μ = μ_ENSTROPHY
t = T_ENSTROPHY
E = ENSTROPHY

indices, vals = findmaxima(E)
tvals = Vector{Float64}([])
peaks = Vector{Float64}([])
pt_idx = Vector{Int64}([])
pt_idx2 = Vector{Int64}([])
breakdown_μ = []
μ_idx = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250]
for (i, (idx, val)) in enumerate(zip(indices, vals))
    if isnan(val)
        push!(breakdown_μ, μ[i])
    else
        if (i ∈ μ_idx) && (i ≠ 1)
            push!(pt_idx, idx)
            push!(pt_idx2, i)
        end
        push!(tvals, t[idx])
        push!(peaks, val)
    end
end
breakdown_μ = breakdown_μ[1] # Not so good

fig = Figure(fontsize=33, resolution=(846.01953f0, 503.74707f0))
ax = Axis(fig[1, 1], width=600, height=400,
    xlabel=L"t", ylabel=L"E(t)",
    xticks=(0:5:5, [L"0", L"5"]),
    yticks=(0:0.5:2, [L"0", L"0.5", L"1", L"1.5", L"2"]))
μ_idx = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250]
colors = cgrad(LINSPECER_12_J, μ[μ_idx]; categorical=false)
lines!(ax, t, E[:, 1], color=colors[1], linestyle=:dash)
[lines!(ax, t, E[:, j], color=colors[i+1]) for (i, j) in enumerate(μ_idx[2:end])]
scatter!(ax, t[pt_idx], vals[pt_idx2], color=:red, markersize=4)
xlims!(ax, 0, 5)
ylims!(ax, 0, 2)
Colorbar(fig[1, 2], limits=(0.0, 0.5), colormap=colors, label=L"\mu", vertical=true,
    ticks=([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], [L"0.0", L"0.1", L"0.2", L"0.3", L"0.4", L"0.5"]))
save("$FIGURES/enstrophy_plot.$EXTENSION", fig)
