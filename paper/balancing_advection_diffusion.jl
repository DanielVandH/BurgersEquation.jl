## Measuring slopes on the real line
μ = [0.0, LinRange(1e-2, 0.5, 250)...]
t = LinRange(1e-6, 5.0, 1000)
slopes = maximum_slope(t, μ)
indices, vals = findmaxima(slopes)
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
breakdown_μ = breakdown_μ[1]

## Tracking the closest pole
t_min = 1e-6
t_max = 5.0
z₀ = 3.7913503689060826 + 3.7531275068755505im
t_vals, pole_locs = tracking_poles_exact(z₀, μ; t_max=t_max, t_min=t_min)

## Plots 
# Real line
fig = Figure(fontsize=26, resolution=(800, 400))
ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"Maximum$ $ absolute slope", title=L"(a):$ $ Slope analysis", titlealign=:left,
    xticks=([0.0, 5.0, 10.0], [L"0", L"5", L"10"]),
    yticks=([0.0, 2.0, 4.0, 6.0, 8.0], [L"0", L"2", L"4", L"6", L"8"]))
colors = cgrad(LINSPECER_12_J, μ[μ_idx]; categorical=false)
lines!(ax, t[1:5:283], slopes[1:5:283, 1], color=colors[1], linestyle=:dash)
[lines!(ax, t, slopes[:, j], color=colors[i+1]) for (i, j) in enumerate(μ_idx[2:end])]
scatter!(ax, t[pt_idx], vals[pt_idx2], color=:red, markersize=4)
xlims!(ax, 0, t_max)

# Closest pole
ax = Axis(fig[1, 2], xlabel=L"t", ylabel=L"Distance$ $ to real line", title=L"(b):$ $ Pole locations", titlealign=:left,
    xticks=([0.0, 5.0, 10.0], [L"0", L"5", L"10"]),
    yticks=([0.0, 2.0, 4.0, 6.0], [L"0", L"2", L"4", L"6"]))
zt_vals = Vector{ComplexF64}([])
for t in t_vals
    _, zs = inviscid_singularities(t)
    if t < 8sqrt(3) / 9
        push!(zt_vals, zs[1])
    else
        push!(zt_vals, NaN + im * NaN)
    end
end
lines!(ax, t_vals, imag.(zt_vals), color=colors[1], linestyle=:dash)
[lines!(ax, t_vals, imag.(pole_locs[j, :]), color=colors[i+1]) for (i, j) in enumerate(μ_idx[2:end])]
xlims!(ax, 0, t_max)
minima_t = Vector{Float64}([])
minima_p = Vector{Float64}([])
plot_t = Vector{Float64}([])
plot_p = Vector{Float64}([])
break_i = Int64[]
μ_dt = Float64[]
for i in 2:length(μ)
    a, b = findminima(imag(pole_locs[i, 2:end]))
    try
        b, idx = findmin(b)
        a = a[idx]
        if i > 2 && abs(t[a] - minima_t[end]) > 0.1
            push!(break_i, i)
            break
        end
        push!(minima_t, t[a])
        push!(minima_p, b)
        if i ∈ μ_idx
            push!(plot_t, t[a])
            push!(plot_p, b)
            push!(μ_dt, μ[i])
        end
        if a == 1 || a == length(t)
            push!(break_i, i)
            continue
        end
    catch
        push!(break_i, i)
        continue
    end
end
μb = μ[break_i]
scatter!(ax, plot_t, plot_p, color=:red, markersize=4)
ylims!(ax, 0, 4)
Colorbar(fig[1, 3], limits=(0.0, 0.5), colormap=colors, label=L"\mu", vertical=true,
    ticks=([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], [L"0.0", L"0.1", L"0.2", L"0.3", L"0.4", L"0.5"]))
save("$FIGURES/advection_diffusion_balance.$EXTENSION", fig)

## When are there no more zeros in the lower complex plane in the small-time solution?
x = 1.0:0.001:1.2
y = -0.6:0.001:-0.4
μ = 0.05
Φ₀_vals = Φ₀(x, y, μ)
#fig = Figure()
#portrait!(fig, x, y, Φ₀_vals, 1, 1, nist=NIST)
z₀ = 1.1 - 0.5im
pole_locs, μ_vals = tracking_poles_Φ₀(z₀, μ)
idx = findfirst(imag(pole_locs) .≥ 0.0)
μˢ = μ_vals[idx]
idxs = 1:380:95000
colors = cgrad(LINSPECER_250_J, μ_vals[idxs], categorical=false)

fig = Figure(fontsize=26, resolution=(800, 400))
ax = Axis(fig[1, 1], xlabel=L"\mu", ylabel=L"\mathrm{Im}(\xi)",
    title=L"(a):$ $ Pole distances", titlealign=:left,
    xticks=([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [L"0.0", L"0.2", L"0.4", L"0.6", L"0.8", L"1"]),
    yticks=([-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5], [L"-0.5", L"0", L"0.5", L"1", L"1.5", L"2", L"2.5"]))
[scatter!(ax, [μ_vals[idxs[i]]], [imag(pole_locs[idxs[i]])], color=colors[i], markersize=4) for i in 1:250]
hlines!(ax, [0.0], color=:black)
vlines!(ax, [μˢ], color=:black)
xlims!(ax, 0, 1)
ylims!(ax, -0.5, 2.5)

ax = Axis(fig[1, 2], xlabel=L"\mathrm{Re}(\xi)", ylabel=L"\mathrm{Im}(\xi)",
    title=L"(b):$ $ Pole locations", titlealign=:left,
    xticks=([1.0, 1.2, 1.4, 1.6], [L"1", L"1.2", L"1.4", L"1.6"]),
    yticks=([-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5], [L"-0.5", L"0", L"0.5", L"1", L"1.5", L"2", L"2.5"]))
[scatter!(ax, [real(pole_locs[idxs[i]])], [imag(pole_locs[idxs[i]])], color=colors[i], markersize=4) for i in 1:250]
hlines!(ax, [0.0], color=:black)
vlines!(ax, [real(pole_locs)[idx]], color=:black)
xlims!(ax, 1, 1.6)
ylims!(ax, -0.5, 2.5)
Colorbar(fig[1, 3], limits=(0.0, 1.0), colormap=colors, label=L"\mu", vertical=true,
    ticks=([0.0, 0.25, 0.5, 0.75, 1.0], [L"0", L"0.25", L"0.5", L"0.75", L"1"]))

save("$FIGURES/advection_diffusion_balance_phi0_solution.$EXTENSION", fig)
