### Small μ
t = [2.0, 2.0, 2.0]
μ = [0.05, 0.1468, 0.5]
z = [2.0713428071028632 + 0.48208285060903094im,
    2.125956790502127 + 1.1121784646136845im,
    2.3122739258833174 + 2.711309573638537im]

t_vals_exact, pole_locs_exact = tracking_poles_exact(z, t, μ)
t_vals_saddle, pole_locs_saddle = tracking_poles_saddle(z, t, μ)

## Plot 
____fig = Figure(fontsize=26, resolution=(800, 400))
ax1 = Axis(____fig[1, 1], xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
    title=L"(a)$:$ Trajectories", titlealign=:left,
    xticks=([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], [L"0", L"0.5", L"1", L"1.5", L"2", L"2.5"]),
    yticks=([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], [L"0", L"0.5", L"1", L"1.5", L"2", L"2.5", L"3"]),
    width=400, height=400)
ax2 = Axis(____fig[1, 2], xlabel=L"t", ylabel=L"\mathrm{Im}(z)",
    title=L"(b)$:$ Distances", titlealign=:left,
    xticks=([0.0, 0.5, 1.0, 1.5, 2.0], [L"0", L"0.5", L"1", L"1.5", L"2"]),
    yticks=([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], [L"0", L"0.5", L"1", L"1.5", L"2", L"2.5", L"3"]),
    width=500, height=400)
for (t, z) in zip(t_vals_exact, pole_locs_exact)
    if z == pole_locs_exact[3]
        global ____a1
    end
    ____a1 = lines!(ax1, real.(z), imag.(z), color=:blue, linewidth=3, linestyle=:solid,
        label=L"$ $Exact", marker=:star4, markersize=7)
    lines!(ax2, t, imag.(z), color=:blue, linewidth=3,
        marker=:star4, markersize=7, linestyle=:solid)
end
for (t, z) in zip(t_vals_saddle, pole_locs_saddle)
    if z == pole_locs_saddle[3]
        global ____a2
    end
    ____a2 = lines!(ax1, real.(z), imag.(z), color=:red, linestyle=:dash,
        linewidth=3, label=L"$ $Saddle Point", marker=:rtriangle, markersize=7)
    lines!(ax2, t, imag.(z), color=:red, linewidth=3, linestyle=:dash,
        marker=:rtriangle, markersize=7)
end
ylims!(ax1, 0, 3)
ylims!(ax2, 0, 3)
xlims!(ax1, 0, 2.5)
xlims!(ax2, 0, 2)
text!(ax1, L"\mu = 0.5", position=(1.0, 1.86), rotation=0.45, color=:black, textsize=36.4)
text!(ax1, L"\mu = 0.1468", position=(1.0, 1.1), rotation=0.0, color=:black, textsize=36.4)
text!(ax1, L"\mu = 0.05", position=(0.5, 0.45), rotation=-0.20, color=:black, textsize=36.4)
text!(ax2, L"\mu = 0.5", position=(1.0, 2.25), rotation=0.3, color=:black, textsize=36.4)
text!(ax2, L"\mu = 0.1468", position=(1.0, 1.1), rotation=0.0, color=:black, textsize=36.4)
text!(ax2, L"\mu = 0.05", position=(0.5, 0.4), rotation=-0.20, color=:black, textsize=36.4)

include("aaa_tracking_poles.jl")

save("$FIGURES/tracking_pole_small_mu.$EXTENSION", ____fig)

### Larger μ
t = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
μ = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
z = [5.5 + 5.0im, 6.0 + 8.7im, 6.3 + 14.0im, 7.2 + 22.1im, 9.5 + 34.8im, 12.2 + 52.9im]
t_vals_exact, pole_locs_exact = tracking_poles_exact(z, t, μ)
t_vals_saddle, pole_locs_saddle = tracking_poles_saddle(z, t, μ)

fig = Figure(fontsize=26, resolution=(800, 400))
ax1 = Axis(fig[1, 1], xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
    title=L"(a)$:$ Trajectories", titlealign=:left,
    xticks=([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [L"0", L"1", L"2", L"3", L"4", L"5"]),
    yticks=([0.0, 5.0, 10.0, 15.0, 20.0, 25.0], [L"0", L"5", L"10", L"15", L"20", L"25"]))
ax2 = Axis(fig[1, 2], xlabel=L"t", ylabel=L"\mathrm{Im}(z)",
    title=L"(b)$:$ Distances", titlealign=:left,
    xticks=([0.0, 5.0, 10.0], [L"0", L"5", L"10"]),
    yticks=([0.0, 20.0, 40.0, 60.0], [L"0", L"20", L"40", L"60"]))
colors = [:darkblue, :red, :orange, :purple, :green, :magenta]
for (j, (t, z)) in enumerate(tuple.(t_vals_exact, pole_locs_exact))
    lines!(ax1, real.(z), imag.(z), color=colors[j], linewidth=3, label=L"%$(μ[j])")
    lines!(ax2, t, imag.(z), color=colors[j], linewidth=3, label=L"%$(μ[j])")
end
xlims!(ax1, 0, 5)
ylims!(ax1, 0, 25)
fig[0, 1:2] = Legend(fig, ax1, L"\mu", orientation=:horizontal, labelsize=26, titlesize=26, titleposition=:left)

save("$FIGURES/tracking_pole_larger_mu.$EXTENSION", fig)