x = LinRange(-75, 75, 5000)
t = [10.0, 25.0, 100.0, 500.0]
μ = [0.1, 1.0]
vals = viscous_solution_large_time(x, t, μ)
exvals = viscous_solution(x, t, μ)

fig = Figure(fontsize=41)

ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"u(x,t)",
    title=L"(a): $\mu = 0.1$", titlealign=:left, width=600, height=600,
    xticks=([-75, -50, -25, 0, 25, 50, 75], [L"-75", L"-50", L"-25", L"0", L"25", L"50", L"75"]),
    yticks=([0.0, 0.2, 0.4, 0.6, 0.8], [L"0", L"0.2", L"0.4", L"0.6", L"0.8"]),aspect=1)
lines!(ax, x, vals[:, 1, 1], color=:blue, linewidth=3)
lines!(ax, x, vals[:, 2, 1], color=:red, linewidth=3)
lines!(ax, x, vals[:, 3, 1], color=:darkgreen, linewidth=3)
lines!(ax, x, vals[:, 4, 1], color=:magenta, linewidth=3)
lines!(ax, x, exvals[:, 1, 1], color=:blue, linestyle=:dash, linewidth=3)
lines!(ax, x, exvals[:, 2, 1], color=:red, linestyle=:dash, linewidth=3)
lines!(ax, x, exvals[:, 3, 1], color=:darkgreen, linestyle=:dash, linewidth=3)
lines!(ax, x, exvals[:, 4, 1], color=:magenta, linestyle=:dash, linewidth=3)

ax = Axis(fig[2, 1], xlabel=L"x", ylabel=L"u(x,t)",
    title=L"(d): $\mu = 1$", titlealign=:left, width=600, height=600,
    xticks=([-50, -25, 0, 25, 50], [L"-50", L"-25", L"0", L"25", L"50"]),
    yticks=([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], [L"0", L"0.05", L"0.1", L"0.15", L"0.2", L"0.25", L"0.3"]),aspect=1)
lines!(ax, x, vals[:, 1, 2], color=:blue, linewidth=3)
lines!(ax, x, vals[:, 2, 2], color=:red, linewidth=3)
lines!(ax, x, vals[:, 3, 2], color=:darkgreen, linewidth=3)
lines!(ax, x, vals[:, 4, 2], color=:magenta, linewidth=3)
lines!(ax, x, exvals[:, 1, 2], color=:blue, linestyle=:dash, linewidth=3)
lines!(ax, x, exvals[:, 2, 2], color=:red, linestyle=:dash, linewidth=3)
lines!(ax, x, exvals[:, 3, 2], color=:darkgreen, linestyle=:dash, linewidth=3)
lines!(ax, x, exvals[:, 4, 2], color=:magenta, linestyle=:dash, linewidth=3)

legendentries = OrderedDict(L"10" => LineElement(linestyle=nothing, linewidth=2.0, color=:blue),
    L"25" => LineElement(linestyle=nothing, linewidth=2.0, color=:red),
    L"50" => LineElement(linestyle=nothing, linewidth=2.0, color=:darkgreen),
    L"100" => LineElement(linestyle=nothing, linewidth=2.0, color=:magenta))
Legend(fig[1:2, 0], [values(legendentries)...], [keys(legendentries)...], L"$t$", orientation=:vertical, labelsize=36, titlesize=36, titleposition=:top)

x = -12:0.01:12 |> collect
y = 0:0.01:24 |> collect
μ = [0.1, 1.0]
vals = viscous_solution_large_time_Ψ(x, y, μ)
ρ = 3:0.01:88
n = 1:10000
poles = viscous_solution_large_time_Ψ_roots_θ(ρ, μ, 1)
poles_n = viscous_solution_large_time_Ψ_roots_ρ(n, μ, 1)
poles_2 = viscous_solution_large_time_Ψ_roots_θ(ρ, μ, 2)
poles_n_2 = viscous_solution_large_time_Ψ_roots_ρ(n, μ, 2)
portrait!(fig, x, y, vals[:, :, 1], 1, 2, width=600, height=600, nist=NIST,
    xlabel=L"\mathrm{Re}(\eta)", ylabel=L"\mathrm{Im}(\eta)",
    xticks=([-12, -8, -4, 0, 4, 8, 12],
        [L"-12", L"-8", L"-4", L"0", L"4", L"8", L"12"]),
    yticks=([0, 4, 8, 12, 16, 20, 24],
        [L"0", L"4", L"8", L"12", L"16", L"20", L"24"]),
    title=L"(b): $\mu = 0.1$", titlealign=:left,aspect=1)
#abline!(fig.content[end], 0, 1, color=:black)
#abline!(fig.content[end], 0, -1, color=:black)
lines!(fig.content[end], real(poles[:, 1]), imag(poles[:, 1]), color=:black, linewidth=6)
scatter!(fig.content[end], real(poles_n[:, 1]), imag(poles_n[:, 1]), color=:white, markersize=16)
lines!(fig.content[end], real(poles_2[:, 1]), imag(poles_2[:, 1]), color=:black, linewidth=6)
scatter!(fig.content[end], real(poles_n_2[:, 1]), imag(poles_n_2[:, 1]), color=:white, markersize=16)
xlims!(fig.content[end], x[1], x[end])
ylims!(fig.content[end], y[1], y[end])
xlims!(fig.content[end], -12, 12)
ylims!(fig.content[end], 0, 24)
portrait!(fig, x, y, vals[:, :, 2], 2, 2, width=600, height=600, nist=NIST,
    xlabel=L"\mathrm{Re}(\eta)", ylabel=L"\mathrm{Im}(\eta)",
    xticks=([-12, -8, -4, 0, 4, 8, 12],
        [L"-12", L"-8", L"-4", L"0", L"4", L"8", L"12"]),
    yticks=([0, 4, 8, 12, 16, 20, 24],
        [L"0", L"4", L"8", L"12", L"16", L"20", L"24"]),
    title=L"(e): $\mu = 1$", titlealign=:left,aspect=1)
#abline!(fig.content[end], 0, 1, color=:black)
#abline!(fig.content[end], 0, -1, color=:black)
lines!(fig.content[end], real(poles[:, 2]), imag(poles[:, 2]), color=:black, linewidth=6)
scatter!(fig.content[end], real(poles_n[:, 2]), imag(poles_n[:, 2]), color=:white, markersize=16)
lines!(fig.content[end], real(poles_2[:, 2]), imag(poles_2[:, 2]), color=:black, linewidth=6)
scatter!(fig.content[end], real(poles_n_2[:, 2]), imag(poles_n_2[:, 2]), color=:white, markersize=16)
xlims!(fig.content[end], -12, 12)
ylims!(fig.content[end], 0, 24)

x = -15:0.1:15 |> collect
y = 0:0.1:15 |> collect
μ = [0.1, 1.0]
vals = viscous_solution_large_time_Ψ(x, y, μ)
vals2 = viscous_solution_large_time_Ψ(x, μ)
landscape!(fig, x, y, vals[:, :, 1], 1, 5, width=600, height=600, nist=NIST,
    xlabel=L"\mathrm{Re}(\eta)", ylabel=L"\mathrm{Im}(\eta)", zlabel=L"|\Psi(\eta)|",
    title=L"(c): $\mu = 0.1$", titlealign=:left, zlims=(0, 10),
    xticks=([-15, 0, 15], [L"-15", L"0", L"15"]),
    yticks=([0, 5, 10, 15], [L"0", L"5", L"10", L"15"]),
    zticks=([0, 5, 10], [L"0", L"5", L"10"]))
lines!(fig.content[end], x, zeros(length(x)), vals2[:, 1], color=:black, linewidth=5)
landscape!(fig, x, y, vals[:, :, 2], 2, 5, width=600, height=600, nist=NIST,
    xlabel=L"\mathrm{Re}(\eta)", ylabel=L"\mathrm{Im}(\eta)", zlabel=L"|\Psi(\eta)|",
    title=L"(f): $\mu = 1$", titlealign=:left, zlims=(0, 10),
    xticks=([-15, 0, 15], [L"-15", L"0", L"15"]),
    yticks=([0, 5, 10, 15], [L"0", L"5", L"10", L"15"]),
    zticks=([0, 5, 10], [L"0", L"5", L"10"]))
lines!(fig.content[end], x, zeros(length(x)), vals2[:, 2], color=:black, linewidth=5)
resize_to_layout!(fig)

save("$FIGURES/large_time_solution_plots.$EXTENSION", fig)

fig

