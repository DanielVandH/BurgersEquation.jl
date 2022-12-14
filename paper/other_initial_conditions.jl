###########################################################################
## 1/(1+x^2)^2
###########################################################################
## Exact solution plots for O(1) time
x = LinRange(-5, 5, 250)
y = LinRange(0, 5, 250)
_tₛ = 54sqrt(5) / 125
t = [0.1, 0.5, 0.75, _tₛ]
μ = [0.1, 0.5, 1.0]
u = viscous_solution(x, y, t, μ; ic=2)
u_real = viscous_solution(x, t, μ; ic=2)
time_labels = ["0.1", "0.5", "0.75", "t_s"]
plot_layout = [[(1, 1), (1, 2)], [(1, 3), (1, 4)], [(2, 1), (2, 2)], [(2, 3), (2, 4)]]
alphabets = [['a', 'b'], ['c', 'd'], ['e', 'f'], ['g', 'h']]
for k = 1:length(μ)
    fig = Figure(fontsize=38, resolution=(2000, 800))
    for j in 1:length(t)
        portrait!(fig, Float64.(x), Float64.(y), Complex{Float64}.(u[:, :, j, k]), plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
            title=L"(%$(alphabets[j][1])): $t = %$(time_labels[j])$", titlealign=:left,
            xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            yticks=([0.0, 2.0, 4.0], [L"0", L"2", L"4"]),
            width=600, height=400)
        scatter!(fig.content[end], [0.0], [1.0], color=:black, markersize=12)
        landscape!(fig, Float64.(x), Float64.(y), Complex{Float64}.(u[:, :, j, k]), plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
            title=L"(%$(alphabets[j][2])): $t = %$(time_labels[j])$", titlealign=:left, zlims=(0, 3),
            xticks=([-4.0, 0.0, 4.0], [L"-4", L"0", L"4"]),
            yticks=([0.0, 2.0, 4.0], [L"0", L"2", L"4"]),
            zticks=([0.0, 1, 2, 3], [L"0", L"1", L"2", L"3"]),
            width=600, height=400)
        lines!(fig.content[end], Float64.(x), zeros(length(y)), Complex{Float64}.(u_real[:, j, k]), color=:black, linewidth=5)
    end
    resize_to_layout!(fig)
    fig
    save("$FIGURES/viscous_solutions_initial_condition_1o1px2a2_$k.$EXTENSION", fig)
end

## Exact solution plots for small time
x = LinRange(-0.015, 0.015, 250)
y = LinRange(0.980, 1.02, 250)
μ1 = 2.0
μ2 = 1.0
μ3 = 0.5
μ4 = 0.1
μ = [μ1, μ2, μ3, μ4]
t = 1e-2
u = viscous_solution(x, y, [t], μ; ic=2)
fig = Figure(fontsize=38, resolution=(2449.9102f0, 604.33514f0))
plot_layout = [[(1, 1), (2, 1), (3, 1), (4, 1)], [(1, 2), (2, 2), (3, 2), (4, 2)],
    [(1, 3), (2, 3), (3, 3), (4, 3)], [(1, 4), (2, 4), (3, 4), (4, 4)]]
reverse!(plot_layout)
for j = 1:4
    portrait!(fig, Float64.(x), Float64.(y), Complex{Float64}.(u[:, :, 1, j]), plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(ALPHABET[4-j+1])): $\mu = %$(Float64(μ[j]))$, $t = 10^{-2}$", titlealign=:left,
        xlabel=L"x", ylabel=L"y", width=450, height=450,
        xticks=([-0.010, 0.0, 0.010], [L"-0.010", L"0", L"0.010"]),
        yticks=([0.980, 1.0, 1.020], [L"0.980", L"1", L"1.020"]))
end
save("$FIGURES/similarity_solutions_initial_condition_1o1px2a2.$EXTENSION", fig) # These aren't similarity solutions, but we'll just keep the prefix to make it easier to compare figures later with the other small time results

## Exact solution plots for large time 
# Larger mu
x = LinRange(-200.0, 200.0, 250)
y = LinRange(0, 200.0, 250)
t = [100.0, 250.0, 500.0, 1000.0]
μ = [1.0]
u_vals = viscous_solution(x, y, t, μ; ic=2)
fig = Figure(fontsize=38, resolution=(2970.0f0, 551.9453f0))
portrait!(fig, x, y, u_vals[:, :, 1, 1], 1, 1;
    xlabel=L"x", ylabel=L"y", title=L"(a): $t = 100$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 2, 1], 1, 2;
    xlabel=L"x", ylabel=L"y", title=L"(b): $t = 250$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 3, 1], 1, 3;
    xlabel=L"x", ylabel=L"y", title=L"(c): $t = 500$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 4, 1], 1, 4;
    xlabel=L"x", ylabel=L"y", title=L"(d): $t = 1000$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
save("$FIGURES/small_to_large_transition_initial_condition_1o1px2a2.$EXTENSION", fig)

# Smaller mu
x = LinRange(-200.0, 200.0, 250)
y = LinRange(0, 200.0, 250)
t = [100.0, 250.0, 500.0, 1000.0]
μ = [0.1]
u_vals = viscous_solution(x, y, t, μ; ic=2)
fig = Figure(fontsize=38, resolution=(2970.0f0, 551.9453f0))
portrait!(fig, x, y, u_vals[:, :, 1, 1], 1, 1;
    xlabel=L"x", ylabel=L"y", title=L"(a): $t = 100$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 2, 1], 1, 2;
    xlabel=L"x", ylabel=L"y", title=L"(b): $t = 250$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 3, 1], 1, 3;
    xlabel=L"x", ylabel=L"y", title=L"(c): $t = 500$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 4, 1], 1, 4;
    xlabel=L"x", ylabel=L"y", title=L"(d): $t = 1000$",
    width=600, height=400,
    xticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    yticks=([-200, -100, 0, 100, 200], [L"-200", L"-100", L"0", L"100", L"200"]),
    titlealign=:left)
save("$FIGURES/small_to_large_transition_initial_condition_1o1px2a2_smaller_mu.$EXTENSION", fig)

## Slope analysis
μ = [0.0, LinRange(1e-2, 0.5, 250)...]
t = LinRange(1e-6, 5.0, 5000)
slopes = maximum_slope(t, μ; ic=2, alg=LBFGS(), init=15.0, num_nodes=500)
x = LinRange(-10, 10, 2000)
f = y -> 1 / (1 + y^2)^2
f′ = y -> -4y / (1 + y^2)^3
u = inviscid_solution(f, f′, x, t)
uₓ = similar(u)
for j in eachindex(t)
    for i in eachindex(x)
        uₓ[i, j] = f′(x[i] - u[i, j] * t[j]) / (1 + t[j] * f′(x[i] - u[i, j] * t[j]))
    end
end
auₓ = abs.(uₓ)
maxauₓ = zeros(length(t))
for i in eachindex(maxauₓ)
    @views maxauₓ[i] = maximum(auₓ[:, i])
    if i > 1 && maxauₓ[i] < maxauₓ[i-1] - 0.7
        maxauₓ[i:end] .= NaN
        break
    end
end
slopes[:, 1] .= maxauₓ
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
t_min = 1e-6
t_max = 5.0
fig = Figure(fontsize=33, resolution=(1700, 400))
ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"Maximum$ $ absolute slope", title=L"(a):$ $ Slope analysis", titlealign=:left,
    xticks=([0.0, 5.0, 10.0], [L"0", L"5", L"10"]), height=400, width=600,
    yticks=([0.0, 2.0, 4.0, 6.0, 8.0], [L"0", L"2", L"4", L"6", L"8"]))
colors = cgrad(LINSPECER_12_J, μ[μ_idx]; categorical=false)
lines!(ax, t, slopes[:, 1], color=colors[1], linestyle=:dash)
[lines!(ax, t, slopes[:, j], color=colors[i+1]) for (i, j) in enumerate(μ_idx[2:end])]
scatter!(ax, t[pt_idx], vals[pt_idx2], color=:red, markersize=4)
xlims!(ax, 0, t_max)

# Enstrophy
μ = μ_ENSTROPHY_IC2
_t = T_ENSTROPHY_IC2
E = ENSTROPHY_IC2

f = x -> 1 / (1 + x^2)^2
f′ = x -> -4x / (1 + x^2)^3
u = (x, t) -> inviscid_solution(f, f′, x, t)
dudx = (x, t) -> ForwardDiff.derivative(x -> u(x, t), x)
dudx2 = (x, t) -> dudx(x, t)^2
nodes, weights = gausslegendre(600)
dudx2gh = (x, t) -> dudx(x, t)^2
Efnc = t -> gauss_legendre_finite(x -> dudx2gh(x, t), -5.0, 5.0, nodes, weights) / 2
Efnc0 = Efnc.(_t[_t.<0.86])
E[:, 1] .= NaN
E[_t.<0.86, 1] .= Efnc0

indices, vals = findmaxima(E)
_tvals = Vector{Float64}([])
_peaks = Vector{Float64}([])
_pt_idx = Vector{Int64}([])
_pt_idx2 = Vector{Int64}([])
_breakdown_μ = []
_μ_idx = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250]
for (i, (idx, val)) in enumerate(zip(indices, vals))
    if isnan(val)
        push!(_breakdown_μ, μ[i])
    else
        if (i ∈ μ_idx) && (i ≠ 1)
            push!(_pt_idx, idx)
            push!(_pt_idx2, i)
        end
        push!(_tvals, _t[idx])
        push!(_peaks, val)
    end
end
_breakdown_μ = _breakdown_μ[1] # Not so good
ax = Axis(fig[1, 2], width=600, height=400,
    xlabel=L"t", ylabel=L"E(t)",
    xticks=(0:5:5, [L"0", L"5"]),
    yticks=(0:3, [L"%$s" for s in 0:3]),
    title=L"(b):$ $ Enstrophy", titlealign=:left)
μ_idx = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250]
colors = cgrad(LINSPECER_12_J, μ[μ_idx]; categorical=false)
lines!(ax, _t, E[:, 1], color=colors[1], linestyle=:dash)
[lines!(ax, _t, E[:, j], color=colors[i+1]) for (i, j) in enumerate(μ_idx[2:end])]
scatter!(ax, _t[_pt_idx], vals[_pt_idx2], color=:red, markersize=4)
xlims!(ax, 0, 5)
ylims!(ax, 0, 3)
Colorbar(fig[1, 3], limits=(0.0, 0.5), colormap=colors, label=L"\mu", vertical=true,
    ticks=([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], [L"0.0", L"0.1", L"0.2", L"0.3", L"0.4", L"0.5"]))
resize_to_layout!(fig)
fig.content[1].yticks = (0:3:9, [L"%$s" for s in 0:3:9])
ylims!(fig.content[1], 0, 9)
save("$FIGURES/slope_enstrophy_plot_1o1px2a2.$EXTENSION", fig)

###########################################################################
## 1/(1+x^2)^(1/2)
###########################################################################
## Exact solution plots for O(1) time
x = LinRange(-5, 5, 250)
y = LinRange(0, 5, 250)
_tₛ = 3sqrt(3) / 2
t = [0.1, 0.5, 1.0, _tₛ]
μ = [0.1, 0.5, 1.0]
u = viscous_solution(x, y, t, μ; ic=3)
u_real = viscous_solution(x, t, μ; ic=3)
time_labels = ["0.1", "0.5", "0.75", "t_s"]
plot_layout = [[(1, 1), (1, 2)], [(1, 3), (1, 4)], [(2, 1), (2, 2)], [(2, 3), (2, 4)]]
alphabets = [['a', 'b'], ['c', 'd'], ['e', 'f'], ['g', 'h']]
for k = 1:length(μ)
    fig = Figure(fontsize=38, resolution=(2000, 800))
    for j in 1:length(t)
        portrait!(fig, Float64.(x), Float64.(y), Complex{Float64}.(u[:, :, j, k]), plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
            title=L"(%$(alphabets[j][1])): $t = %$(time_labels[j])$", titlealign=:left,
            xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            yticks=([0.0, 2.0, 4.0], [L"0", L"2", L"4"]),
            width=600, height=400)
        scatter!(fig.content[end], [0.0], [1.0], color=:black, markersize=12)
        landscape!(fig, Float64.(x), Float64.(y), Complex{Float64}.(u[:, :, j, k]), plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
            title=L"(%$(alphabets[j][2])): $t = %$(time_labels[j])$", titlealign=:left, zlims=(0, 3),
            xticks=([-4.0, 0.0, 4.0], [L"-4", L"0", L"4"]),
            yticks=([0.0, 2.0, 4.0], [L"0", L"2", L"4"]),
            zticks=([0.0, 1, 2, 3], [L"0", L"1", L"2", L"3"]),
            width=600, height=400)
        lines!(fig.content[end], Float64.(x), zeros(length(y)), Complex{Float64}.(u_real[:, j, k]), color=:black, linewidth=5)
    end
    resize_to_layout!(fig)
    fig
    save("$FIGURES/viscous_solutions_initial_condition_1osqrt1px2_$k.$EXTENSION", fig)
end

## Exact solution plots for small time
x = LinRange(-0.01, 0.01, 250)
y = LinRange(0.990, 1.01, 250)
μ1 = 2.0
μ2 = 1.0
μ3 = 0.5
μ4 = 0.1
μ = [μ1, μ2, μ3, μ4]
t = 1e-6
u = viscous_solution(x, y, [t], μ; ic=3)
fig = Figure(fontsize=38, resolution=(2449.9102f0, 604.33514f0))
plot_layout = [[(1, 1), (2, 1), (3, 1), (4, 1)], [(1, 2), (2, 2), (3, 2), (4, 2)],
    [(1, 3), (2, 3), (3, 3), (4, 3)], [(1, 4), (2, 4), (3, 4), (4, 4)]]
reverse!(plot_layout)
for j = 1:4
    portrait!(fig, Float64.(x), Float64.(y), Complex{Float64}.(u[:, :, 1, j]), plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(ALPHABET[4-j+1])): $\mu = %$(Float64(μ[j]))$, $t = 10^{-6}$", titlealign=:left,
        xlabel=L"x", ylabel=L"y", width=450, height=450,
        xticks=([-0.010, 0.0, 0.010], [L"-0.010", L"0", L"0.010"]),
        yticks=([0.990, 1.0, 1.010], [L"0.990", L"1", L"1.010"]))
end
save("$FIGURES/similarity_solutions_initial_condition_1osqrt1px2.$EXTENSION", fig) # These aren't similarity solutions, but we'll just keep the prefix to make it easier to compare figures later with the other small time results

## Exact solution plots for large time 
# Larger mu 
x = LinRange(-150.0, 150.0, 250)
y = LinRange(0, 150.0, 250)
t = [100.0, 250.0, 500.0, 1000.0]
μ = [1.0]
u_vals = viscous_solution(x, y, t, μ; ic=3)
fig = Figure(fontsize=38, resolution=(2970.0f0, 551.9453f0))
portrait!(fig, x, y, u_vals[:, :, 1, 1], 1, 1;
    xlabel=L"x", ylabel=L"y", title=L"(a): $t = 100$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 2, 1], 1, 2;
    xlabel=L"x", ylabel=L"y", title=L"(b): $t = 250$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 3, 1], 1, 3;
    xlabel=L"x", ylabel=L"y", title=L"(c): $t = 500$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 4, 1], 1, 4;
    xlabel=L"x", ylabel=L"y", title=L"(d): $t = 1000$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
save("$FIGURES/small_to_large_transition_initial_condition_1osqrt1px2.$EXTENSION", fig)

# Smaller mu
x = LinRange(-150.0, 150.0, 250)
y = LinRange(0, 150.0, 250)
t = [100.0, 250.0, 500.0, 1000.0]
μ = [0.1]
u_vals = viscous_solution(x, y, t, μ; ic=3)
fig = Figure(fontsize=38, resolution=(2970.0f0, 551.9453f0))
portrait!(fig, x, y, u_vals[:, :, 1, 1], 1, 1;
    xlabel=L"x", ylabel=L"y", title=L"(a): $t = 100$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 2, 1], 1, 2;
    xlabel=L"x", ylabel=L"y", title=L"(b): $t = 250$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 3, 1], 1, 3;
    xlabel=L"x", ylabel=L"y", title=L"(c): $t = 500$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
portrait!(fig, x, y, u_vals[:, :, 4, 1], 1, 4;
    xlabel=L"x", ylabel=L"y", title=L"(d): $t = 1000$",
    width=600, height=400,
    xticks=([-150,-100,-50,0,50,100,150], [L"-150",L"-100",L"-50",L"0",L"50",L"100",L"150"]),
    yticks=([0,50,100,150], [L"0",L"50",L"100",L"150"]),
    titlealign=:left)
save("$FIGURES/small_to_large_transition_initial_condition_1osqrt1px2_smaller_mu.$EXTENSION", fig)

## Slope analysis
μ = [0.0, LinRange(1e-2, 0.5, 250)...]
t = LinRange(1e-6, 25.0, 5000)
slopes = maximum_slope(t, μ; ic=3, alg=LBFGS(), init=15.0, num_nodes=500)
x = LinRange(-10, 10, 2000)
f = y -> 1 / (1 + y^2)^(1 / 2)
f′ = y -> -y / (1 + y^2)^(3 / 2)
u = inviscid_solution(f, f′, x, t)
uₓ = similar(u)
for j in eachindex(t)
    for i in eachindex(x)
        uₓ[i, j] = f′(x[i] - u[i, j] * t[j]) / (1 + t[j] * f′(x[i] - u[i, j] * t[j]))
    end
end
auₓ = abs.(uₓ)
maxauₓ = zeros(length(t))
for i in eachindex(maxauₓ)
    @views maxauₓ[i] = maximum(auₓ[:, i])
    if i > 1 && maxauₓ[i] < maxauₓ[i-1] - 0.7
        maxauₓ[i:end] .= NaN
        break
    end
end
slopes[:, 1] .= maxauₓ
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
t_min = 1e-6
t_max = 25.0
fig = Figure(fontsize=33, resolution=(1700, 400))
ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"Maximum$ $ absolute slope", title=L"(a):$ $ Slope analysis", titlealign=:left,
    xticks=([0.0, 5.0, 10.0], [L"0", L"5", L"10"]), height=400, width=600,
    yticks=([0.0, 2.0, 4.0, 6.0, 8.0], [L"0", L"2", L"4", L"6", L"8"]))
colors = cgrad(LINSPECER_12_J, μ[μ_idx]; categorical=false)
lines!(ax, t, slopes[:, 1], color=colors[1], linestyle=:dash)
[lines!(ax, t, slopes[:, j], color=colors[i+1]) for (i, j) in enumerate(μ_idx[2:end])]
scatter!(ax, t[pt_idx], vals[pt_idx2], color=:red, markersize=4)
xlims!(ax, 0, t_max)

μ = μ_ENSTROPHY_IC3
_t = T_ENSTROPHY_IC3
E = ENSTROPHY_IC3

f = x -> 1 / (1 + x^2)^(1 / 2)
f′ = x -> -x / (1 + x^2)^(3 / 2)
u = (x, t) -> inviscid_solution(f, f′, x, t)
dudx = (x, t) -> ForwardDiff.derivative(x -> u(x, t), x)
dudx2 = (x, t) -> dudx(x, t)^2
nodes, weights = gausslegendre(600)
dudx2gh = (x, t) -> dudx(x, t)^2
Efnc = t -> gauss_legendre_finite(x -> dudx2gh(x, t), -5.0, 5.0, nodes, weights) / 2
Efnc0 = Efnc.(_t[_t.<2.5])
E[:, 1] .= NaN
E[_t.<2.5, 1] .= Efnc0

indices, vals = findmaxima(E)
_tvals = Vector{Float64}([])
_peaks = Vector{Float64}([])
_pt_idx = Vector{Int64}([])
_pt_idx2 = Vector{Int64}([])
_breakdown_μ = []
_μ_idx = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250]
for (i, (idx, val)) in enumerate(zip(indices, vals))
    if isnan(val)
        push!(_breakdown_μ, μ[i])
    else
        if (i ∈ μ_idx) && (i ≠ 1)
            push!(_pt_idx, idx)
            push!(_pt_idx2, i)
        end
        push!(_tvals, _t[idx])
        push!(_peaks, val)
    end
end
_breakdown_μ = _breakdown_μ[1] # Not so good
ax = Axis(fig[1, 2], width=600, height=400,
    xlabel=L"t", ylabel=L"E(t)",
    xticks=(0:5:25, [L"%$s" for s in 0:5:25]),
    yticks=(0:3, [L"%$s" for s in 0:3]),
    title=L"(b):$ $ Enstrophy", titlealign=:left)
μ_idx = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250]
ylims!(ax, 0, 2)
xlims!(ax, 0,t_max)
colors = cgrad(LINSPECER_12_J, μ[μ_idx]; categorical=false)
lines!(ax, _t, E[:, 1], color=colors[1], linestyle=:dash)
[lines!(ax, _t, E[:, j], color=colors[i+1]) for (i, j) in enumerate(μ_idx[2:end])]
scatter!(ax, _t[_pt_idx], vals[_pt_idx2], color=:red, markersize=4)
Colorbar(fig[1, 3], limits=(0.0, 0.5), colormap=colors, label=L"\mu", vertical=true,
    ticks=([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], [L"0.0", L"0.1", L"0.2", L"0.3", L"0.4", L"0.5"]))
resize_to_layout!(fig)
save("$FIGURES/slope_enstrophy_plot_1osqrt1px2.$EXTENSION", fig)

###########################################################################
## 1/(|x|^(1/2)(1+|x|))
###########################################################################
## Solution on the real line 
x = -12:0.1:12
t = [0.0, 0.5, 1.0, 5.0]
μ = [0.1, 1.0]
u_viscous = viscous_solution(x, t, μ; ic=4)

fig = Figure(fontsize=38, resolution=(800, 400))
ab = ['a', 'b']
for k in 1:length(μ)
    ax = Axis(fig[1, k], xlabel=L"x", ylabel=L"u(x,t)",
        title=L"(%$(ab[k])): $\mu = %$(μ[k])$", titlealign=:left,
        width=600, height=400,
        xticks=([-20, -10, 0, 10, 20], [L"-20", L"-10", L"0", L"10", L"20"]),
        yticks=([0, 0.25, 0.50, 0.75, 1.00], [L"0", L"0.25", L"0.5", L"0.75", L"1"]))
    for j in 1:length(t)
        lines!(ax, x, u_viscous[:, j, k], color=:black, linestyle=:solid, linewidth=2.3)
    end
end
resize_to_layout!(fig)

###########################################################################
## β = 1/4, 1/2, 1, 2
###########################################################################
x = [-20:0.1:20, -250:0.1:250.0]
t = [[0.0, 2.0, 4.0, 10.0], [0.0, 250.0, 500.0, 1000.0]]
μ = [0.1, 1.0]

β = [1 / 4, 1 / 2, 1, 2]
u_viscous = [zeros(length(β), length(x[i]), length(t[i]), length(μ)) for i in 1:2]
for j in 1:2
    for (i, b) in pairs(β)
        u_viscous[j][i, :, :, :] .= viscous_solution(x[j], t[j], μ; ic=(NaN, b))
    end
end

fig = Figure(fontsize=38, resolution=(2974.0586f0, 1095.8906f0))
ab = [["a", "b", "c", "d"], ["e", "f", "g", "h"]]
for r in 1:2
    for i in 1:length(β)
        ax = Axis(fig[r, i],
            xlabel=L"x", ylabel=L"u(x,t)",
            title=L"(%$(ab[r][i])): $\beta = %$(β[i])$", titlealign=:left,
            width=600, height=400,
            xticks=(r == 1 ? [-20, -10, 0, 10, 20] : -250:100:250, r == 1 ? [L"-20", L"-10", L"0", L"10", L"20"] : [L"%$s" for s in -250:100:250]),
            yticks=([0, 0.25, 0.50, 0.75, 1.00], [L"0", L"0.25", L"0.5", L"0.75", L"1"]))
        for k = 1:length(μ)
            for j in 1:length(t[r])
                lines!(ax, x[r], u_viscous[r][i, :, j, k], color=(:blue, :black, :red, :darkgreen)[j], linewidth=4, linestyle=(:solid, :dash)[k])
            end
        end
    end
end
resize_to_layout!(fig)
save("$FIGURES/all_viscous_plots_for_more_beta.$EXTENSION", fig)