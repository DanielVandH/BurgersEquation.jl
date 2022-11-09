####### 1/(1+x^2)^(1/2)
####### The file name 1osqrt1px2 is 1 over sqrt 1 plus x squared

## Exact solution plots for O(1) time
x = LinRange(-5, 5, 250)
y = LinRange(0, 5, 250)
_tₛ = 3.0sqrt(3) / 2.0
t = [0.1, 0.5, 1.0, _tₛ]
μ = [0.1, 1.0]

# Compute solution 
u = viscous_solution(x, y, t, μ; ic=3)
u_real = viscous_solution(x, t, μ; ic=3)

# Plot the solutions at μ = 0.1, 1.0
time_labels = ["0.1", "0.5", "1.0", "t_s"]
plot_layout = [[(1, 1), (1, 2)], [(1, 3), (1, 4)], [(2, 1), (2, 2)], [(2, 3), (2, 4)]]
alphabets = [['a', 'b'], ['c', 'd'], ['e', 'f'], ['g', 'h']]
for k = 1:length(μ)
    fig = Figure(fontsize=38, resolution=(2000, 800))
    for j in 1:length(t)
        portrait!(fig, x, y, u[:, :, j, k], plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
            title=L"(%$(alphabets[j][1])): $t = %$(time_labels[j])$", titlealign=:left,
            xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            yticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            width=600, height=400)
        scatter!(fig.content[end], [0.0], [1.0], color=:black, markersize=12)
        landscape!(fig, x, y, u[:, :, j, k], plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
            title=L"(%$(alphabets[j][2])): $t = %$(time_labels[j])$", titlealign=:left, zlims=(0, 3),
            xticks=([-4.0, 0.0, 4.0], [L"-4", L"0", L"4"]),
            yticks=([-4.0, 0.0, 4.0], [L"-4", L"0", L"4"]),
            zticks=([0.0, 1, 2, 3], [L"0", L"1", L"2", L"3"]),
            width=600, height=400)
        lines!(fig.content[end], x, zeros(length(y)), u_real[:, j, k], color=:black, linewidth=5)
    end
    resize_to_layout!(fig)
    fig
    save("$FIGURES/viscous_solutions_initial_condition_1osqrt1px2_$k.$EXTENSION", fig)
end

## Exact solution plots for small time
x = LinRange(-0.008, 0.008, 250)
y = LinRange(0.990, 1.010, 250)
μ1 = 2.0
μ2 = 1.0
μ3 = 0.5
μ4 = 0.1
μ = [μ1, μ2, μ3, μ4]
t = 1e-6

# Compute solution 
u = viscous_solution(x, y, [t], μ; ic=3)

# Plot the solutions at each μ
fig = Figure(fontsize=38, resolution=(2449.9102f0, 604.33514f0))
plot_layout = [[(1, 1), (2, 1), (3, 1), (4, 1)], [(1, 2), (2, 2), (3, 2), (4, 2)],
    [(1, 3), (2, 3), (3, 3), (4, 3)], [(1, 4), (2, 4), (3, 4), (4, 4)]]
reverse!(plot_layout)
for j = 1:4
    portrait!(fig, Float64.(x), Float64.(y), Complex{Float64}.(u[:, :, 1, j]), plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(ALPHABET[4-j+1])): $\mu = %$(Float64(μ[j]))$, $t = 10^{-6}$", titlealign=:left,
        xlabel=L"x", ylabel=L"y", width=450, height=450,
        xticks=([-0.005, 0.0, 0.005], [L"-0.005", L"0", L"0.005"]),
        yticks=([0.995, 1.0, 1.005], [L"0.995", L"1", L"1.005"]))
    xlims!(fig.content[end], -0.006, 0.006)
    ylims!(fig.content[end], 0.994, 1.006)
end
save("$FIGURES/similarity_solutions_initial_condition_1osqrt1px2_.$EXTENSION", fig) # These aren't similarity solutions, but we'll just keep the prefix to make it easier to compare figures later with the other small time results

## Exact solution plots for large time 
x = LinRange(-200.0, 200.0, 250)
y = LinRange(0, 200.0, 250)
t = [100.0, 250.0, 500.0, 1000.0]
μ = [1.0]
u_vals = viscous_solution(x, y, t, μ; ic=3)

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

fig

save("$FIGURES/small_to_large_transition_initial_condition_1osqrt1px2.$EXTENSION", fig)
