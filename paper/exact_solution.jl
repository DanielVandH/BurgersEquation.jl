x = LinRange(-5, 5, 250)
y = LinRange(0, 5, 250)
t = [0.1, 0.5, 1.0, tₛ]
μ = [0.1, 1.0]

# Compute solution 
u = viscous_solution(x, y, t, μ)
u_real = viscous_solution(x, t, μ)

# Plot the solutions at μ = 0.1
time_labels = ["0.1", "0.5", "1", "t_s"]
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
        landscape!(fig, x, y, u[:, :, j, k], plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
            title=L"(%$(alphabets[j][2])): $t = %$(time_labels[j])$", titlealign=:left, zlims=(0, 3),
            xticks=([-4.0, 0.0, 4.0], [L"-4", L"0", L"4"]),
            yticks=([-4.0, 0.0, 4.0], [L"-4", L"0", L"4"]),
            zticks=([0.0, 1, 2, 3], [L"0", L"1", L"2", L"3"]),
            width=600, height=400)
        lines!(fig.content[end], x, zeros(length(y)), u_real[:, j, k], color=:black, linewidth=5)
    end
    resize_to_layout!(fig)
    save("$FIGURES/viscous_solutions_$k.$EXTENSION", fig)
end