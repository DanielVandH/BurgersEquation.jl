L = 30
T = 2.0
N = 1000
t = [0.1, 0.5, 1.0, 8sqrt(3) / 9]
μ = 0.1
X = -4:0.01:4
Y = 0:0.01:4
U = viscous_solution_aaa(N, L, T, t, μ, X, Y)

time_labels = ["0.1", "0.5", "1", "t_s"]
plot_layout = [[(1, 1), (1, 2)], [(1, 3), (1, 4)], [(2, 1), (2, 2)], [(2, 3), (2, 4)]]
alphabets = [['a', 'b'], ['c', 'd'], ['e', 'f'], ['g', 'h']]
fig = Figure(fontsize=38, resolution=(2000, 800))
for j in 1:length(t)
    portrait!(fig, X, Y, U[:, :, j], plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(alphabets[j][1])): $t = %$(time_labels[j])$", titlealign=:left,
        xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
        yticks=([0.0, 2.0, 4.0], [L"0", L"2", L"4"]))
    landscape!(fig, X, Y, U[:, :, j], plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
        title=L"(%$(alphabets[j][2])): $t = %$(time_labels[j])$", titlealign=:left, zlims=(0, 5),
        xticks=([-4.0, 0.0, 4.0], [L"-4", L"0", L"4"]),
        yticks=([0.0, 2.0, 4.0], [L"0", L"2", L"4"]),
        zticks=([0.0, 2.5, 5.0], [L"0", L"2.5", L"5"]),
        zlabel=L"|r(z, t)|")
    lines!(fig.content[end], X, zeros(length(X)), U[:, 1, j], color=:black, linewidth=5)
end
resize_to_layout!(fig)
save("$FIGURES/AAA_viscous_solutions.$EXTENSION", fig)
