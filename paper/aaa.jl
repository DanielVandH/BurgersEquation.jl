# Below is one way to make these plots in Julia, but we just use the MATLAB results since they seem to be more accurate. See "constants.jl" and "burger_aaa.m".
#=
L = 30
T = 2.0
N = 1000
t = [0.1, 0.5, 1.0, 8sqrt(3) / 9]
μ = 0.1
X = -4:0.01:4
Y = 0:0.01:4
U = viscous_solution_aaa(N, L, T, t, μ, X, Y)
=#

time_labels = ["0.1", "0.5", "1", "t_s"]
plot_layout = [[(1, 1), (1, 2)], [(1, 3), (1, 4)], [(2, 1), (2, 2)], [(2, 3), (2, 4)]]
alphabets = [['a', 'b'], ['c', 'd'], ['e', 'f'], ['g', 'h']]
fig = Figure(fontsize=48, resolution=(3200.0f0, 1156.5531f0))
for j in 1:length(T_AAA)
    portrait!(fig, X_AAA, Y_AAA, U_AAA[j], plot_layout[j][1][1], plot_layout[j][1][2]; nist=NIST,
        title=L"(%$(alphabets[j][1])): $t = %$(time_labels[j])$", titlealign=:left,
        xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
        width = 600, height = 400,
        yticks=([0.0, 2.0, 4.0], [L"0", L"2", L"4"]))
    landscape!(fig, X_AAA, Y_AAA, U_AAA[j], plot_layout[j][2][1], plot_layout[j][2][2]; nist=NIST,
        title=L"(%$(alphabets[j][2])): $t = %$(time_labels[j])$", titlealign=:left, zlims=(0, 4),
        xticks=([-4.0, 0.0, 4.0], ["-4", "0", "4"]),
        yticks=([0.0, 2.0, 4.0], ["0", "2", "4"]),
        zticks=([0.0, 2, 4], ["0", "2", "4"]),
        width = 600, height = 400,
        zlabel=L"|r(z, t)|")
    lines!(fig.content[end], X_AAA, zeros(length(X_AAA)), U_AAA[j][:, 1], color=:black, linewidth=5)
end 
resize_to_layout!(fig)
save("$FIGURES/AAA_viscous_solutions.$EXTENSION", fig)
