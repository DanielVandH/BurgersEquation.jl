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
t = T_AAA[3]
u = viscous_solution(X_AAA, Y_AAA, t, μ_AAA)

fig = Figure(fontsize=33)
portrait!(fig, X_AAA, Y_AAA, U_AAA_3, 1, 1; nist=NIST, width=600, height=400,
    title=L"(a):$ $ AAA approximant", titlealign=:left,
    xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
    yticks=([ 0.0, 2.0, 4.0], [L"0", L"2", L"4"]))
portrait!(fig, X_AAA, Y_AAA, u, 1, 2; nist=NIST,
    title=L"(b):$ $ Exact solution", titlealign=:left, width=600, height=400,
    xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
    yticks=([.0, 2.0, 4.0], [L"0", L"2", L"4"]))

text!(fig.content[1], L"r_0^{(1)}", position=(0.7, 0.4), color=:white, textsize=36)
text!(fig.content[1], L"r_0^{(2)}", position=(-2.3, 2.0), color=:black, textsize=36)
text!(fig.content[2], L"z_0^{(1)}", position=(1.0, 0.4), color=:white, textsize=36)
text!(fig.content[2], L"z_0^{(2)}", position=(-2.0, 1.9), color=:black, textsize=36)

resize_to_layout!(fig)

save("$FIGURES/AAA_viscous_solutions_close1.$EXTENSION", fig)