L = 30
T = 2.0
N = 1000
t = [0.1, 0.5, 1.0, 8sqrt(3) / 9]
μ = 0.1
X = -4:0.01:4
Y = -4:0.01:4
U = viscous_solution_aaa(N, L, T, t, μ, X, Y)
t = 1.0
μ = 0.1
u = viscous_solution(X, Y, t, μ)

fig = Figure(fontsize=25)
portrait!(fig, X, Y, U[:, :, 3], 1, 1; nist=NIST, width=400, height=400,
    title=L"(a):$ $ AAA approximant", titlealign=:left,
    xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
    yticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]))
portrait!(fig, X, Y, u, 1, 2; nist=NIST,
    title=L"(b):$ $ Exact solution", titlealign=:left, width=400, height=400,
    xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
    yticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]))

text!(fig.content[1], L"r_1", position=(1.0, 0.4), color=:white, textsize=36)
text!(fig.content[1], L"r_2", position=(-2.0, 2.0), color=:black, textsize=36)
text!(fig.content[1], L"r_3", position=(-2.0, -2.0), color=:white, textsize=36)
text!(fig.content[1], L"r_4", position=(1.0, -0.6), color=:white, textsize=36)
text!(fig.content[2], L"z_1", position=(1.0, 0.4), color=:white, textsize=36)
text!(fig.content[2], L"z_2", position=(-2.0, 2.0), color=:black, textsize=36)
text!(fig.content[2], L"z_3", position=(-2.0, -2.0), color=:white, textsize=36)
text!(fig.content[2], L"z_4", position=(1.0, -0.6), color=:white, textsize=36)

resize_to_layout!(fig)

save("$FIGURES/AAA_viscous_solutions_close1.$EXTENSION", fig)

#=
ua, xa, ta = viscous_solution_finite_diff(N, L, T, t, μ)
r = aaa(xa, ua[:, 1])
pr, rr, zr = prz(r)
z₁ = 1.5 + 0.8im
z₂ = -1.5 + 2.4im
z₃ = -1.5 - 2.4im
z₄ = 1.5 - 0.8im
r₁ = pr[argmin(abs.(pr .- z₁))]
r₂ = pr[argmin(abs.(pr .- z₂))]
r₃ = pr[argmin(abs.(pr .- z₃))]
r₄ = pr[argmin(abs.(pr .- z₄))]
num_nodes = 250
nodes, weights = gausshermite(num_nodes)
glnodes, glweights = gausslegendre(num_nodes)
z₁ = locate_pole(z₁, t, μ, nodes, weights, glnodes, glweights; newton=false)
z₂ = locate_pole(z₂, t, μ, nodes, weights, glnodes, glweights; newton=false)
z₃ = locate_pole(z₃, t, μ, nodes, weights, glnodes, glweights; newton=false)
z₄ = locate_pole(z₄, t, μ, nodes, weights, glnodes, glweights; newton=false)
[r₁ z₁
    r₂ z₂
    r₃ z₃
    r₄ z₄]
[100abs(r₁ - z₁) / abs(z₁),
    100abs(r₂ - z₂) / abs(z₂),
    100abs(r₃ - z₃) / abs(z₃),
    100abs(r₄ - z₄) / abs(z₄)]
=#