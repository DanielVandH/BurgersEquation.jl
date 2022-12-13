## Prepare
x = LinRange(-5, 5, 250)
y = LinRange(0, 5, 250)
f₁ = z -> 1 / (1 + z^2)
f₁′ = z -> -2z / (1 + z^2)^2
f₂ = z -> 1 / (1 + z^2)^2
f₂′ = z -> -4z / (1 + z^2)^3
f₃ = z -> 1 / (1 + z^2)^(1 / 2)
f₃′ = z -> -z / (1 + z^2)^(3 / 2)
tₛ₁ = 8sqrt(3) / 9
tₛ₂ = 54sqrt(5) / 125
tₛ₃ = 3sqrt(3) / 2
t1 = [0.1, 0.5, 1.0, tₛ₁]
t2 = [0.1, 0.5, 0.75, tₛ₂]
t3 = [0.1, 0.5, 1.0, tₛ₃]
t = [t1, t2, t3]

## Obtain the solutions 
u₁ = inviscid_solution(f₁, f₁′, x, y, t1; Δt=t1 / 1000)
u₂ = inviscid_solution(f₂, f₂′, x, y, t2; Δt=t1 / 1000)
u₃ = inviscid_solution(f₃, f₃′, x, y, t3; Δt=t1 / 1000)
u = [u₁, u₂, u₃]
u₁_real = inviscid_solution(f₁, x, t1)
u₂_real = inviscid_solution(f₂, x, t2)
u₃_real = inviscid_solution(f₃, x, t3)
u_real = [u₁_real, u₂_real, u₃_real]

## Figures 
time_labels = [
    ["0.1", "0.5", "1", "8/27^{1/2}"],
    ["0.1", "0.5", "0.75", "54/5^{5/2}"],
    ["0.1", "0.5", "1", "27^{1/2}/2"]
]
ic_titles = [
    "1/(1+x^2)",
    "1/(1+x^2)^2",
    "1/(1+x^2)^{1/2}"
]
plot_layout = [
    [(1, 1), (1, 2), (1, 3), (1, 4)],
    [(2, 1), (2, 2), (2, 3), (2, 4)],
    [(3, 1), (3, 2), (3, 3), (3, 4)]
]
alphabets = [
    ["a", "b", "c", "d"],
    ["e", "f", "g", "h"],
    ["i", "j", "k", "l"]
]

for k in eachindex(t)
    fig = Figure(fontsize=38, resolution=(2779.8281f0, 555.0031f0))
    for j in eachindex(t[k])
        portrait!(fig, x, y, u[k][:, :, j],
            plot_layout[k][j][1], plot_layout[k][j][2]; nist=NIST,
            title=L"(%$(alphabets[k][j])): $t = %$(time_labels[k][j])$, $u(x, 0) = %$(ic_titles[k])$", titlealign=:left,
            xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            yticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            width=600, height=400)
    end
    save("$FIGURES/inviscid_solutions_for_each_initial_condition_portrait_row_$k.$EXTENSION", fig)
end

for k in eachindex(t)
    fig = Figure(fontsize=38, resolution=(2744.0f0, 492.0f0))
    for j in eachindex(t[k])
        landscape!(fig, x, y, u[k][:, :, j],
            plot_layout[k][j][1], plot_layout[k][j][2]; nist=NIST,
            title=L"(%$(alphabets[k][j])): $t = %$(time_labels[k][j])$, $u(x, 0) = %$(ic_titles[k])$", titlealign=:left,
            xticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            yticks=([-4.0, -2.0, 0.0, 2.0, 4.0], [L"-4", L"-2", L"0", L"2", L"4"]),
            zlims=(0, 1),
            zticks=([0.0, 0.5, 1.0], [L"0", L"0.5", L"1.0"]),
            width=600, height=400)
        lines!(fig.content[end], x, zeros(length(y)), u_real[k][:, j], color=:black, linewidth=5)
    end
    save("$FIGURES/inviscid_solutions_for_each_initial_condition_landscape_row_$k.$EXTENSION", fig)
end

#### Real line 
## Obtain the solutions 
u₁_real = inviscid_solution(f₁, x, t1)
u₂_real = inviscid_solution(f₂, x, t2)
u₃_real = inviscid_solution(f₃, x, t3)
u_real = [u₁_real, u₂_real, u₃_real]

## Figures 
time_labels = [
    ["0.1", "0.5", "1", "8/27^{1/2}"],
    ["0.1", "0.5", "0.75", "54/5^{5/2}"],
    ["0.1", "0.5", "1", "27^{1/2}/2"]
]
ic_titles = [
    "1/(1+x^2)",
    "1/(1+x^2)^2",
    "1/(1+x^2)^{1/2}"
]
plot_layout = [
    [(1, 1), (1, 2), (1, 3), (1, 4)],
    [(2, 1), (2, 2), (2, 3), (2, 4)],
    [(3, 1), (3, 2), (3, 3), (3, 4)]
]
alphabets = [
    ["a", "b", "c", "d"],
    ["e", "f", "g", "h"],
    ["i", "j", "k", "l"]
]

fig = Figure(fontsize=38, resolution=(2779.8281f0, 555.0031f0))
for k in eachindex(t)
    for j in eachindex(t[k])
        ax = Axis(fig[plot_layout[k][j][1], plot_layout[k][j][2]])
        lines!(ax, x, u_real[k][:, j])
    end
    #save("$FIGURES/inviscid_solutions_for_each_initial_condition_portrait_row_$k.$EXTENSION", fig)
end
