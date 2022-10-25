using CairoMakie
using DelimitedFiles

x₀ = vec(readdlm("seminar/data/xt0.dat"))
xₛ = vec(readdlm("seminar/data/xts.dat"))
x₅ = vec(readdlm("seminar/data/xt5.dat"))
x₁₀ = vec(readdlm("seminar/data/xt10.dat"))
x₅_shock = vec(readdlm("seminar/data/xt5_shock.dat"))
x₁₀_shock = vec(readdlm("seminar/data/xt10_shock.dat"))
u₀ = vec(readdlm("seminar/data/ut0.dat"))
uₛ = vec(readdlm("seminar/data/uts.dat"))
u₅ = vec(readdlm("seminar/data/ut5.dat"))
u₁₀ = vec(readdlm("seminar/data/ut10.dat"))
u₅_shock = vec(readdlm("seminar/data/ut5_shock.dat"))
u₁₀_shock = vec(readdlm("seminar/data/ut10_shock.dat"))

fig = Figure(fontsize=38)
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"u(x, t)",
    width=600, height=400,
    xticks=([-20, -10, 0, 10, 20], [L"-20", L"-10", L"0", L"10", L"20"]),
    yticks=([0, 0.25, 0.50, 0.75, 1.00], [L"0", L"0.25", L"0.5", L"0.75", L"1"]))
lines!(ax, x₀, u₀, color=:black)
lines!(ax, xₛ, uₛ, color=:black)
lines!(ax, x₅, u₅, color=:black)
lines!(ax, x₁₀, u₁₀, color=:black)
lines!(ax, x₅_shock, u₅_shock, color=:red, linestyle=:dash)
lines!(ax, x₁₀_shock, u₁₀_shock, color=:red, linestyle=:dash)

save("seminar/figures/inviscid_solutions.pdf", fig)
save("seminar/figures/inviscid_solutions.png", fig)