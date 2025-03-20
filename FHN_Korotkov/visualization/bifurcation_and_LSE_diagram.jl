username = "irrito"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")
include(pathtorepo * "/system.jl")

using JLD2, CairoMakie, GLMakie;

path_to_load = "/home/irrito/MEGA2/dynamical_systems/FHN/data/";
file_name_dia_Poincare = "dia_Poincare_k1=0.09_change_k2_from_0_to_100.jld2";
file_name_matrix_LSE = "matrix_LSE_fix_k1=0.09_change_k2_from_0_to_100.jld2";

poincare_diagram = load(path_to_load*file_name_dia_Poincare)["output"];
matrix_LSE = load(path_to_load*file_name_matrix_LSE)["matrix_LSE"];

k2_start = 0.0
k2_end = 100.0
len = 1000
rangek2 = range(k2_start, k2_end, length = len)

markersize = 1.25;
lbsize = 50;
ticksize = 28;

CairoMakie.activate!();

fig = Figure(size = (1200, 350))

axis_poincare = Axis(fig[1,1],
        ylabel = L"x_1",
        xlabelsize = lbsize, ylabelsize = lbsize,
        xticklabelsize = ticksize,yticklabelsize = ticksize,
        xgridvisible = true, ygridvisible = true, xticks = [0, 25, 50, 75, 100],xticklabelsvisible = false);

axis_LSE = Axis(fig[2,1],
        xlabel = L"k_2", ylabel = L"x_1",
        xlabelsize = lbsize, ylabelsize = lbsize,
        xticklabelsize = ticksize,yticklabelsize = ticksize,
        xgridvisible = true, ygridvisible = true, xticks = [0, 25, 50, 75, 100]);



for (j, p) in enumerate(rangek2)
scatter!(axis_poincare, fill(p, length(poincare_diagram[j])), poincare_diagram[j]; color = ("black", 0.25), markersize = markersize)
end

lines!(axis_LSE, rangek2, matrix_LSE[:, 1], linewidth = 1.5, color = :red);
lines!(axis_LSE, rangek2, matrix_LSE[:, 2], linewidth = 1.5, color = :blue);
hlines!(axis_LSE, 0, rangek2, linestyle = :dash, linewidth = 1.5, color = :black)

display(GLMakie.Screen(), fig)

path_to_save_image = "/home/irrito/MEGA2/dynamical_systems/FHN/images/";
file_name_image = "Poincare_LSE_diagram.eps";

save(path_to_save_image*file_name_image, fig)