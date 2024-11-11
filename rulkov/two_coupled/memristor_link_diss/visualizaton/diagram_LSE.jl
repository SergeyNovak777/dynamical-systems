username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using JLD2, CairoMakie, GLMakie


array_LSEs = load("/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/g1=5_change_g2_LSE_diagram_LSEs.jld2");
array_LSEs = array_LSEs["array_LSEs"];
array_EEs = load("/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/g1=5_change_g2_LSE_diagram_EEs.jld2");
array_EEs = array_EEs["array_EEs"];
array_u0s = load("/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/g1=5_change_g2_LSE_diagram_u0s.jld2");
array_u0s = array_u0s["array_u0s"]

range_g2 = range(0.0, 10.0, length = 2000);
window_height, window_width = 400, 1000;
xlabel = L"g_2";
ylabel_LSEs, ylabelEEs = L"LSE", L"EE_{count}"
label_size = 35;
tickssize = 25;
linewidth = 2.0;
color_LLE1 = :red;
color_LLE2 = :green;
color_LLE3 = :blue;
color_LLE4 = :black;
xticks = [0, 2, 5, 8, 10];
yticks_EEs = [0, 2000, 4000, 6000];


CairoMakie.activate!();
figure = Figure(figsize = (window_height, window_width));
ax_LSE = Axis(figure[1,1],
        xlabel = xlabel, ylabel = ylabel_LSEs,
        xlabelsize = label_size, ylabelsize = label_size,
        xticklabelsize = tickssize, yticklabelsize = tickssize,
        xticks = xticks
);
ax_EE = Axis(figure[2, 1],
        xlabel = xlabel, ylabel = ylabelEEs,
        xlabelsize = label_size, ylabelsize = label_size,
        xticklabelsize = tickssize, yticklabelsize = tickssize,
        xticks = xticks, yticks = yticks_EEs
);

lines!(ax_LSE, range_g2, array_LSEs[:, 1], linewidth = linewidth, color = color_LLE1);
lines!(ax_LSE, range_g2, array_LSEs[:, 2], linewidth = linewidth, color = color_LLE2);
lines!(ax_LSE, range_g2, array_LSEs[:, 3], linewidth = linewidth, color = color_LLE3);

lines!(ax_EE, range_g2, array_EEs, linewidth = linewidth, color = color_LLE1);
display(GLMakie.Screen(), figure);

path_to_save_image = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/images/";
name_image = "diagram_LSE.eps";
full_path = path_to_save_image * name_image;
save(full_path, figure);
