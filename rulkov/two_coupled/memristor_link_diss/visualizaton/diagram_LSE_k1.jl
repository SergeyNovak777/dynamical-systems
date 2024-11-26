username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using JLD2, CairoMakie, GLMakie


array_LSEs = load("/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/g1=5.0_g2=4.0_change_k1_LSE_diagram_LSEs.jld2");
array_LSEs = array_LSEs["array_LSEs"];
array_EEs = load("/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/g1=5.0_g2=4.0_change_k1_LSE_diagram_EEs.jld2");
array_EEs = array_EEs["array_EEs"];


range_g2 = range(0.0, 0.2, length = 3000);
window_height, window_width = 400, 1000;
xlabel = L"k_1";
ylabel_LSEs, ylabelEEs = L"λ_1", L"EE_{count}"
label_size = 30;
tickssize = 20;
linewidth = 2.0;
color_LLE1 = :red;
color_LLE2 = :green;
color_LLE3 = :blue;
color_LLE4 = :black;


CairoMakie.activate!();
figure = Figure(figsize = (window_height, window_width));
ax_LSE = Axis(figure[1,1],
        xlabel = xlabel, ylabel = ylabel_LSEs,
        xlabelsize = label_size, ylabelsize = label_size,
        xticklabelsize = tickssize, yticklabelsize = tickssize,
        xgridvisible = false, ygridvisible = false);
ax_EE = Axis(figure[2, 1],
        xlabel = xlabel, ylabel = ylabelEEs,
        xlabelsize = label_size, ylabelsize = label_size,
        xticklabelsize = tickssize, yticklabelsize = tickssize,
        xgridvisible = false, ygridvisible = false);

lines!(ax_LSE, range_g2, array_LSEs[:, 1], linewidth = linewidth, color = color_LLE1);

lines!(ax_EE, range_g2, array_EEs, linewidth = linewidth, color = :black);
display(GLMakie.Screen(), figure);

path_to_save_image = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/images/";
name_image = "g1=5_g2=4_diagram_LSE.eps";
full_path = path_to_save_image * name_image;
save(full_path, figure);
