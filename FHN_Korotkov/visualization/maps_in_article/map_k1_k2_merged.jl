if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
end

using JLD2, CairoMakie, MAT, GLMakie


Λs_k2_0_80 = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/LSE_400x400_k_1_k_2.jld2")["λs"]
length_range_k2_0_80 = 400;
k1range_k2_0_80 = range( 0.0, 0.12, length = length_range_k2_0_80);
k2range_k2_0_80 = range(0.0, 80.0, length = length_range_k2_0_80);

Λs_k2_0_10 = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/k1_k2_default_scale/LSE_300x300_k_1_k_2.jld2")["λs"]
length_range_k2_0_10 = 300;
k1range_k2_0_10 = range( 0.0, 0.1, length = length_range_k2_0_10);
k2range_k2_0_10 = range(0.0, 10.0, length = length_range_k2_0_10);

Λs_k2_0_2 = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/map_LSE_k1_k2_to_2_g_0_1 ___/LSE_350x350_k_1_k_2.jld2")["λs"]
length_range_k2_0_2 = 350;
k1range_k2_0_2 = range( 0.0, 0.1, length = length_range_k2_0_2);
k2range_k2_0_2 = range(0.0, 2.0, length = length_range_k2_0_2);

Λs_k2_0_9_1_2 = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/zoom_chaos_k1_k2_g=0.1/LSE_350x350_k_1_k_2.jld2")["λs"]
length_range_k2_0_9_1_2 = 350;
k1range_k2_0_9_1_2 = range( 0.010, 0.016, length = length_range_k2_0_9_1_2);
k2range_k2_0_9_1_2 = range(0.9, 1.2, length = length_range_k2_0_9_1_2);


index = 1
absmax_k2_0_80 = maximum(abs.(Λs_k2_0_80[:, :, index]))
mn_k2_0_80, mx_k2_0_80 =  -absmax_k2_0_80, absmax_k2_0_80

absmax_k2_0_10 = maximum(abs.(Λs_k2_0_10[:, :, index]))
mx_k2_0_10 = maximum(Λs_k2_0_10[:, :, index])
mn_k2_0_10 = -mx_k2_0_10

mx_k2_0_2 = maximum(Λs_k2_0_2[:, :, index])
mn_k2_0_2 = -mx_k2_0_2

absmax_k2_0_9_1_2 = maximum(abs.(Λs_k2_0_9_1_2[:, :, index]))
mn_k2_0_9_1_2, mx_k2_0_9_1_2 =  -absmax_k2_0_9_1_2, absmax_k2_0_9_1_2

CairoMakie.activate!()

ticksize = 35
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12
color_r1 = :mediumorchid3;
lw_quad = 5.0;
window_width =  350;
window_height = 300;
f = Figure(size=(1150,1150))

ax_k2_0_80 = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xticks = [0, 20, 40, 60, 78],
            height = window_height, width = window_width)

ax_k2_0_10 = Axis(f[1, 2], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xticks = [0, 4 , 8],
            height = window_height, width = window_width)

ax_k2_0_2 = Axis(f[2, 2], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xticks = [0.0, 0.5, 1.0, 1.5, 1.9],
            height = window_height, width = window_width)

ax_k2_0_9_1_2 = Axis(f[2, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xtickformat = "{:.2f}",
            xticks = [0.9, 1.0, 1.1, 1.17],
            yticks = [0.010, 0.012, 0.014, 0.016],
            height = window_height, width = window_width)

Λs_k2_0_80[324:400, :, 1] .= -1
Λs_k2_0_10[291:300, :, 1] .= -1
Λs_k2_0_2[340:350, :, 1] .= -1

hm = heatmap!(ax_k2_0_80, k2range_k2_0_80, k1range_k2_0_80, transpose(Λs_k2_0_80[:, :, index]), colormap = :seismic,
                colorrange = (mn_k2_0_80, mx_k2_0_80))
# R1
lines!(ax_k2_0_80, [0.0, 10.0], [0.0, 0.0], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_80, [0.0, 10.0], [0.1, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_80, [0.0, 0.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_80, [10.0, 10.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)
text!(ax_k2_0_80, 11, 0.05, text = L"R_1", fontsize = ticksize*1.5)
# R3 
# k1 0.085, 0.094
# k2 60 80
lines!(ax_k2_0_80, [60.0, 80.0], [0.085, 0.085], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_80, [60.0, 80.0], [0.094, 0.094], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_80, [60.0, 60.0], [0.085, 0.094], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_80, [80.0, 80.0], [0.085, 0.094], linewidth = lw_quad, color = color_r1)
text!(ax_k2_0_80, 67, 0.052, text = L"R_3", fontsize = ticksize * 1.5)

hm = heatmap!(ax_k2_0_10, k2range_k2_0_10, k1range_k2_0_10, transpose(Λs_k2_0_10[:, :, index]), colormap = :seismic,
                colorrange = (mn_k2_0_10, mx_k2_0_10))
lines!(ax_k2_0_10, [0.0, 2.0], [0.0, 0.0], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_10, [0.0, 2.0], [0.1, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_10, [0.0, 0.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_10, [2.0, 2.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)

hm = heatmap!(ax_k2_0_2, k2range_k2_0_2, k1range_k2_0_2, transpose(Λs_k2_0_2[:, :, index]), colormap = :seismic,
                colorrange = (mn_k2_0_2, mx_k2_0_2))

color_r1 = :mediumorchid3;
lw_quad = 5.0;
lines!(ax_k2_0_2, [0.0, 1.0], [0.090, 0.090], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_2, [0.0, 1.0], [0.095, 0.095], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_2, [1.0, 1.0], [0.090, 0.095], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_2, [0.0, 0.0], [0.090, 0.095], linewidth = lw_quad, color = color_r1)
text!(ax_k2_0_2, 0.5, 0.06, text = L"R_2", fontsize = ticksize * 1.5)

#=
k1range = range( 0.010, 0.016, length = length_range);
k2range = range(0.9, 1.2, length = length_range);
=#
lines!(ax_k2_0_2, [0.9, 1.2], [0.010, 0.010], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_2, [0.9, 1.2], [0.016, 0.016], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_2, [0.9, 0.9], [0.010, 0.016], linewidth = lw_quad, color = color_r1)
lines!(ax_k2_0_2, [1.2, 1.2], [0.010, 0.016], linewidth = lw_quad, color = color_r1)

hm = heatmap!(ax_k2_0_9_1_2, k2range_k2_0_9_1_2, k1range_k2_0_9_1_2, transpose(Λs_k2_0_9_1_2[:, :, index]), colormap = :seismic,
                colorrange = (mn_k2_0_9_1_2, mx_k2_0_9_1_2))

rowgap!(f.layout, 100)
colgap!(f.layout, 100)

display(GLMakie.Screen(), f);
display(f);
#= pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/k1_k2_extended.eps"
fullpath = pathtosave * filename 
save(fullpath, f) =#