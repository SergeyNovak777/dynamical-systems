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


Λs_R3_norm = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/R3/LSE_350x350_k_1_k_2.jld2")["λs"]
length_range_R3_norm = 350;
k1range_R3_norm = range( 0.085, 0.094, length = length_range_R3_norm);
k2range_R3_norm = range(60.0, 80.0, length = length_range_R3_norm);

Λs_R3_zoom = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/R3_zoom/LSE_350x350_k_1_k_2.jld2")["λs"]
length_range_R3_zoom = 350;
k1range_R3_zoom = range( 0.09075, 0.092, length = length_range_R3_zoom);
k2range_R3_zoom = range(62.5, 74.0, length = length_range_R3_zoom);

index = 1;

absmax_R3_norm = maximum(abs.(Λs_R3_norm[:, :, index]));
mn_R3_norm, mx_R3_norm =  -absmax_R3_norm, absmax_R3_norm;

absmax_R3_zoom = maximum(abs.(Λs_R3_zoom[:, :, index]));
mn_R3_zoom, mx_R3_zoom =  -absmax_R3_zoom, absmax_R3_zoom;

CairoMakie.activate!()  

ticksize = 30
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12

color_r1 = :black
lw_quad = 5.0;
window_width =  280;
window_height = 200;
f = Figure(size=(920,340))


ax_R3_norm = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            height = window_height, width = window_width,
            xticks = [60, 70, 80],
            yticks = [0.086, 0.090, 0.094])

hm = heatmap!(ax_R3_norm, k2range_R3_norm, k1range_R3_norm, transpose(Λs_R3_norm[:, :, index]), colormap = :seismic,
                colorrange = (mn_R3_norm, mx_R3_norm));

lines!(ax_R3_norm, [62.5, 74.0], [0.09075, 0.09075], linewidth = lw_quad, color = color_r1)
lines!(ax_R3_norm, [62.5, 74.0], [0.092, 0.092], linewidth = lw_quad, color = color_r1)
lines!(ax_R3_norm, [62.5, 62.5], [0.09075, 0.092], linewidth = lw_quad, color = color_r1)
lines!(ax_R3_norm, [74.0, 74.0], [0.09075, 0.092], linewidth = lw_quad, color = color_r1)

ax_R3_zoom = Axis(f[1, 2], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            height = window_height, width = window_width,
            xticks = [63, 68, 73])

hm = heatmap!(ax_R3_zoom, k2range_R3_zoom, k1range_R3_zoom, transpose(Λs_R3_zoom[:, :, index]), colormap = :seismic,
                colorrange = (mn_R3_zoom, mx_R3_zoom))

display(f);

pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/map_R3_merged.pdf"
fullpath = pathtosave * filename 
save(fullpath, f)