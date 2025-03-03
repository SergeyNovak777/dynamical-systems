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


Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/LSE_400x400_k_1_k_2.jld2")["λs"]
u0s = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/u0s_400x400_k_1_k_2.jld2")
init_point = u0s["init_points"]
last_point = u0s["last_points"]

path_to_hopf_curve = "/home/sergey/MEGA/MatCont7p5/Systems/FHN/diagram/control_params:k1_k2_hopf_curve.mat"
file_hopf_curve = matopen(path_to_hopf_curve)
hopf_curve = read(file_hopf_curve, "x")
close(file_hopf_curve)

path_to_lc1_for = "/home/sergey/MEGA/MatCont7p5/Systems/FHN/diagram/limit_cycle_1_forward.mat"
file_lp1_for = matopen(path_to_lc1_for)
lp1_forward = read(file_lp1_for, "x") # index 645 k1, index 646 k2
close(file_lp1_for)

path_to_lc1_back = "/home/sergey/MEGA/MatCont7p5/Systems/FHN/diagram/limit_cycle_1_backward.mat"
file_lp1_back = matopen(path_to_lc1_back)
lp1_back = read(file_lp1_back, "x")
close(file_lp1_back)

length_range = 400;
k1range = range( 0.0, 0.12, length = length_range);
k2range = range(0.0, 80.0, length = length_range);

index = 1
absmax = maximum(abs.(Λs[:, :, index]))

mn, mx =  -absmax, absmax
CairoMakie.activate!()  

ticksize = 35
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12

f = Figure()
ax = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xticks = [0, 20, 40, 60, 78])
Λs[324:400, :, 1] .= -1
hm = heatmap!(ax, k2range, k1range, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))

# R1
color_r1 = :mediumorchid3;
lw_quad = 5.0;
lines!(ax, [0.0, 10.0], [0.0, 0.0], linewidth = lw_quad, color = color_r1)
lines!(ax, [0.0, 10.0], [0.1, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax, [0.0, 0.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax, [10.0, 10.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)
text!(ax, 11, 0.05, text = L"R_1", fontsize = ticksize*1.5)

# R3 
# k1 0.085, 0.094
# k2 60 80
lines!(ax, [60.0, 80.0], [0.085, 0.085], linewidth = lw_quad, color = color_r1)
lines!(ax, [60.0, 80.0], [0.094, 0.094], linewidth = lw_quad, color = color_r1)
lines!(ax, [60.0, 60.0], [0.085, 0.094], linewidth = lw_quad, color = color_r1)
lines!(ax, [80.0, 80.0], [0.085, 0.094], linewidth = lw_quad, color = color_r1)
text!(ax, 68, 0.062, text = L"R_3", fontsize = ticksize * 1.5)

# hopf curve
lines!(ax, hopf_curve[6, :], hopf_curve[5, :], linewidth = 3.0, color = :green3)
scatter!(ax, 39.2265, 0.0969, markersize = 12.0, color = :black)
text!(ax, 39.2265, 0.0969, text = L"GH", fontsize = 30, color = :black, align = (:center, :top))


# limit cycle 1
lines!(ax, lp1_forward[646, :], lp1_forward[645, :], linewidth = 1.5, color = :black)
lines!(ax, lp1_back[646, :], lp1_back[645, :], linewidth = 1.5, color = :black)
scatter!(ax, 38.8684, 0.0026, markersize = 12.0, color = :black)
text!(ax, 38.8684, 0.0026, text = L"NS", fontsize = 30, color = :black, align = (:center, :bottom))

xlims!(ax, 0.0, 80.0);
ylims!(ax, 0.0, 0.1)

display(GLMakie.Screen(), f);

pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/k1_k2_extended_with_curve.eps"
fullpath = pathtosave * filename 
save(fullpath, f)

#= index_p1 = 324
index_p2 = 324
println("k1: $(k1range[index_p1]); k1 index: $index_p1")
println("k2: $(k2range[index_p2]); k2 index: $index_p2")
println("u0: $(init_point[index_p1, index_p2, :])")
println("λs: $(Λs[index_p1, index_p2, index]) ")
println("last pont: $(last_point[index_p1,index_p2,:])") =#