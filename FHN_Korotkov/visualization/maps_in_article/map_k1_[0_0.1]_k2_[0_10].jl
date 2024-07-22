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

using JLD2, CairoMakie, GLMakie, MAT


Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/k1_k2_default_scale/LSE_300x300_k_1_k_2.jld2")["λs"]
u0s = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/k1_k2_default_scale/u0s_300x300_k_1_k_2.jld2")
init_point = u0s["init_points"]
last_point = u0s["last_points"]

length_range = 300;
k1range = range( 0.0, 0.1, length = length_range);
k2range = range(0.0, 10.0, length = length_range);

index = 1
absmax = maximum(abs.(Λs[:, :, index]))

mx = maximum(Λs[:, :, index])
mn = -mx
CairoMakie.activate!()  
f = Figure()    
ticksize = 35
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12

ax = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xticks = [0, 4 , 8])
# 288
Λs[291:300, :, 1] .= -1
hm = heatmap!(ax, k2range, k1range, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))
color_r1 = :mediumorchid3;
lw_quad = 5.0;
lines!(ax, [0.0, 2.0], [0.0, 0.0], linewidth = lw_quad, color = color_r1)
lines!(ax, [0.0, 2.0], [0.1, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax, [0.0, 0.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)
lines!(ax, [2.0, 2.0], [0.0, 0.1], linewidth = lw_quad, color = color_r1)

display(GLMakie.Screen(), f);

#= pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/g=0_1_k1_k2.eps"
fullpath = pathtosave * filename 
save(fullpath, f) =#

#= index_p1 = 291
index_p2 = 134
println("k1: $(k1range[index_p1])")
println("k2: $(k2range[index_p2])")
println("u0: $(init_point[index_p1, index_p2, :])")
println("λs: $(Λs[index_p1, index_p2, index]) ")
println("last pont: $(last_point[index_p1,index_p2,:])") =#