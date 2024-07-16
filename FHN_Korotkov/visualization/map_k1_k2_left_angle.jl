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


Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/map_LSE_k1_k2_to_2_g_0_1 debil/LSE_350x350_k_1_k_2.jld2")["λs"]
u0s = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/map_LSE_k1_k2_to_2_g_0_1 debil/u0s_350x350_k_1_k_2.jld2")
init_point = u0s["init_points"]
last_point = u0s["last_points"]

length_range = 350;
k1range = range( 0.0, 0.01, length = length_range);
k2range = range(0.0, 2.0, length = length_range);

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

Λs[340:350, :, 1] .= -1

ax = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad)

hm = heatmap!(ax, k2range, k1range, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))

display(GLMakie.Screen(), f);

#= pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/k1_k2_left_angle.pdf"
fullpath = pathtosave * filename 
save(fullpath, f) =#

findall(x->x==-1,Λs[:, :, index] )

index_p1 = 339
index_p2 = 134
println("k1: $(k1range[index_p1])")
println("k2: $(k2range[index_p2])")
println("u0: $(init_point[index_p1, index_p2, :])")
println("λs: $(Λs[index_p1, index_p2, index]) ")
println("last pont: $(last_point[index_p1,index_p2,:])")