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


Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/change_k2[0_80]_g_[0.095_0.25]_fix_k1/LSE_350x350_g_k_2.jld2")["λs"]
u0s = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/change_k2[0_80]_g_[0.095_0.25]_fix_k1/u0s_350x350_g_k_2.jld2")
init_point = u0s["init_points"]
last_point = u0s["last_points"]


length_range = 350;
grange = range( 0.095, 0.25, length = length_range);
k2range = range(0.0, 80.0, length = length_range);


index = 1
absmax = maximum((Λs[:, :, index]))

mn, mx =  -absmax, absmax
CairoMakie.activate!()  

ticksize = 35
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12

f = Figure()
ax = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"g", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad)
#Λs[324:400, :, 1] .= -1
hm = heatmap!(ax, k2range, grange, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))
#ylims!(0.0, 0.097)
display(GLMakie.Screen(), f);

#= pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/map_k2_g_fix_k1.pdf"
fullpath = pathtosave * filename 
save(fullpath, f) =#

#= index_p1 = 324
index_p2 = 324
println("k1: $(grange[index_p1]); g index: $index_p1")
println("k2: $(k2range[index_p2]); k2 index: $index_p2")
println("u0: $(init_point[index_p1, index_p2, :])")
println("λs: $(Λs[index_p1, index_p2, index]) ")
println("last pont: $(last_point[index_p1,index_p2,:])") =#