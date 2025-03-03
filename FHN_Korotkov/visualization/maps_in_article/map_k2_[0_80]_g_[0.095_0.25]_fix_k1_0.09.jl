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


# curve of hopf
path_to_hopf_curve_for = "/home/sergey/MEGA/MatCont7p5/Systems/FHN/diagram/hopf_curve_control_params:k2_g_forward.mat"
file_hopf_curve_for = matopen(path_to_hopf_curve_for)
hopf_curve_for = read(file_hopf_curve_for, "x")
close(file_hopf_curve_for)

path_to_hopf_curve_back = "/home/sergey/MEGA/MatCont7p5/Systems/FHN/diagram/hopf_curve_control_params:k2_g_backward.mat"
file_hopf_curve_back = matopen(path_to_hopf_curve_back)
hopf_curve_back = read(file_hopf_curve_back, "x")
close(file_hopf_curve_back)

length_range = 350;
grange = range( 0.095, 0.25, length = length_range);
k2range = range(0.0, 80.0, length = length_range);


index = 1

mn, mx =  minimum((Λs[:, :, index])), maximum((Λs[:, :, index]))
CairoMakie.activate!()  

#= index_p1 = 262
index_p2 = 324
println("g: $(grange[index_p1]); g index: $index_p1")
println("k2: $(k2range[index_p2]); k2 index: $index_p2")
println("u0: $(init_point[index_p1, index_p2, :])")
println("λs: $(Λs[index_p1, index_p2, index]) ")
println("last pont: $(last_point[index_p1,index_p2,:])") =#

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

hm = heatmap!(ax, k2range, grange, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (-0.3, 0.3))

# hopf
lines!(ax, hopf_curve_for[6, :], hopf_curve_for[5, :], linewidth = 3.0, color = :magenta)
lines!(ax, hopf_curve_back[6, :], hopf_curve_back[5, :], linewidth = 3.0, color = :magenta)

# GH
scatter!(ax, 43.373919, 0.23276985, markersize = 12.0, color = :black)
text!(ax, 43.373919, 0.23276985, text = L"GH", fontsize = 30, color = :black, align = (:center, :top), offset = (0, -10))

xlims!(ax, 0.0, 80);
ylims!(ax, 0.095, 0.25);

display(GLMakie.Screen(), f);

pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/map_k2_g_fix_k1_extended_with_curve.pdf"
fullpath = pathtosave * filename 
#save(fullpath, f)



#findall(x->x== maximum((Λs[:, :, index])), Λs[:, :, index])

#findall(x-> x == mx, Λs[:, :, 1])