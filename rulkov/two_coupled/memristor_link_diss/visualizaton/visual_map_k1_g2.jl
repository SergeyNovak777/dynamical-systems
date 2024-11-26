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

using JLD2, CairoMakie, GLMakie


Λs = load("/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/map_LSE/LSE_200x200_g_2_k_2.jld2")["λs"]
#load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_inh_from_colab/LSE_350x350_g_1_g_2.jld2")["λs"]
#load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma/LSE_350x350_g_1_g_2.jld2")["λs"]

length_map = 200;
range_p1 = range( 0.0, 10.0, length = length_map);
range_p2 = range( 0.0, 0.1, length = length_map);

index = 1
absmax = maximum((Λs[:, :, index]))

mn, mx =  -absmax, absmax
CairoMakie.activate!()  
f = Figure()    
ticksize = 35
labelsize = 45;
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12

ax = Axis(f[1, 1], xlabel = L"g_2",ylabel = L"k_2", xlabelsize = labelsize, ylabelsize = labelsize,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad, xtickformat = "{:.0f}", ytickformat = "{:.2f}",
            xticks = [0, 2, 4, 6, 8.5], yticks = [0, 0.02, 0.04, 0.06, 0.085])

hm = heatmap!(ax, range_p1, range_p2, Λs[:, :, index], colormap = :seismic,
                colorrange = (-0.05, 0.05))

display(GLMakie.Screen(), f);

pathtosave = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/images/"
filename = "/fix_g1=5_k1=0.1_map_k1_g2.pdf"
fullpath = pathtosave * filename 
save(fullpath, f)