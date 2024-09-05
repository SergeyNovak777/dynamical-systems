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


Λs = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/LSE_400x400_g_1_g_2.jld2")["λs"]
#load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_inh_from_colab/LSE_350x350_g_1_g_2.jld2")["λs"]
#load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma/LSE_350x350_g_1_g_2.jld2")["λs"]

length_range = 400;
g_1_range = range( 0.0, 10.0, length = length_range);
g_2_range = range( 0.0, 10.0, length = length_range);

index = 1
absmax = maximum(abs.(Λs[:, :, index]))

mn, mx =  -absmax, absmax
CairoMakie.activate!()  
f = Figure()    
ticksize = 35
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12

vector_ticks = [0, 2, 4, 6, 8, 10]

ax = Axis(f[1, 1], xlabel = L"g_2",ylabel = L"g_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xticks = vector_ticks, yticks = vector_ticks)

hm = heatmap!(ax, g_1_range, g_2_range, Λs[:, :, index], colormap = :seismic,
            colorrange = (-0.08, 0.08))
Colorbar(f[:, end+1], hm)
display(GLMakie.Screen(), f);

#= pathtosave = "/home/sergey/MEGA/dynamical-systems/Rulkov/Images/maps/"
filename = "/g=LSE_without_gamma.pdf"
fullpath = pathtosave * filename 
save(fullpath, f) =#