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


EEs = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/with_gamma_broke_first_iter/EEs/EEs_400x400_g_1_g_2.jld2")["matrix_EEs"]


length_range = 400;
g_1_range = range( 0.0, 10.0, length = length_range);
g_2_range = range( 0.0, 10.0, length = length_range);

absmax = maximum(abs.(EEs[:, :]))

mn, mx =  -100, 100
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

hm = heatmap!(ax, g_1_range, g_2_range, EEs, colormap = :seismic,
            colorrange = (mn, mx))

display(GLMakie.Screen(), f);

#= pathtosave = "/home/sergey/MEGA/dynamical-systems/Rulkov/Images/maps/"
filename = "/g=LSE_without_gamma.pdf"
fullpath = pathtosave * filename 
save(fullpath, f) =#