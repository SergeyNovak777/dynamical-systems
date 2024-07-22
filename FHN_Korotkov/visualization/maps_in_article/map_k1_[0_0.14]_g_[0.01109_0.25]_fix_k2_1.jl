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

Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/k1_g_k2=1/LSE_300x300_k_1_g.jld2")["λs"]

length_range = 300;
k1range = range( 0.0, 0.14, length = length_range);
grange = range(0.01109, 0.25, length = length_range);

index = 1
absmax = maximum((Λs[:, :, index]))

mn, mx =  -absmax, absmax
GLMakie.activate!()  
f = Figure()    
ticksize = 35
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12

ax = Axis(f[1, 1], xlabel = L"g",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad,
            xtickformat = "{:.2f}", xticks = [0.05, 0.10, 0.15, 0.22])#,
            #xticks = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25], yticks = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14])

hm = heatmap!(ax, grange, k1range, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))

display(GLMakie.Screen(), f);

pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/k2=1_k1_g.pdf"
fullpath = pathtosave * filename 
save(fullpath, f)