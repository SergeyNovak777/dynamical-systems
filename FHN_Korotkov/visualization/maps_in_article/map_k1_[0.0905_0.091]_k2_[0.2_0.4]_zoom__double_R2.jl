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


Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/strip_chaos_zoom_double/LSE_350x350_k_1_k_2.jld2")["λs"]

length_range = 350;
k1range = range( 0.0905, 0.091, length = length_range);
k2range = range(0.2, 0.4, length = length_range);

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
            xtickformat = "{:.1f}",)
            #xticks = [0.0, 0.5, 0.97])

hm = heatmap!(ax, k2range, k1range, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))

display(GLMakie.Screen(), f);

#= pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/map_R2.pdf"
fullpath = pathtosave * filename 
save(fullpath, f) =#