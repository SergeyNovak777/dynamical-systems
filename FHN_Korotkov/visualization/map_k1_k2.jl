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

using JLD2, CairoMakie, MAT


Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/LSE_400x400_k_1_k_2.jld2")["λs"]

length_range = 400;
k1range = range( 0.0, 0.12, length = length_range);
k2range = range(0.0, 80.0, length = length_range);

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

ax = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad)

hm = heatmap!(ax, k2range, k1range, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))

display(f);

#=pathtosave = "/home/sergey/MEGA/dynamical-systems/Tsodyks Markram/images/RJND/"
filename = "/hom_center.pdf"
fullpath = pathtosave * filename 
save(fullpath, f)=#