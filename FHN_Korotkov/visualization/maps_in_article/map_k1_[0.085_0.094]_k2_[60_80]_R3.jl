if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "irrito"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
end

using JLD2, CairoMakie, MAT, GLMakie


Λs = load("/home/irrito/MEGA downloads/dynamical-systems/FHN_Korotkov/data/maps_LSE/R3/LSE_350x350_k_1_k_2.jld2")["λs"]
u0s = load("/home/irrito/MEGA downloads/dynamical-systems/FHN_Korotkov/data/maps_LSE/R3/u0s_350x350_k_1_k_2.jld2")
init_point = u0s["init_points"]
last_point = u0s["last_points"]


length_range = 350;
k1range = range( 0.085, 0.094, length = length_range);
k2range = range(60.0, 80.0, length = length_range);


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
            xticklabelpad = tickpad, yticklabelpad = tickpad)#,
            #xticks = [63, 70, 79],
            #yticks = [0.086, 0.090, 0.093])

hm = heatmap!(ax, k2range, k1range, transpose(Λs[:, :, index]), colormap = :seismic,
                colorrange = (mn, mx))
Colorbar(f[:, end+1], hm)
display(GLMakie.Screen(), f);


index = 2
absmax = maximum(abs.(Λs[:, :, index]))

mn, mx = -0.01, 0.01 #-absmax, absmax
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
            xticklabelpad = tickpad, yticklabelpad = tickpad)#,
            #xticks = [63, 70, 79],
            #yticks = [0.086, 0.090, 0.093])

hm = heatmap!(ax, k2range, k1range, transpose(Λs[:, :, index]), colormap = :cyclic_tritanopic_wrwc_70_100_c20_n256,
                colorrange = (mn, mx))
Colorbar(f[:, end+1], hm)
display(GLMakie.Screen(), f);


#= pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/maps"
filename = "/map_R3.pdf"
fullpath = pathtosave * filename 
save(fullpath, f) =#

function check_condition(spectrum)
    checknull(spectrum[1]) && spectrum[2]<0  ? color = 5 :
    checknull(spectrum[1]) && checknull(spectrum[2])  && spectrum[3]<0 ? color = 4 :
    checknull(spectrum[1]) && checknull(spectrum[2])  && checknull(spectrum[3]) && spectrum[4]<0 ? color = 3 :
    spectrum[1]>0  && checknull(spectrum[2])  && spectrum[3]<0 ? color = 2 :
    spectrum[1]>0  && spectrum[2]>0   && checknull(spectrum[3]) && spectrum[4]<0 ? color = 1 :
    color = 6
end

function checknull(value)
    isapprox(value, 0.0; atol = 1e-2)
end

color_matrix00 = zeros(350, 350)

for (i, k1) in enumerate(k2range)
    
    for (j, k2) in enumerate(k1range)
        spectrum = Λs[i, j, :]
        color_matrix00[i, j] = check_condition(spectrum)
    end
end


f = Figure()
ax = Axis(f[1, 1], xlabel = L"k_2",ylabel = L"k_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad)#,
            #xticks = [63, 70, 79],
            #yticks = [0.086, 0.090, 0.093])

#hm = heatmap!(ax, k2range, k1range, transpose(color_matrix00), colormap = :RdBu_6)
hm = heatmap!(ax, k1range, k2range, color_matrix00, colormap = :RdBu_6)
Colorbar(f[:, end+1], hm)
display(GLMakie.Screen(), f);
