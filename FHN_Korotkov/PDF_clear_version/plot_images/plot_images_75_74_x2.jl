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
    include("/home/sergey/work/repo/dynamical-systems/system.jl")
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl")
end

using StaticArrays, Statistics, CairoMakie, GLMakie, JLD2

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

path_to_load = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/local_maxs/"
filename_data_local_max = "data_local_max_x2.jld2"
filename_data_local_min = "data_local_min_x2.jld2"

data_local_max = load(path_to_load*filename_data_local_max)["data_local_max"]
data_local_min = load(path_to_load*filename_data_local_min)["data_local_min"]


all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
mean_amplitudes = Statistics.mean(all_amplitudes)
println("mean amplitude: $mean_amplitudes")
peaks_spikes, t_peaks_spikes, amplitudes_above_mean = select_spikes(data_local_min[1], data_local_max, mean_amplitudes) # Statistics.std(all_amplitudes))

data_local_max = nothing; data_local_min = nothing; GC.gc();

println("count spikes: $(length(peaks_spikes))")
Hs_x = Hs(amplitudes_above_mean, 7.0)

index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
peaks_EEs = peaks_spikes[index_EEs]
t_EEs = t_peaks_spikes[index_EEs]
println("count EEs: $(length(peaks_EEs))")


array_IEI = get_IEI(t_EEs)
#array_IEI = sort(array_IEI)
array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 10)

Hs_IEI_coeff_8 = Hs(array_IEI, 7)
Hs_IEI_coeff_6 = Hs(array_IEI, 5)

labelsize = 40;
ticksize = 30;

CairoMakie.activate!(type="png")

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/EEs/"
filename_hist = "PDF_IEE_hist_x2.eps"
filename_tEE_IEI = "t_EE_IEI_x2.eps"
filename_peaks_EE_IEI = "peaks_EE_IEI_x2.eps"

f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI_{x2}", ylabel = L"PDF_{IEI_{x2}}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_IEI, weights = array_PDF_IEI, bins = 100)
vlines!(ax, Hs_IEI_coeff_8, linewidth = 3.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_IEI_coeff_6, linewidth = 3.0, linestyle = :dash, color = :green)
display(f)
save(path_to_save*filename_hist, f)


f = Figure()
ax = Axis(f[1, 1], xlabel = L"t_{EE_{x2}}", ylabel = L"IEI_{x2}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_EEs[2:end], array_IEI, linewidth = 1.0)
hlines!(ax, Hs_IEI_coeff_8, linewidth = 3.0, linestyle = :dash, color = :red)
hlines!(ax, Hs_IEI_coeff_6, linewidth = 3.0, linestyle = :dash, color = :green)
display(f)
save(path_to_save*filename_tEE_IEI, f)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"peaks_{EE_{x2}}", ylabel = L"IEI_{x2}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
scatter!(ax, peaks_EEs[2:end], array_IEI, markersize = 1.0)
hlines!(ax, Hs_IEI_coeff_8, linewidth = 3.0, linestyle = :dash, color = :red)
hlines!(ax, Hs_IEI_coeff_6, linewidth = 3.0, linestyle = :dash, color = :green)
display(f)
save(path_to_save*filename_peaks_EE_IEI, f)