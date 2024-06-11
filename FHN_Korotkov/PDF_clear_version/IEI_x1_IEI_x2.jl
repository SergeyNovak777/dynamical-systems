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

path_to_folder = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/local_maxs/"

filename_data_local_max_x1 = "data_local_max_x1.jld2"
filename_data_local_min_x1 = "data_local_min_x1.jld2"
filename_data_local_max_x2 = "data_local_max_x2.jld2"
filename_data_local_min_x2 = "data_local_min_x2.jld2"


data_local_max_x1 = load(path_to_folder*filename_data_local_max_x1)["data_local_max"]
data_local_min_x1 = load(path_to_folder*filename_data_local_min_x1)["data_local_min"]
data_local_max_x2 = load(path_to_folder*filename_data_local_max_x2)["data_local_max"]
data_local_min_x2 = load(path_to_folder*filename_data_local_min_x2)["data_local_min"]

local_max_x1, t_local_max_x1 = data_local_max_x1
local_min_x1, t_local_min_x1 = data_local_min_x1

local_max_x2, t_local_max_x2 = data_local_max_x2
local_min_x2, t_local_min_x2 = data_local_min_x2

println("count local max 1: $(length(local_max_x1))")
println("count local min x1: $(length(local_min_x1))")
println("different between count local max and local min x1: $(abs( length(local_min_x1) - length(local_max_x1) ))")
println("-------------------------")
println("count local max x2: $(length(local_max_x2))")
println("count local min x2: $(length(local_min_x2))")
println("different between count local max and local min x2: $(abs( length(local_min_x2) - length(local_max_x2) ))")


t_plot_start = 1; t_plot_end = 10000;
f = Figure()
ax = Axis(f[1, 1])
scatter!(ax, t_local_max_x1[t_plot_start:t_plot_end], local_max_x1[t_plot_start:t_plot_end], markersize = 5.0, color = :red)
scatter!(ax, t_local_min_x1[t_plot_start:t_plot_end], local_min_x1[t_plot_start:t_plot_end], markersize = 5.0, color = :blue)
display(GLMakie.Screen(), f)

#=

all_amplitudes_x1 = get_amplitudes_all_events(data_local_max_x1[1], data_local_min_x1[1])
mean_amplitudes_x1 = Statistics.mean(all_amplitudes_x1)
println("mean amplitude x1: $mean_amplitudes_x1")
peaks_spikes_x1, t_peaks_spikes_x1, amplitudes_above_mean_x1 = select_spikes(data_local_min_x1[1], data_local_max_x1, mean_amplitudes_x1)

all_amplitudes_x2 = get_amplitudes_all_events(data_local_max_x2[1], data_local_min_x2[1])
mean_amplitudes_x2 = Statistics.mean(all_amplitudes_x2)
println("mean amplitude x2: $mean_amplitudes_x2")
peaks_spikes_x2, t_peaks_spikes_x2, amplitudes_above_mean_x2 = select_spikes(data_local_min_x2[1], data_local_max_x2, mean_amplitudes_x2)

data_local_max_x1 = nothing; data_local_min_x1 = nothing;
data_local_max_x2 = nothing; data_local_min_x2 = nothing;
GC.gc();

Hs_x1 = Hs(amplitudes_above_mean_x1, 8);
Hs_x2 = Hs(amplitudes_above_mean_x2, 8);
println("Hs_x1: $Hs_x1 ;");
println("Hs_x2: $Hs_x2 ;");

peaks_EE_x1 = peaks_spikes_x1[peaks_spikes_x1 .>= Hs_x1];
peaks_EE_x2 = peaks_spikes_x2[peaks_spikes_x2 .>= Hs_x2];

peaks_EE_x1

#= 
all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
mean_amplitudes = Statistics.mean(all_amplitudes)
println("mean amplitude: $mean_amplitudes")
peaks_spikes, t_peaks_spikes, amplitudes_above_mean = select_spikes(data_local_min[1], data_local_max, mean_amplitudes) # Statistics.std(all_amplitudes))

data_local_max = nothing; data_local_min = nothing; GC.gc();

println("count spikes: $(length(peaks_spikes))")
Hs_x = Hs(amplitudes_above_mean, 8)

index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
peaks_EEs = peaks_spikes[index_EEs]
t_EEs = t_peaks_spikes[index_EEs]
println("count EEs: $(length(peaks_EEs))")


array_IEI = get_IEI(t_EEs)
#array_IEI = sort(array_IEI)
array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 10)

Hs_IEI_coeff_8 = Hs(array_IEI, 8)
Hs_IEI_coeff_6 = Hs(array_IEI, 6)

labelsize = 40;
ticksize = 30;

CairoMakie.activate!()

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/EEs/"
filename_hist = "PDF_IEE_hist_x1.eps"
filename_tEE_IEI = "t_EE_IEI_x1.eps"
filename_peaks_EE_IEI = "peaks_EE_IEI_x1.eps"

f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI_{x1}", ylabel = L"PDF_{IEI_{x1}}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_IEI, weights = array_PDF_IEI, bins = 100)
vlines!(ax, Hs_IEI_coeff_8, linewidth = 3.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_IEI_coeff_6, linewidth = 3.0, linestyle = :dash, color = :green)
display(f)
save(path_to_save*filename_hist, f)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"t_{EE_{x1}}", ylabel = L"IEI{x1}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_EEs[2:end], array_IEI, linewidth = 1.0)
hlines!(ax, Hs_IEI_coeff_8, linewidth = 3.0, linestyle = :dash, color = :red)
hlines!(ax, Hs_IEI_coeff_6, linewidth = 3.0, linestyle = :dash, color = :green)
display(f)
save(path_to_save*filename_tEE_IEI, f)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"peaks_{EE_{x1}}", ylabel = L"IEI_{x1}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
scatter!(ax, peaks_EEs[2:end], array_IEI, markersize = 1.0)
hlines!(ax, Hs_IEI_coeff_8, linewidth = 3.0, linestyle = :dash, color = :red)
hlines!(ax, Hs_IEI_coeff_6, linewidth = 3.0, linestyle = :dash, color = :green)
display(f)
save(path_to_save*filename_peaks_EE_IEI, f) =#