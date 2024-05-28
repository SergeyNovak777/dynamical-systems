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

path_to_folder = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/"
filename_data_local_max = "data_local_max.jld2"
filename_data_local_min = "data_local_min.jld2"

data_local_max = load(path_to_folder*filename_data_local_max)["data_local_max"]
data_local_min = load(path_to_folder*filename_data_local_min)["data_local_min"]


all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
mean_amplitudes = Statistics.mean(all_amplitudes)

peaks_spikes, t_peaks_spikes, amplitudes_above_mean = select_spikes(data_local_min[1], data_local_max, mean_amplitudes) # Statistics.std(all_amplitudes))
println("count spikes: $(length(peaks_spikes))")
Hs_x = Hs(amplitudes_above_mean, 8)

index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
peaks_EEs = peaks_spikes[index_EEs]
t_EEs = t_peaks_spikes[index_EEs]
println("count EEs: $(length(peaks_EEs))")

array_IEI = get_IEI(t_EEs)
#array_IEI = sort(array_IEI)
array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 10)

Hs_IEI = Hs(array_IEI, 8)

labelsize = 40;
ticksize = 30;
path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/t_EE_IEI.eps"
CairoMakie.activate!()
f = Figure()
ax = Axis(f[1, 1], xlabel = L"t_{EE}", ylabel = L"IEI",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_EEs[2:end], array_IEI, linewidth = 1.0)
hlines!(ax, Hs_IEI, linewidth = 3.0, linestyle = :dash, color = :red)
save(path_to_save, f)
display(f)
