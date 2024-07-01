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

using StaticArrays, Statistics, JLD2, CairoMakie, GLMakie

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

path_to_folder = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/local_maxs/"

filename_data_local_max_x1 = "data_local_max_x1.jld2"
filename_data_local_min_x1 = "data_local_min_x1.jld2"
filename_data_local_max_x2 = "data_local_max_x2.jld2"
filename_data_local_min_x2 = "data_local_min_x2.jld2"

local_max_x1, t_local_max_x1 = load(path_to_folder*filename_data_local_max_x1)["data_local_max"]
local_min_x1, t_local_min_x1 = load(path_to_folder*filename_data_local_min_x1)["data_local_min"]

println("count local max 1: $(length(local_max_x1))")
println("count local min x1: $(length(local_min_x1))")
println("different between count local max and local min x1: $(abs( length(local_min_x1) - length(local_max_x1) ))")
println("-------------------------")

all_amplitudes_x1 = get_amplitudes_all_events(local_max_x1, local_min_x1);
mean_amplitudes_x1 = Statistics.mean(all_amplitudes_x1);
println("mean amplitude x1: $mean_amplitudes_x1");
peaks_spikes_x1, t_peaks_spikes_x1, amplitudes_above_mean_x1 = select_spikes(local_min_x1, [local_max_x1, t_local_max_x1], mean_amplitudes_x1);

local_max_x1 = nothing; t_local_max_x1 = nothing;
local_min_x1 = nothing; t_local_min_x1 = nothing;
GC.gc();

local_max_x2, t_local_max_x2 = load(path_to_folder*filename_data_local_max_x2)["data_local_max"]
local_min_x2, t_local_min_x2  = load(path_to_folder*filename_data_local_min_x2)["data_local_min"]

println("count local max x2: $(length(local_max_x2))")
println("count local min x2: $(length(local_min_x2))")
println("different between count local max and local min x2: $(abs( length(local_min_x2) - length(local_max_x2) ))")

all_amplitudes_x2 = get_amplitudes_all_events(local_max_x2, local_min_x2)
mean_amplitudes_x2 = Statistics.mean(all_amplitudes_x2)
println("mean amplitude x2: $mean_amplitudes_x2")
peaks_spikes_x2, t_peaks_spikes_x2, amplitudes_above_mean_x2 = select_spikes(local_min_x2, [local_max_x2, t_local_max_x2], mean_amplitudes_x2)

local_max_x2 = nothing; t_local_max_x2 = nothing;
local_min_x2 = nothing; t_local_min_x2 = nothing;
GC.gc();

Hs_x1 = Hs(amplitudes_above_mean_x1, 8);
Hs_x2 = Hs(amplitudes_above_mean_x2, 8);
println("Hs_x1: $Hs_x1 ;");
println("Hs_x2: $Hs_x2 ;");

index_EEs_x1 = findall(x-> x >= Hs_x1, peaks_spikes_x1)
peaks_EEs_x1 = peaks_spikes_x1[index_EEs_x1]
t_EEs_x1 = t_peaks_spikes_x1[index_EEs_x1]

index_EEs_x2 = findall(x-> x >= Hs_x2, peaks_spikes_x2)
peaks_EEs_x2 = peaks_spikes_x2[index_EEs_x2]
t_EEs_x2 = t_peaks_spikes_x2[index_EEs_x2]

labelsize = 40;
ticksize = 30;

f = Figure();
ax1 = Axis(f[1, 1]);
lines!(ax1, amplitudes_above_mean_x1[1:1200_000], linewidth = 0.5, color = :black);
display(GLMakie.Screen(), f)

f = Figure();
ax2 = Axis(f[1, 1]);
lines!(ax2, amplitudes_above_mean_x2[1:1200_000], linewidth = 0.5, color = :black);
display(GLMakie.Screen(), f)

#= array_IEI_x1 = get_IEI(t_EEs_x1)
array_PDF_IEI_x1 = get_PDF_IEI(array_IEI_x1; shift = 10)

array_IEI_x2 = get_IEI(t_EEs_x2)
array_PDF_IEI_x2 = get_PDF_IEI(array_IEI_x2; shift = 10)

f = Figure()
ax_x1 = Axis(f[1, 1], xlabel = L"peaks_EEs_x1", ylabel = L"IEI_{x1}")
ax_x2 = Axis(f[2, 1], xlabel = L"peaks_EEs_x2", ylabel = L"IEI_{x2}")
scatter!(ax_x1, peaks_EEs_x1[2:end], array_IEI_x1, markersize = 1.0, color = :blue)
scatter!(ax_x2, peaks_EEs_x2[2:end], array_IEI_x2, markersize = 1.0, color = :green)
display(GLMakie.Screen(), f) =#