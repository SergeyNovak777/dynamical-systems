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
using Base.Threads

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

path_to_folder = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/local_maxs/k2=75/"

filename_data_local_max = "data_local_max_x1.jld2"
filename_data_local_min = "data_local_min_x1.jld2"
data_local_max = load(path_to_folder*filename_data_local_max)["data_local_max"]
data_local_min = load(path_to_folder*filename_data_local_min)["data_local_min"]

array_PDF_ISI_parallel = load(path_to_folder*"PDF_ISI.jld2")["array_PDF_ISI_parallel"];
array_ISI_par = load(path_to_folder*"ISI.jld2")["array_ISI_par"];


all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
mean_amplitudes = Statistics.mean(all_amplitudes)
all_amplitudes = nothing;
println("mean amplitude: $mean_amplitudes")
peaks_spikes, t_peaks_spikes, _ = select_spikes(data_local_min[1], data_local_max, mean_amplitudes) # Statistics.std(all_amplitudes))

data_local_max = nothing; data_local_min = nothing; GC.gc();

println("count spikes: $(length(peaks_spikes))")
#Hs_x = Hs(amplitudes_above_mean, 8)


#= index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
peaks_EEs = peaks_spikes[index_EEs]
t_EEs = t_peaks_spikes[index_EEs]
println("count EEs: $(length(peaks_EEs))") =#

labelsize = 40;
ticksize = 30;

#= f = Figure()
ax = Axis(f[1, 1], xlabel = L"t", ylabel = L"peaks_{x1}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_peaks_spikes[1:20_000], peaks_spikes[1:20_000], linewidth = 1.0)
hlines!(ax, Hs_x, linewidth = 3.0, linestyle = :dash, color = :red)
display(GLMakie.Screen(), f) =#


function get_PDF_IEI_parallel(IEI; shift = 10)

    total_count_IEI = length(IEI)
    array_PDF = zeros(total_count_IEI)
    array_IEI_par = zeros(total_count_IEI)

    @threads for index in 1:total_count_IEI
        count_IEI_i = count(IEI[index]-shift .<= IEI .<= IEI[index]+shift)
        PDF_IEI_i = count_IEI_i / total_count_IEI
        array_PDF[index] = PDF_IEI_i
        array_IEI_par[index] = IEI[index]
    end

    return array_PDF, array_IEI_par
end


count_peaks = length(t_peaks_spikes);
println("t_end: $(t_peaks_spikes[count_peaks])");
#array_ISI = get_IEI(t_peaks_spikes[1:count_peaks]);
#array_IEI = sort(array_IEI)

#array_PDF_ISI = get_PDF_IEI(array_ISI; shift = 0.0)

#= array_PDF_ISI_parallel, array_ISI_par = get_PDF_IEI_parallel(array_ISI; shift = 0.0)

jldsave(path_to_folder*"PDF_ISI.jld2"; array_PDF_ISI_parallel)
jldsave(path_to_folder*"ISI.jld2"; array_ISI_par) =#

#= Hs_ISI_coeff_8 = Hs(array_ISI, 8);
Hs_ISI_coeff_6 = Hs(array_ISI, 6); =#

Hs_ISI_coeff_8_par = Hs(array_ISI_par, 8);
Hs_ISI_coeff_6_par = Hs(array_ISI_par, 6);

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/EEs/"
filename_hist = "k2=75_PDF_ISI_hist_x1.pdf"
filename_tEE_IEI = "k2=75_t_spikes_ISI_x1.pdf"
filename_peaks_EE_IEI = "k2=75_peaks_spikes_ISI_x1.pdf"

GLMakie.activate!()

# linear version
#= f1 = Figure()
ax = Axis(f1[1, 1], xlabel = L"ISI", ylabel = L"PDF_{ISI}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_ISI, weights = array_PDF_ISI, bins = 100)
vlines!(ax, Hs_ISI_coeff_8, linewidth = 5.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_ISI_coeff_6, linewidth = 5.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f1)

f2 = Figure()
ax = Axis(f2[1, 1], xlabel = L"t_{spike}", ylabel = L"ISI",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_peaks_spikes[2:count_peaks], array_ISI[1:end], linewidth = 1.0)
hlines!(ax, Hs_ISI_coeff_8, linewidth = 5.0, linestyle = :dash, color = :red, label = L"H_s=4612")
hlines!(ax, Hs_ISI_coeff_6, linewidth = 5.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f2) =#

#-----------------------------------------------------------------------------------
# parallel version
CairoMakie.activate!();
f1 = Figure()
ax = Axis(f1[1, 1], xlabel = L"ISI", ylabel = L"PDF_{ISI}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_ISI_par, weights = array_PDF_ISI_parallel, bins = 100)
vlines!(ax, Hs_ISI_coeff_8_par, linewidth = 5.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_ISI_coeff_6_par, linewidth = 5.0, linestyle = :dash, color = :green)
text!(ax, 25.5, 0.001, text = L"Hs_8", fontsize = labelsize)
text!(ax, 18, 0.001, text = L"Hs_6", fontsize = labelsize)
display(f1)
save(path_to_save*filename_hist, f1)


t_interval = 8_000_000;
f2 = Figure()
ax = Axis(f2[1, 1], xlabel = L"t_{spike}", ylabel = L"ISI",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_peaks_spikes[2:t_interval], array_ISI_par[1:t_interval-1], linewidth = 1.0)
hlines!(ax, Hs_ISI_coeff_8_par, linewidth = 5.0, linestyle = :dash, color = :red)
hlines!(ax, Hs_ISI_coeff_6_par, linewidth = 5.0, linestyle = :dash, color = :green)
display(f2)
save(path_to_save*filename_tEE_IEI, f2)