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
filename_data_local_max = "data_local_max_x1.jld2"
filename_data_local_min = "data_local_min_x1.jld2"

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
array_IEI = sort(array_IEI)


shift_ = 0.03
sorted_amplitudes = sort(amplitudes_above_mean[1:10_000_000])

array_PDF_amplitudes = get_PDF_IEI(sorted_amplitudes; shift = shift_)
array_PDF_amplitudes_without_less_shift = get_PDF_without_less_shift(sorted_amplitudes, shift_)

f = Figure()
ax = Axis(f[1, 1], title = "PDF shift to both side")
lines!(ax, sorted_amplitudes, array_PDF_amplitudes, linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)
f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10, title = "PDF shift to both side")
hist!(ax, sorted_amplitudes, weights = array_PDF_amplitudes, bins = 30)
display(GLMakie.Screen(), f)


f = Figure()
ax = Axis(f[1, 1], title = "PDF shift to up side")
lines!(ax, sorted_amplitudes, array_PDF_amplitudes_without_less_shift, linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)
f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10, title = "PDF shift to up side")
hist!(ax, sorted_amplitudes, weights = array_PDF_amplitudes_without_less_shift, bins = 30)
display(GLMakie.Screen(), f)

#= f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
lines!(ax, array_IEI, array_PDF_IEI, linewidth = 1.0)
display(GLMakie.Screen(), f) =#

#= f = Figure() 
ax = Axis(f[1, 1])
lines!(ax, all_amplitudes, linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f) =#