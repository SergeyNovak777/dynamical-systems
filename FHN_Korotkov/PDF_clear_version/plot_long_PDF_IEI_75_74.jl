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
filename_peaks = "spikes_peaks.jld2"
filename_t_peaks = "spikes_t_peaks.jld2"

peaks_spikes = load(path_to_folder*filename_peaks)
t_peaks_spikes = load(path_to_folder*filename_t_peaks)

println("count spikes: $(length(peaks_spikes))")


Hs_x = Hs(amplitudes_above_mean, 8)

index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
peaks_EEs = peaks_spikes[index_EEs]
t_EEs = t_peaks_spikes[index_EEs]
println("count EEs: $(length(peaks_EEs))")

array_IEI = get_IEI(t_EEs)
array_IEI = sort(array_IEI)
array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 10)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
hist!(ax, array_IEI, weights = array_PDF_IEI, bins = 25)#length(array_IEI))
display(GLMakie.Screen(), f)

# all amplitudes
#=f = Figure(size = (width_window, height_window))
ax = Axis(f[1, 1], ylabel = L"amplitudes")
lines!(ax, all_amplitudes)
hlines!(ax, Hs_x, linewidth = 1.0, linestyle = :dash, color = :red)
display(GLMakie.Screen(), f)=#

#=f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
lines!(ax, array_IEI, array_PDF_IEI, linewidth = 1.0)
display(GLMakie.Screen(), f)=#