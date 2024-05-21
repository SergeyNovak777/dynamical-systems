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

function vcat_arrays(sol, sol_times, len_charts)
    for iteration in range(1, len_charts)

        namefile_sol_x1 = "$(iteration)_sol_x1.jld"
        full_path_to_save_sol_x1 = path_to_save*namefile_sol_x1

        sol_x1, sol_t = load(full_path_to_save_sol_x1)["data"]
        sol = vcat(sol, sol_x1)
        sol_times = vcat(sol_times, sol_t)
        push!(len_part_timeseris, length(sol_x1))
    end

    return sol, sol_times
end
Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/timeseries_k2_75_72/"

len_part_timeseris = Int64[]
sol = Float64[]
sol_times = Float64[]
len_charts = 1
sol, sol_times = vcat_arrays(sol, sol_times, len_charts)

data = [sol, sol_times]
println("len sol: $(length(data[1]))")

data_local_max = get_local_max(data)
data_local_min = get_local_min(data)

drop_artifacts(data_local_max, data_local_min)

println("length local maxs: $(length(data_local_max[1]))")
println("length local mins: $(length(data_local_min[1]))")

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
array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 10)

#PLOTTING
height_window = 400; width_window = 1100;
# timeseries
t_start =  length(data[2])-10000; t_end = length(data[2])
f = Figure(size = (width_window, height_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_1")
lines!(ax, data[2][t_start:t_end], data[1][t_start:t_end])
scatter!(ax, data_local_max[2], data_local_max[1], color = :green, markersize = 10)   
#scatter!(ax, t_peaks_spikes, peaks_spikes, color = :orange, markersize = 10)                                                                                                                                        
scatter!(ax, data_local_min[2], data_local_min[1], color = :blue, markersize = 10)
#hlines!(ax, Hs_x, linewidth = 5.0, linestyle = :dash, color = :red)
#scatter!(ax, t_EEs, peaks_EEs, color = :deeppink, markersize = 10)
xlims!(ax, data[2][t_start], data[2][t_end])
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
#=
f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
hist!(ax, array_IEI, weights = array_PDF_IEI, bins = 30)#length(array_IEI))
display(GLMakie.Screen(), f)=#