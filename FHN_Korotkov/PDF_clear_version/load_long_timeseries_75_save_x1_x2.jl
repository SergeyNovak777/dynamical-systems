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

function vcat_arrays(array_local_max, array_t_local_max,
                    array_local_min, array_t_local_min,
                    len_charts)

    for iteration in range(1, len_charts)

        namefile_sol_x1 = "$(iteration)_sol_x1.jld2"
        full_path_to_save_sol_x1 = path_to_save*namefile_sol_x1

        sol_x1, sol_t = load(full_path_to_save_sol_x1)["datax1"]
        data = [sol_x1, sol_t]
        
        array_local_max_lc, array_t_local_max_lc = get_local_max(data)
        array_local_min_lc, array_t_local_min_lc = get_local_min(data)

        data = nothing; GC.gc()

        array_local_max = vcat(array_local_max, array_local_max_lc)
        array_t_local_max = vcat(array_t_local_max, array_t_local_max_lc)
        array_local_min = vcat(array_local_min, array_local_min_lc)
        array_t_local_min = vcat(array_t_local_min, array_t_local_min_lc)

        push!(len_part_timeseris, length(sol_x1))
        println("iteration: $iteration"); flush(stdout)
        println("----------------------");flush(stdout)
    end

    return  array_local_max, array_t_local_max, array_local_min, array_t_local_min
    
end

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

path_to_save = "/home/sergey/timeseries_k2_75_save_x1_x2/"

len_part_timeseris = Int64[]

array_local_max = Float64[]
array_t_local_max = Float64[]
array_local_min = Float64[]
array_t_local_min = Float64[]

len_charts = 1000

array_local_max, array_t_local_max, array_local_min, array_t_local_min = vcat_arrays(array_local_max, array_t_local_max,array_local_min, array_t_local_min, len_charts)

data_local_max = [array_local_max, array_t_local_max]
data_local_min = [array_local_min, array_t_local_min]

drop_artifacts(data_local_max, data_local_min)

println("length local maxs: $(length(data_local_max[1]))")
println("length local mins: $(length(data_local_min[1]))")

#= all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
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
array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 10) =#

path_to_folder = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/local_maxs/k2=75/"
filename_data_local_max = "data_local_max_x1.jld2"
filename_data_local_min = "data_local_min_x1.jld2"

jldsave(path_to_folder*filename_data_local_max; data_local_max)
jldsave(path_to_folder*filename_data_local_min; data_local_min)