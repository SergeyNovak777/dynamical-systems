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
end

using StaticArrays, Statistics, GLMakie, JLD2

include("/home/sergey/work/repo/dynamical-systems/system.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/plotter.jl")

function concat_vectors(sol, sol_times, len_charts)
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

function analys_timeseries(sol, sol_times;
    shift_PDF = 10,
    plot_pdf_hist = true, bins = 25,
    plot_pdf_scatter = false,
    plot_timesereis = false, window_time_series = "nothing", t_start = 0, t_end = 0)

    data = [sol, sol_times]
    len_data = length(sol)
    sol = nothing; sol_times = nothing; GC.gc();
    data_local_max = get_local_max(data)
    data_local_min = get_local_min(data)

    println("len timeseries: $(length(data))")
    println("length local maximums before drop artifacts: $(length(data_local_max[1]))")
    println("length local minimums before drop artifacts: $(length(data_local_min[1]))")

    drop_artifacts(data_local_max, data_local_min)

    println("length local maximums after drop artifacts: $(length(data_local_max[1]))")
    println("length local minimums after drop artifacts: $(length(data_local_min[1]))")

    all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
    amplitudes_threshold = Statistics.mean(all_amplitudes)
    peaks_spikes, t_peaks_spikes, amplitudes_above_mean = select_spikes(data_local_min[1], data_local_max, amplitudes_threshold)
    
    Hs_x = Hs(amplitudes_above_mean, 8)
    index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
    peaks_EEs = peaks_spikes[index_EEs]
    t_EEs = t_peaks_spikes[index_EEs]

    println("count spikes: $(length(peaks_spikes))")
    println("minumal amplitude: $(minimum(all_amplitudes))")
    println("maximal amplitude: $(maximum(all_amplitudes))")
    println("threshold for detect spike: $(amplitudes_threshold)")
    println("count EEs: $(length(peaks_EEs))")

    array_IEI = get_IEI(t_EEs)
    array_IEI = sort(array_IEI)
    array_PDF_IEI = get_PDF_IEI(array_IEI; shift = shift_PDF)

    if plot_pdf_hist == true
        plot_PDF_hist(array_IEI, array_PDF_IEI, bins)
    end
    if plot_pdf_scatter == true
        plot_PDF_scatter(array_IEI, array_PDF_IEI)
    end

    if plot_timesereis == true
        
        if window_time_series == "end"
            t_start = len_data - 10000;
            t_end = len_data
        elseif window_time_series == "start"
            t_start = 1;
            t_end = 10000;
        elseif window_time_series == "full"
            t_start = 1;
            t_end = len_data;
        end
        plot_timeseries(data, data_local_max, data_local_min,
        t_peaks_spikes, peaks_spikes, Hs_x,
        t_EEs, peaks_EEs,
        t_start, t_end)
    end
end

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/timeseries_k2_75_74_change_tol/"
len_charts_timeseries = 15

len_part_timeseris = Int64[]
sol = Float64[]
sol_times = Float64[]
sol, sol_times = concat_vectors(sol, sol_times, len_charts_timeseries)

analys_timeseries(sol, sol_times; plot_pdf_hist = false, plot_timesereis = true, window_time_series = "end")