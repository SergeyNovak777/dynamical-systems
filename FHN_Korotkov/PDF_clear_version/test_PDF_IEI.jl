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

using StaticArrays, DifferentialEquations, Statistics, Plots

include("/home/sergey/work/repo/dynamical-systems/system.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl")
#include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/plotter.jl")

function clear_work_space(sol, len_sol, ttr)
    sol = nothing
    len_sol = nothing
    ttr = nothing
    return sol, len_sol, ttr
end

t_truncate(t) = floor(Int64, t / 2)
Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function solver_IEI(prob, integrator_setting;
    plot_pdf_hist = true, bins = 25,
    plot_pdf_scatter = false,
    plot_timesereis = false, window_time_series = "nothing", t_start = 0, t_end = 0)

    sol = solve(prob, integrator_setting.alg, adaptive = true, abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
    maxiters = integrator_setting.max_iters, dense = integrator_setting.dense);
    println("sol_t_end: $(sol.t[end])")
    len_sol = length(sol.t)
    ttr = t_truncate(len_sol)
    data = [sol[1, ttr:len_sol], sol.t[ttr:len_sol]]
    len_data = length(data[1])
    println("len sol: $len_sol")
    println("len sol after truncate: $(len_sol - ttr)")

    sol, len_sol, ttr = clear_work_space(sol, len_sol, ttr)
    GC.gc()

    data_local_max = get_local_max(data)
    data_local_min = get_local_min(data)
    drop_artifacts(data_local_max, data_local_min)

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
    array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 20)

    if plot_pdf_hist == true
        plot_PDF_hist(array_IEI, array_PDF_IEI, bins)
    end
    if plot_pdf_scatter == true
        plot_PDF_scatter(array_IEI, array_PDF_IEI)
    end

    if plot_timesereis == true
        
        if window_time_series == "end"
            t_start = len_data - 100000;
            t_end = len_data
        elseif window_time_series == "start"
            t_start = 1;
            t_end = 10000;
        elseif window_time_series == "full"
            t_start = 1;
            t_end = len_data;
        end
        println("t_end: $(data[2][end])")
        plot_timeseries(data, data_local_max, data_local_min,
        t_peaks_spikes, peaks_spikes, Hs_x,
        t_EEs, peaks_EEs,
        t_start, t_end; clear_plot = true)
    end
end

parameters = FHN2_try3_params()
tspan = (0.0, 500_000.0)
parameters[7] = 0.09
parameters[8] = 75.74
u0 = [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525] 
u0 = SVector{5}(u0)
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)

alg = Vern9();
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol, max_iters = max_iters, dense = false)

solver_IEI(prob, integrator_setting;
plot_pdf_hist = true,
plot_pdf_scatter = false, plot_timesereis = true, window_time_series = "end"); GC.gc()                                                                                                                                                                                                                                                                                                                                                                                                                                                                        