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

using StaticArrays, DifferentialEquations, Statistics, CairoMakie, GLMakie, Distributions

function clear_work_space(sol, len_sol, ttr)
    sol = nothing
    len_sol = nothing
    ttr = nothing
    return sol, len_sol, ttr
end

t_truncate(t) = floor(Int64, t / 2)
Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function plot_PDF_hist(array_IEI, array_PDF_IEI, bins_)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
    hist!(ax, array_IEI, weights = array_PDF_IEI, bins = bins_) #length(array_IEI))
    display(GLMakie.Screen(), f)
end

function plot_amplitudes(all_amplitudes, Hs_x; width_window = 1000, height_window = 300)
    f = Figure(size = (width_window, height_window))
    ax = Axis(f[1, 1], ylabel = L"amplitudes")
    lines!(ax, all_amplitudes)
    hlines!(ax, Hs_x, linewidth = 1.0, linestyle = :dash, color = :red)
    display(GLMakie.Screen(), f)
end

function plot_PDF_scatter(array_IEI, array_PDF_IEI; width_window = 1000, height_window = 300)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
    scatter!(ax, array_IEI, array_PDF_IEI)
    display(GLMakie.Screen(), f)
end

function plot_timeseries(data, data_local_max, data_local_min, t_peaks_spikes, peaks_spikes,
    t_start, t_end;
    width_window = 1000, height_window = 300)
    f = Figure(size = (width_window, height_window))
    ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_1")
    lines!(ax, data[2][t_start:t_end], data[1][t_start:t_end])
    scatter!(ax, data_local_max[2], data_local_max[1], color = :green, markersize = 10)   
    scatter!(ax, t_peaks_spikes, peaks_spikes, color = :orange, markersize = 10)                                                                                                                                        
    scatter!(ax, data_local_min[2], data_local_min[1], color = :blue, markersize = 10)
    hlines!(ax, Hs_x, linewidth = 1.0, linestyle = :dash, color = :red)
    scatter!(ax, t_EEs, peaks_EEs, color = :deeppink, markersize = 10)
    xlims!(ax, data[2][t_start], data[2][t_end])
    display(GLMakie.Screen(), f)
end


function solver_IEI(prob, integrator_setting;
    plot_pdf_hist = true, bins = 25,
    plot_pdf_scatter = false,
    plot_timesereis = false, t_start = 1, t_end = 10000)
    sol = solve(prob, integrator_setting.alg, adaptive = true, abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
    maxiters = integrator_setting.max_iters);

    len_sol = length(sol.t)
    ttr = t_truncate(len_sol)
    data = [sol[1, ttr:len_sol], sol.t[ttr:len_sol]]
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
#=     t_start = length(data[2])-10000; t_end = length(data[2])

    plot_timeseries(data,
                    data_local_max, data_local_min,
                    t_peaks_spikes, peaks_spikes,
                    t_start, t_end)
    plot_PDF_hist(array_IEI, array_PDF_IEI; bins_ = 25) =#
end

parameters = FHN2_try3_params()
tspan = (0.0, 75000.0)
parameters[7] = 0.09
parameters[8] = 75.74
u0 = [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525] 
u0 = SVector{5}(u0)
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)

alg = Vern9();
abs_tol = 1e-14;
rel_tol = 1e-14;
max_iters = 1e8;
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol, max_iters = max_iters)

solver_IEI(prob, integrator_setting; plot_pdf_scatter = true)