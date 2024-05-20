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
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/pdf_function.jl")
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/preprocessing_spike.jl")
end

using StaticArrays, DifferentialEquations, DynamicalSystems, Statistics, CairoMakie, GLMakie, JLD2

function clear_workspace()
    sol=nothing
    data_x1 = nothing
    array_spikes_max = nothing
    array_t_spikes_max = nothing
    array_spikes_thresholds = nothing
    array_t_spikes_thresholds = nothing
    EE = nothing
    t_EE = nothing
    GC.gc()
    sleep(3)
end

function save_data(iteration, data_x1, path_to_save)
    namefile_sol_x1 = "$(iteration)_sol_x1.jld"
    #namefile_EE = "$(iteration)_EE.jld"
    full_path_to_save_sol_x1 = path_to_save*namefile_sol_x1
    #full_path_to_save_EE = path_to_save*namefile_EE
    jldsave(full_path_to_save_sol_x1; data_x1)
    #jldsave(full_path_to_save_EE, t_EE = t_EE, peaks_EE = EE)
end

function first_iteration(u0, parameters, tspan, path_to_save)

    prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
    sol = solve(prob, alg, adaptive = true, abstol = abs_tl, reltol = rel_tl, maxiters = max_iters);
    last_point = sol[end]
    println("last point: $(last_point)"); flush(stdout)
    len_sol = length(sol.t)
    tstart = t_truncate(len_sol); tend = len_sol

    data_x1 = [sol[1, tstart:tend], sol.t[tstart:tend]]
    println("ITERATION 1"); flush(stdout)
    #=array_spikes_max, array_t_spikes_max = get_peaks(data_x1; level_zero = "nothing")
    array_spikes_thresholds, array_t_spikes_thresholds = get_peaks_neg(data_x1)
    println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))"); flush(stdout)

    if check_timeseries(array_spikes_max, array_spikes_thresholds) == true
        println("good timeseries"); flush(stdout)
        array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
    else
        println("bad timeseries; drop false start or end spike"); flush(stdout)
        println("DROPED: length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))"); flush(stdout)
        drop_false_start_end(array_t_spikes_max, array_spikes_max,  array_t_spikes_thresholds)
        array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
    end
    
    println("count correct spike: $(length(array_spikes_max_correct))"); flush(stdout)
    
    Hs_x  = Hs_above(amplitudes, 8)
    index_EE = findall(x-> x >= Hs_x, array_spikes_max_correct)
    t_EE = array_t_spikes_max_correct[index_EE]
    EE = array_spikes_max_correct[index_EE]
    println("count EE: $(length(EE))")=#
    save_data(1, data_x1, path_to_save)

    return last_point
end

function calculate_timeseris(u0_start, parameters, t_point, count_iteration, path_to_save)

    global tspan = (0.0, t_point)
    global u0 = u0_start

    u0 = first_iteration(u0, parameters, tspan, path_to_save)
    println("CALL GC"); flush(stdout)
    GC.gc(); sleep(3)
    println("----------------------------------------------"); flush(stdout)
    println("");flush(stdout)
    for iteration in range(2, count_iteration)

        println("ITERATION: $(iteration)")

        tspan = (t_point*(iteration-1), t_point*iteration)
        prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
        println("u0: $(u0)"); flush(stdout)
        sol = solve(prob, alg, adaptive = true, abstol = abs_tl, reltol = rel_tl, maxiters = max_iters);
        u0 = sol[end]
        println("last point: $(u0)"); flush(stdout)
        data_x1 = [sol[1, :], sol.t]
        #=array_spikes_max, array_t_spikes_max = get_peaks(data_x1; level_zero = "nothing")
        array_spikes_thresholds, array_t_spikes_thresholds = get_peaks_neg(data_x1)
        println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))"); flush(stdout)

        if check_timeseries(array_spikes_max, array_spikes_thresholds) == true
            println("good timeseries"); flush(stdout)
            array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
        else
            println("bad timeseries; drop false start or end spike"); flush(stdout)
            println("DROPED: length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))"); flush(stdout)
            drop_false_start_end(array_t_spikes_max, array_spikes_max,  array_t_spikes_thresholds)
            array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
        end
        
        println("count correct spike: $(length(array_spikes_max_correct))"); flush(stdout)
        
        Hs_x  = Hs_above(amplitudes, 8)
        index_EE = findall(x-> x >= Hs_x, array_spikes_max_correct)
        t_EE = array_t_spikes_max_correct[index_EE]
        EE = array_spikes_max_correct[index_EE]
        println("count EE: $(length(EE))")=#
        save_data(iteration, data_x1, path_to_save)

        println("CALL clear_workspace"); flush(stdout)
        clear_workspace()
        println("----------------------------------------------"); flush(stdout)
        println("");flush(stdout)
    end
end

alg = Vern9()
max_iters = 1e8
abs_tl = 1e-14; rel_tl = 1e-14
path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/timeseries_k2_75/"
parameters = FHN2_try3_params()
parameters[7] = 0.09
parameters[8] = 75.73
u0_start = [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35]
u0_start = SVector{5}(u0_start)
t_point = 50000.0
count_iteration = 30

calculate_timeseris(u0_start, parameters, t_point, count_iteration, path_to_save)