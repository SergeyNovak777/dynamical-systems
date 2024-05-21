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

using StaticArrays, DifferentialEquations, Statistics, JLD2 

function clear_work_space(sol, len_sol, ttr)
    sol = nothing
    len_sol = nothing
    ttr = nothing
    return sol, len_sol, ttr
end

t_truncate(t) = floor(Int64, t / 2)

function save_data(iteration, data, path_to_save)
    file_name = "$(iteration)_sol_x1.jld"
    full_path = path_to_save * file_name
    jldsave(full_path; data)
end

function first_iteration(u0, parameters, tspan, integrator_setting, path_to_save)

    prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
    sol = solve(prob, integrator_setting.alg, adaptive = true,
    abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
    maxiters = integrator_setting.max_iters);
    last_point = sol[end]
    println("last point: $(last_point)"); flush(stdout)
    len_sol = length(sol.t)

    ttr = t_truncate(len_sol); tend = len_sol
    data_x1 = [sol[1, ttr:tend], sol.t[ttr:tend]]
    sol, len_sol, ttr = clear_work_space(sol, len_sol, ttr)
    GC.gc()

    save_data(1, data_x1, path_to_save)
    return last_point
end

function calculate_timeseris(u0_start, parameters, integrator_setting,
    t_point, count_iteration, path_to_save)

    global tspan = (0.0, t_point)
    global u0 = u0_start

    u0 = first_iteration(u0, parameters, tspan, integrator_setting, path_to_save)
    println("----------------------------------------------"); flush(stdout)
    println("");flush(stdout)

    for iteration in range(2, count_iteration)

        println("ITERATION: $(iteration)")

        tspan = (t_point*(iteration-1), t_point*iteration)
        prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
        println("u0: $(u0)"); flush(stdout)
        sol = solve(prob, integrator_setting.alg, adaptive = true,
        abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
        maxiters = integrator_setting.max_iters);
        u0 = sol[end]

        data_x1 = [sol[1, :], sol.t]
        sol = nothing; GC.gc();
        save_data(iteration, data_x1, path_to_save)

        println("last point: $(u0)"); flush(stdout)
        println("----------------------------------------------"); flush(stdout)
        println("");flush(stdout)
    end
end

alg = Vern9()
max_iters = 1e8
abs_tol = 1e-14; rel_tol = 1e-14
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol, max_iters = max_iters)

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/timeseries_k2_75_72/"
parameters = FHN2_try3_params()
parameters[7] = 0.09
parameters[8] = 75.72

u0_start = [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35]
u0_start = SVector{5}(u0_start)

t_point = 50000.0
count_iteration = 30

calculate_timeseris(u0_start, parameters, integrator_setting,
    t_point, count_iteration, path_to_save)