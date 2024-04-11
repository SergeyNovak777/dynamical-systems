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
include("/home/sergey/work/repo/dynamical-systems/system.jl")


#using Pkg
#Pkg.activate("D:\\work\\dynamical-systems\\env\\integrate")
#include("D:\\work\\dynamical-systems\\system.jl")


using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, JLD2, BenchmarkTools

function get_u0(x1, y1, x2, y2)
    z = y1 - y2
    return SVector{5}([x1, y1, x2, y2, z])
end

function get_percent(number, percent)
    return floor(Int64, (number / 100) * percent )
end

function solver(sys, u0, params, integrator_setting, tspan)

    prob = ODEProblem(sys, u0, tspan, params)

    if integrator_setting.adaptive == true
        sol = solve(prob, integrator_setting.alg, adaptive = true, abstol = integrator_setting.abstol, reltol  = integrator_setting.reltol)
    else
        sol = solve(prob, integrator_setting.alg, adaptive = false, dt = integrator_setting.dt, maxiters = integrator_setting.maxiters)
    end    
    return prob, sol
end

function calculate_LSE(sys, u0, params, integrator_setting, time_calculate)
    ds = CoupledODEs(sys, u0, params,
    diffeq = integrator_setting);
    LSE = lyapunovspectrum(ds, time_calculate)
    return LSE
end

function calculate_LLE(sys, u0, params, integrator_setting, time_calculate)
    ds = CoupledODEs(sys, u0, params,
    diffeq = integrator_setting);
    LSE = lyapunov(ds, time_calculate)
    return LSE
end

function get_timeseries(solution, index_variable, percent)
    len_solution = length(solution)
    start_time_plot_timeseries =  get_percent( len_solution, 100 - percent )
    end_time_plot_timeseries = len_solution
    trange = solution.t[start_time_plot_timeseries:end_time_plot_timeseries]
    x_range = solution[index_variable, start_time_plot_timeseries:end_time_plot_timeseries]
    return trange, x_range
end

function plot_timeseries(trange, x_range, ylabel_, labelsize, ticksize, lw, resolution_)
    fig = Figure(size = resolution_)
    ax = Axis(fig[1, 1], xlabel = L"time", ylabel = ylabel_, xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize)
    lines!(ax, trange, x_range, linewidth = lw)
    display(fig)
end

function plot_phase_space(X, labels, labelsize, ticksize, lw, resolution)
    fig = Figure(size = resolution)
    ax = Axis3(fig[1, 1], xlabel = labels[1], ylabel = labels[2], zlabel = labels[3],
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize)
    lines!(ax, X[1], X[2], X[3], linewidth = lw)
    display(fig)
end

function get_phase_space(solution, percent_plotting_phase_space, indexs)
    len_solution = length(solution)
    start_time =  get_percent( len_solution, 100 - percent_plotting_phase_space )
    end_time = len_solution
    X = [solution[indexs[1], start_time:end_time], solution[indexs[2], start_time:end_time], solution[indexs[3], start_time:end_time]]
    return X
end

function difference_between_points(solution, len_matrix)
    matrix_difference = zeros(len_matrix, 5)
    for index in range(1, len_matrix - 1, step = 2)
        difference = solution[index] - solution[index+1]
        difference = abs.(collect(difference))
        matrix_difference[index, :] = difference
    end
    return matrix_difference
end
#---------------------------------------------------------------------------------------

function show_mode(g, k1, k2, u0, static)
    sys = FHN2_try3
    params = FHN2_try3_params()
    tspan_for_solver = (0.0, 3000.0)
    time_calculate_LSE = 20000
    params[3] = g
    params[7] = k1
    params[8] = k2
    #u0 = get_u0(-5.0,-0.7,-10.0,-0.3)
    #u0 = SVector{5}(u0)
    integrator_setting = (alg = DP8(), adaptive = true, abstol = 1e-10, reltol = 1e-10)

    if static == false
        prob, solution = solver(sys, SVector{5}(u0), params, integrator_setting, tspan_for_solver)
    else
        prob, solution = solver(sys, u0, params, integrator_setting, tspan_for_solver)
    end
    
    ϵ=1e-9
    len_sol = length(solution)
    percent = get_percent(len_sol, 100-50)
    len_matrix = length(solution[:, percent:end])
    percent_sol = solution[:, percent:end]
    if mod(len_matrix, 2) != 0
        len_matrix = len_matrix - 1
    end

    matrix_difference = difference_between_points(percent_sol, len_matrix)
    if all(matrix_difference .<= ϵ) == true
        LSE = -1
        LLE = -1
    else
        LSE = calculate_LSE(sys, u0, params, integrator_setting, time_calculate_LSE)
        LLE = calculate_LLE(sys, u0, params, integrator_setting, time_calculate_LSE)
    end
    println("LSE: $LSE")
    println("LLE: $LLE")

    index_variable_timeseries = 1
    percent_ploting_timeseries = 10

    width_window, height_window = [1000, 300]
    resolution_timeseries = (width_window, height_window)
    ylabel, labelsize, ticksize, lw = [L"x_1", 40, 25, 1.0]

    trange, x_range = get_timeseries(solution, index_variable_timeseries, percent_ploting_timeseries)
    plot_timeseries(trange, x_range, ylabel, labelsize, ticksize, lw, resolution_timeseries)

    percent_plotting_phase_space = 50
    indexs = [1, 3, 4]
    labels = [L"x_1", L"x_2", L"y_1"]
    resolution_phase_space = (600, 900)
    X = get_phase_space(solution, percent_plotting_phase_space, indexs)
    plot_phase_space(X, labels, labelsize, ticksize, lw, resolution_phase_space)

    return prob, solution, matrix_difference
end


g = 0.03522232323232323 #0.01109
k1 = 0.13717171717171717
k2 = 1.0
#u0 = [-5.0,-0.7,-10.0,-0.3, -0.7+0.3]
u0 = [1.9075684503044907, -0.3904496392396389, 1.9413139567288633, -0.4925421548500994, 0.10209251561045568]
prob, solution, matrix_difference = show_mode(g, k1, k2, u0, false);
;

#=
index_cycle_k_1: 98; value_k_1: 0.13717171717171717
index_cycle_g: 11; value_g: 0.03522232323232323
init_point: [1.9075684503044907, -0.3904496392396389, 1.9413139567288633, -0.4925421548500994, 0.10209251561045568]
last_point: [-1.0099999998314062, -0.6584764912780838, -1.0100000001685954, -0.6584764912567943, -2.1574027921234565e-11]
LSE: [0.0011149179338305387, -1.0885250606064384e-6, -4.303153192615448, -4.303166333295496, -6.959117692338693]
=#


#=ϵ=1e-9
len_sol = length(sol)
percent = get_percent(len_sol, 100-20)
len_matrix = length(sol[:, percent:end])
percent_sol = sol[:, percent:end]
if mod(len_matrix, 2) != 0
    len_matrix = len_matrix - 1
end

matrix_difference = difference_between_points(percent_sol, len_matrix)
if all(matrix_difference .<= ϵ) == true
    println("FP!!!!!")
end
=#
#=
CYCLE
g = 0.15
k1 = 0.14
k2 = 1.0
u0 = [1.9075684503044907, -0.3904496392396389, 1.9413139567288633, -0.4925421548500994, 0.10209251561045568]

g = 0.01109
k1 = 0.0
k2 = 1.0
u0 = [-1.0083078648152173, -0.6646900583246136, -1.0118497540232405, -0.6636637416235039, -0.0010263167011098397]
=#

