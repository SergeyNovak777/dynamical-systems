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
using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, JLD2

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
    return sol
end

function calculate_LSE(sys, u0, params, integrator_setting, time_calculate)
    ds = CoupledODEs(sys, u0, params,
    diffeq = integrator_setting);
    LSE = lyapunovspectrum(ds, time_calculate)
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

function plot_phase_space(X, labels, lebelsize, ticksize, lw, resolution)
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
#---------------------------------------------------------------------------------------

sys = FHN2_try3
params = FHN2_try3_params()
tspan_for_solver = (0.0, 5000.0)
time_calculate_LSE = 10000
params[3] = 0.1
params[7] = 0.09
params[8] = 0.05
#u0 = get_u0(-5.0,-0.7,-10.0,-0.3)
u0 = (1.0, 0.0, -1.5, 0.4, -0.3)
u0 = SVector{5}(u0)
integrator_setting = (alg = Vern6(), adaptive = true, abstol = 1e-9, reltol = 1e-9)

solution = solver(sys, u0, params, integrator_setting, tspan_for_solver)
LSE = calculate_LSE(sys, solution[end], params, integrator_setting, time_calculate_LSE)

index_variable_timeseries = 1
percent_ploting_timeseries = 50

width_window, height_window = [1000, 300]
resolution_timeseries = (width_window, height_window)
ylabel, labelsize, ticksize, lw = [L"x_1", 40, 25, 1.0]

trange, x_range = get_timeseries(solution, index_variable_timeseries, percent_ploting_timeseries)
plot_timeseries(trange, x_range, ylabel, labelsize, ticksize, lw, resolution_timeseries)

percent_plotting_phase_space = 99
indexs = [1, 3, 4]
labels = [L"x_1", L"x_2", L"y_1"]
resolution_phase_space = (600, 900)
X = get_phase_space(solution, percent_plotting_phase_space, indexs)
plot_phase_space(X, labels, labelsize, ticksize, lw, resolution_phase_space)