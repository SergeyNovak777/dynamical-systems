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

# INCLUDE
#---------------------------------------------------------------------
include("/home/sergey/work/repo/dynamical-systems/system.jl")
#----------------------------------------------------------------------
#PACKAGES
using StaticArrays, DifferentialEquations, DynamicalSystems
using CairoMakie
using JLD2
#----------------------------------------------------------------------

global ticksize = 25
global labelsize = 40
global lw = 2.0

function plot_two_timeseries(sol, tstart, tend, indexs, labels, path_to_save)
    indexx, indexy = indexs
    f = Figure(size = (1000, 250))
    ax = Axis(f[1,1], xlabel = L"time", ylabel = labels, xlabelsize = labelsize, ylabelsize = labelsize, xticklabelsize = ticksize, yticklabelsize = ticksize)
    lines!(ax, sol.t[tstart:tend], sol[indexx, tstart:tend], color = :red, linewidth = lw)
    lines!(ax, sol.t[tstart:tend], sol[indexy, tstart:tend], color = :green, linewidth = lw)
    display(f)
    save(path_to_save, f)
end

function plot_phase_space(sol, tstart, tend, indexs, labels, path_to_save)
    f = Figure(size = (900, 600))
    ax = Axis3(f[1, 1], xlabel = labels[1], ylabel = labels[2], zlabel = labels[3], azimuth = -0.3pi,
        xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize, xticklabelsize= ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize)
    lines!(ax, sol[indexs[1], tstart:tend], sol[indexs[2], tstart:tend], sol[indexs[3], tstart:tend])
    display(f)
    save(path_to_save, f)
end

function get_u0(x1, y1, x2, y2)
    z = y1 - y2
    return SVector{5}([x1, y1, x2, y2, z])
end

function main(x1, y1, x2, y2)
    parameters = FHN2_try3_params()
    parameters[7] = 0.0
    parameters[8] = 0.0

    tspan = (0.0, 1000)

    integrator_setting = (alg = DP8(), adaptive = true, abstol = 1e-11, reltol = 1e-11)

    u0 = get_u0(1.0, 0.0, 0.01, -1.0)
    prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
    sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
    abstol = integrator_setting.abstol, reltol = integrator_setting.abstol)

    len_sol = length(sol.u[:,1])
    tstart = floor(Int64, len_sol / 2); tend = len_sol

    pathtosave = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/figure 3/"
    filename_timeseries = "timeseries_figure_3.eps"
    filename_phase_space = "phase_space_figure_3.eps"

    plot_two_timeseries(sol, tstart, tend, [1, 3], L"x_i", pathtosave * filename_timeseries)
    plot_phase_space(sol, tstart, tend, [1,3,4], [L"x_1", L"x_2", L"y_2"], pathtosave * filename_phase_space)
end