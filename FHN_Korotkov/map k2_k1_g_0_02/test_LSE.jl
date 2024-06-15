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
using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

function solver(prob, integrator_setting)
    if integrator_setting.adaptive == true
        sol = solve(prob, integrator_setting.alg, adaptive = true, abstol = integrator_setting.abstol, reltol  = integrator_setting.reltol,
        maxiters = integrator_setting.maxiters)
    else
        sol = solve(prob, integrator_setting.alg, adaptive = false, dt = integrator_setting.dt, maxiters = integrator_setting.maxiters)
    end
    return sol
end

sys = FHN2_try3

#=
butterfly

u0 = [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525] 

params[3] = 0.02
params[7] = 0.005
params[8] = 1.0

params[3] = 0.02
params[7] = 0.005
params[8] = 0.1 
=#

params = FHN2_try3_params()
params[3] = 0.02
params[7] = 0.004
params[8] =  0.75

tspan = (0.0, 10000.0)
time_calculate_LSE = 5000

u0 = [-1.122685745876286, -0.5655795683079851, -1.2279038564970448, -0.538136579127275, -0.027442989180487945]
#[-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525] 
u0 = SVector{5}(u0)

integrator_setting = (alg = DP8(), adaptive = true, abstol = 1e-10, reltol = 1e-10, maxiters = 5e7)

prob = ODEProblem(sys, u0, tspan, params)
sol = solver(prob, integrator_setting)

ds = CoupledODEs(sys, sol[end], params,
diffeq = integrator_setting);
LSE = lyapunovspectrum(ds, time_calculate_LSE)
println("LSE: $LSE")

t_plot_start = 30000; t_plot_end = 50000;
f = Figure()
ax = Axis(f[1, 1], xlabel = L"x_1", ylabel = L"x_2")
lines!(ax, sol[1, t_plot_start:t_plot_end], sol[3, t_plot_start:t_plot_end], linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)

f = Figure()
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_1")
lines!(ax, sol[1, t_plot_start:t_plot_end], sol[3, t_plot_start:t_plot_end], sol[2, t_plot_start:t_plot_end], linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)