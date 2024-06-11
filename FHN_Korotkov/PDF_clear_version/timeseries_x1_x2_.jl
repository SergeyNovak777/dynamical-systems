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
end

using StaticArrays, DifferentialEquations, JLD2, CairoMakie, GLMakie

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9()
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
println("alg: $alg"); println("abstol: $abs_tol; reltol: $(rel_tol)")
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol, max_iters = max_iters)

path_to_save = "/home/sergey/timeseries_k2_75_74_save_x1_x2/"
parameters = FHN2_try3_params()
parameters[7] = 0.09
parameters[8] = 75.74

u0_start = [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35]   
u0_start = SVector{5}(u0_start)

t_point = 200_000.0
tspan = (0.0, t_point)


prob = ODEProblem(FHN2_try3, u0_start, tspan, parameters)

sol = solve(prob, integrator_setting.alg, adaptive = true,
    abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
    maxiters = integrator_setting.max_iters);

len_sol = length(sol.t)
ttr = t_truncate(len_sol); tend = len_sol

labelsize = 25;
ticksize = 15;

f = Figure()

ax1 = Axis(f[1, 1], xlabel = L"t", ylabel = L"x_2",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)

#lines!(ax1, sol.t[ttr:tend], sol[1, ttr:tend], linewidth = 1.0, color = :red)
lines!(ax1, sol.t[ttr:tend], sol[3, ttr:tend], linewidth = 1.0, color = :green)

display(GLMakie.Screen(), f)