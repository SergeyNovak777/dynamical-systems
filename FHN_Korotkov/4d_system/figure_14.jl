username = "irrito"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")
include(pathtorepo * "/system.jl")


using StaticArrays, SciMLBase, OrdinaryDiffEqVerner, CairoMakie, GLMakie

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9();
adaptive = true;
abs_tol = 1e-13;
rel_tol = 1e-13;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[3] = 0.1;
parameters[7] = 0.09;
parameters[8] = 75.74;
filename_hist = "route1_h_k2=$(parameters[8]).eps"

u0_start = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];
u0_start = SVector{4}(u0_start);

t_end = 30_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

path_to_save = "/home/irrito/MEGA2/dynamical_systems/FHN/images/"

labelsize = 75;
ticksize = 45;
CairoMakie.activate!();
length_sol = length(sol);
ttr = t_truncate(length_sol)
t_plot_start =  ttr;
t_plot_end = t_plot_start +130_000; # 150_000
xlimits = [-2, -1, 0, 1, 2] #[-1.1, -0.9, -0.75]
zlimits = [-0.5, 0, 0.5] #[-0.63, -0.61, -0.59]
indexx = 1; indexy  = 3; indexz = 4;
f = Figure(size = (1000 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_2",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 150, protrusions = (20, 20,100, 20),
    xticks = xlimits, yticks = xlimits, zticks = zlimits)

lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 0.7, color = :black);
display(GLMakie.Screen(), f);

save(path_to_save*filename_hist, f)