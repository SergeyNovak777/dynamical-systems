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

using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

function FHN2_4d(u, p ,t)
    x1, y1, x2, y2 = u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + (k1 + k2 * (y1 - y2)^2) * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + (k1 + k2 * (y1 - y2)^2) * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    return SVector(dx1dt, dy1dt, dx2dt, dy2dt)
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
parameters[8] = 76.8;
filename_hist = "fig_15_e.eps"

u0_start = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];
u0_start = SVector{4}(u0_start);

t_end = 10_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/rewrite_images/"


labelsize = 80;
ticksize = 50;
CairoMakie.activate!();
length_sol = length(sol);
ttr = t_truncate(length_sol)
t_plot_start =  ttr
t_plot_end = t_plot_start + 50_000;
           
xticks = [-1.13, -1.0, -0.88];
yticks = [-1.12, -0.97, -0.80];
zticks = [-0.64, -0.62, -0.60];
indexx = 1; indexy  = 3; indexz = 4;
f = Figure(size = (1100 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_2",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 160, protrusions = (30, 30, 120, 30),
    xticks = xticks, yticks = yticks, zticks = zticks)#,
    #xticks = [-0.7, -0.9, -1.12], yticks = [-0.7, -0.9, -1.12], zticks = [-0.635, -0.61, -0.58])

lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);

save(path_to_save*filename_hist, f)