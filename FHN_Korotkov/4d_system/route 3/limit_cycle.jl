if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "irrito"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
    include("/home/irrito/work/repo/dynamical-systems/system.jl")
end

using StaticArrays, SciMLBase, OrdinaryDiffEqVerner, ChaosTools, CairoMakie, GLMakie

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

alg = Vern9();
adaptive = true;
abs_tol = 1e-12;
rel_tol = 1e-12;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[7] = 0.0935
parameters[8] = 64.76190476190476

u0_start = [-1.0836728460611933, -0.6318417392022484, -0.9017528537331925, -0.624049721609583]

u0_start = SVector{4}(u0_start);

t_end = 50_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_4d, sol[end], parameters,
diffeq = integrator_setting);

x1 = x2 = interval(-1.5, 1.5)
y1 = y2 = interval(-1.5, 1.5)
box = [x1, y1, x2, y2]

fixed_point, eigs, _ = fixedpoints(ds, box)
fixed_point = fixed_point[1];
labelsize = 85;
ticksize = 50;

t_plot_start = 1;
t_plot_end = 15_000; #len_sol;

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/scenario/"

CairoMakie.activate!();

f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"y_1", ylabel = L"y_2", zlabel = L"x_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 85, ylabeloffset = 85, zlabeloffset = 115,
    protrusions = (30, 30,120, 30),
    xticks = [-0.635, -0.622], yticks = [-0.635, -0.622], zticks = [-1.05, -0.95]);
lines!(ax, sol[2, t_plot_start:t_plot_end], sol[4, t_plot_start:t_plot_end],
        sol[1, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[2], fixed_point[4], fixed_point[1], markersize = 15, color = :red)
text!(ax, fixed_point[2], fixed_point[4], fixed_point[1], text = L"O_1", fontsize = labelsize, align = (:center, :top), offset = (0, -23))
display(GLMakie.Screen(), f);

#save(path_to_save * "stable_fixed_point.eps", f)