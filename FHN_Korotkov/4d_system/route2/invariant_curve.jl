username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")
include(pathtorepo * "/system.jl")

using StaticArrays, SciMLBase, OrdinaryDiffEqVerner, ChaosTools, CairoMakie, GLMakie

t_truncate(t) = floor(Int64, t / 2)
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
parameters[7] =  0.09325
parameters[8] = 64.76190476190476

u0_start = [-0.9753300387911031, -0.6357434525860388, -1.0780505297640295, -0.6189964286441518];

u0_start = SVector{4}(u0_start);

t_end = 2_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_4d, sol[end], parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 50_000);
println("LSE: $(LSE)")

x1 = x2 = interval(-1.5, 1.5)
y1 = y2 = interval(-1.5, 1.5)
box = [x1, y1, x2, y2]

fixed_point, eigs, _ = fixedpoints(ds, box)
fixed_point = fixed_point[1];
labelsize = 75 #85;
ticksize = 30 #50;

t_plot_start = t_truncate(length(sol.t));
t_plot_end = t_plot_start + 15_000; #len_sol;

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/scenario/"

CairoMakie.activate!();

indexx = 1; indexy = 3; indexz = 2;

f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"y_1", ylabel = L"y_2", zlabel = L"x_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 85, ylabeloffset = 85, zlabeloffset = 115,
    protrusions = (30, 30,120, 30))#,
    #xticks = [-0.635, -0.622], yticks = [-0.635, -0.622], zticks = [-1.05, -0.95]);
lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[indexx], fixed_point[indexy], fixed_point[indexz], markersize = 15, color = :red)
text!(ax, fixed_point[indexx], fixed_point[indexy], fixed_point[indexz], text = L"O_1", fontsize = labelsize, align = (:center, :top), offset = (0, -23))
display(GLMakie.Screen(), f);

pmap = PoincareMap(ds, (1,  -1.013))
tr, trange = trajectory(pmap, 150_000)

len_tr_map = length(trange);
ttr_map = t_truncate(len_tr_map);

t_plot_start_map = ttr_map;
t_plot_end_map = t_plot_start_map + 75_000;

f = Figure(size = (700 ,700));
ax = Axis(f[1, 1], xlabel = L"x_2", ylabel = L"y_2", xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false,
    xticks = [-1.041, -1.038, -1.035]);
scatter!(ax, tr[t_plot_start_map:t_plot_end_map, 3], tr[t_plot_start_map:t_plot_end_map, 4], markersize = 1.0, color = :black);
display(GLMakie.Screen(), f);

filename = "route1_a_k1=$(parameters[7]).eps";
path_to_save = "/home/sergey/work/images/FHN/route2/"
#"/home/sergey/MEGA2/dynamical_systems/FHN/images/route1/";
save(path_to_save * filename, f)