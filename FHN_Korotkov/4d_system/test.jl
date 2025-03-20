username = "irrito"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")
include("/home/irrito/work/repo/dynamical-systems/system.jl")


using StaticArrays, DifferentialEquations, DynamicalSystems
using CairoMakie, GLMakie

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9()
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
adaptive = true;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[7] =  0.0968
parameters[8] = 64.76190476190476

u0_start = [-1.0836728460611933, -0.6318417392022484, -0.9017528537331925, -0.624049721609583];
u0_start = SVector{4}(u0_start);

t_end = 15_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters)

sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_4d, sol[end], parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 5000);
println("LSE: $(LSE)");

labelsize = 50;
ticksize = 25;

length_sol = length(sol);
ttr = t_truncate(length_sol)
t_plot_start =  ttr
t_plot_end = t_plot_start + 150_000

indexx = 1; indexy  = 3; indexz = 2;
f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false)#,
lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);


# Poincare Map
pmap = PoincareMap(ds, (1, -1.01))
tr, trange = trajectory(pmap, 20000)

len_tr_map = length(trange);
ttr_map = t_truncate(len_tr_map);

t_plot_start_map = ttr_map;
t_plot_end_map = t_plot_start_map + 10_000;

f = Figure(size = (900 ,900));
ax = Axis(f[1, 1], xlabel = L"x_2", ylabel = L"y_2", xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false)
scatter!(ax, tr[t_plot_start_map:t_plot_end_map, 3], tr[t_plot_start_map:t_plot_end_map, 4], markersize = 1.0, color = :black);
display(GLMakie.Screen(), f);