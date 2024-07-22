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

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9();
adaptive = true;
abs_tol = 1e-12;
rel_tol = 1e-12;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[7] =  0.0932626
parameters[8] = 64.76190476190476

u0_start = [-1.0836728460611933, -0.6318417392022484, -0.9017528537331925, -0.624049721609583, -0.0077920175930263]
u0_start = SVector{5}(u0_start);
fixed_point = [-1.01, -0.6367552038435214, -1.01, -0.6367552038435204, -3.620190802273638e-13];
t_end = 50_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_try3, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_try3, u0_start, parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 100_000);
println(LSE);

len_sol = length(sol.t);
ttr = t_truncate(len_sol);

labelsize = 70;
ticksize = 30;
t_plot_start = ttr;
t_plot_end = t_plot_start + 150_000; #len_sol;

CairoMakie.activate!();

f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"y_1", ylabel = L"y_2", zlabel = L"x_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false);
lines!(ax, sol[2, t_plot_start:t_plot_end], sol[4, t_plot_start:t_plot_end],
        sol[1, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);

display(GLMakie.Screen(), f)
display(f)
#scatter!(ax, fixed_point[2], fixed_point[4], fixed_point[1], markersize = 15, color = :red)
#text!(ax, fixed_point[2], fixed_point[4], fixed_point[1], text = L"O_1", fontsize = labelsize, align = (:center, :top), offset = (0, -23))

pmap = PoincareMap(ds, (1, -1.01))
tr, trange = trajectory(pmap, 700000)

len_tr_map = length(trange);
ttr_map = t_truncate(len_tr_map);

t_plot_start_map = 500000;
t_plot_end_map = t_plot_start_map + 200_000;

f = Figure(size = (900 ,900));
ax = Axis(f[1, 1], xlabel = L"x_2", ylabel = L"y_2", xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false,
    xtickformat = "{:.5f}",
    ytickformat = "{:.7f}");
scatter!(ax, tr[t_plot_start_map:t_plot_end_map, 3], tr[t_plot_start_map:t_plot_end_map, 4], markersize = 1.0, color = :black);
display(GLMakie.Screen(), f);
display(f);