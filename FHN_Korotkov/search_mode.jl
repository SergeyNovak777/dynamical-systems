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

function conv_sol_u(sol_u)
    len_sol = length(sol_u);
    dim = length(sol_u[1]);
    sol = zeros(len_sol, dim);

    for index in range(1, dim)
        sol[:, index] = [v[index] for v in sol_u];
    end

    return sol;
end

function get_sol(prob, integrator_setting)

    sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol,
                maxiters = integrator_setting.maxiters);
    sol_t = sol.t;
    sol_u = conv_sol_u(sol.u);
    return sol_t, sol_u;
end

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9();
adaptive = true;
abs_tol = 1e-12;
rel_tol = 1e-12;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[7] =  0.09203007518796992
parameters[8] = 64.76190476190476

u0_start = [-1.0073393282360215, -0.6392350435710693, -1.022720802290101, -0.6276963254782901, -0.011538718093154163]
u0_start = SVector{5}(u0_start);
fixed_point = [-1.01, -0.6367552038435214, -1.01, -0.6367552038435204, -3.620121549323092e-13]
t_end = 50_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_try3, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_try3, u0_start, parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 10000);
println(LSE);

len_sol = length(sol.t);
ttr = t_truncate(len_sol);

labelsize = 85;
ticksize = 50;
t_plot_start = ttr;
t_plot_end = t_plot_start + 50_000; #len_sol;

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/scenario/"

CairoMakie.activate!();

f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"y_1", ylabel = L"y_2", zlabel = L"x_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xticks = [-0.63, -0.6], yticks = [-0.63, -0.6], zticks = [-1.07, -0.8],
    xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 115,
    protrusions = (30, 30,120, 30),
    elevation = 0.03pi);
lines!(ax, sol[2, t_plot_start:t_plot_end], sol[4, t_plot_start:t_plot_end],
        sol[1, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[2], fixed_point[4], fixed_point[1], markersize = 15, color = :red)
text!(ax, fixed_point[2], fixed_point[4], fixed_point[1], text = L"O_1", fontsize = labelsize, align = (:center, :top), offset = (0, -50))
display(f);

save(path_to_save * "narrowing.eps", f)


#= pmap = PoincareMap(ds, (1, -1.01))
tr, trange = trajectory(pmap, 700000)

len_tr_map = length(trange);
ttr_map = t_truncate(len_tr_map);

t_plot_start_map = 500000;
t_plot_end_map = t_plot_start_map + 80_000;

f = Figure(size = (900 ,900));
ax = Axis(f[1, 1], xlabel = L"x_2", ylabel = L"y_2", xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false,
    xtickformat = "{:.1f}",
    xticks = [-1.1, -1.0, -0.9, -0.81]);
scatter!(ax, tr[t_plot_start_map:t_plot_end_map, 3], tr[t_plot_start_map:t_plot_end_map, 4], markersize = 2.0, color = :black);
display(f);
save(path_to_save * "narrowing_poincare.eps", f) =#