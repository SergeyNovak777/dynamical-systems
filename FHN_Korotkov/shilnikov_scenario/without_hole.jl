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
using LinearAlgebra

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
parameters[7] =  0.091695
parameters[8] = 64.76190476190476

fixed_point = [-1.01, -0.636755203843491, -1.01, -0.6367552038434884, -3.4203139675502085e-13];
jac = Matrix(jac_FHN(fixed_point, parameters, 0.0));
eigen_vals = eigvals(jac);


u0_start = [-1.0221779287339454, -0.6359770136715965, -1.0119938655496694, -0.6322079892358413, -0.003769024436094622]
u0_start = SVector{5}(u0_start);

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
t_plot_end = t_plot_start + 100_000; #len_sol;

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/scenario/"

CairoMakie.activate!();
indexx = 1; indexy  = 2; indexz = 3;
f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"y_1", zlabel = L"x_2",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 115,   
    protrusions = (30, 30,120, 30))#,
    #xticks = [-0.5, 0.5], yticks = [-0.5, 0.5], zticks = [-1.5, -1.0],
    #elevation = 0.04pi);
lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[indexx], fixed_point[indexy], fixed_point[indexz], markersize = 15, color = :red)
text!(ax, fixed_point[2], fixed_point[4], fixed_point[1], text = L"O_1", fontsize = labelsize, align = (:center, :top), offset = (0, -30))
display(GLMakie.Screen(), f);

#save(path_to_save * "homoclinic.eps", f)

function d(p1, p2)
    p1x1, p1y1, p1x2, p1y2, p1z = p1
    p2x1, p2y1, p2x2, p2y2, p2z = p2
    dist = sqrt( (p1x1-p2x1)^2 + (p1y1-p2y1)^2 + (p1x2-p2x2)^2 + (p1y2-p2y2)^2 + (p1z-p2z)^2)
    return dist
end

distance = zeros(length(sol[1, ttr:end]))

for i in range(1, length(distance), step = 1)
    distance[i] = d(sol[i], fixed_point)
end

minimum(distance)

pmap = PoincareMap(ds, (1, -1.01))
tr, trange = trajectory(pmap, 700000)

len_tr_map = length(trange);
ttr_map = t_truncate(len_tr_map);

t_plot_start_map = 500000;
t_plot_end_map = t_plot_start_map + 80_000;

f = Figure(size = (900 ,900));
ax = Axis(f[1, 1], xlabel = L"x_2", ylabel = L"y_2", xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false);
scatter!(ax, tr[t_plot_start_map:t_plot_end_map, 3], tr[t_plot_start_map:t_plot_end_map, 4], markersize = 2.0, color = :black);
display(f);
#save(path_to_save * "homoclinic_poincare.eps", f)