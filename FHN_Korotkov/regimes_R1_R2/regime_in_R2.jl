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

using JLD2, CairoMakie, GLMakie
using StaticArrays, DifferentialEquations, DynamicalSystems
using LinearAlgebra
function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

t_truncate(t) = floor(Int64, t / 2)

Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/strip_chaos/LSE_350x350_k_1_k_2.jld2")["λs"]
u0s = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/strip_chaos/u0s_350x350_k_1_k_2.jld2")
init_point = u0s["init_points"]
last_point = u0s["last_points"]

length_range = 350;
k1range = range( 0.090, 0.095, length = length_range);
k2range = range(0.0, 1.0, length = length_range);

index_p1 = 300
index_p2 = 180
println("k1: $(k1range[index_p1])")
println("k2: $(k2range[index_p2])")
println("u0: $(init_point[index_p1, index_p2, :])")
println("λs: $(Λs[index_p1, index_p2, :]) ")
println("last pont: $(last_point[index_p1,index_p2,:])")

alg = Vern9();
adaptive = true;
abs_tol = 1e-12;
rel_tol = 1e-12;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[7] = k1range[index_p1]
parameters[8] = k2range[index_p2]

u0_start = init_point[index_p1, index_p2, :]
u0_start = SVector{5}(u0_start);
fixed_point = [-1.01, -0.6367552038435674, -1.01, -0.6367552038435417, 3.0674543491693373e-13]
jac = Matrix(jac_FHN(fixed_point, parameters, 0));
eigen_values = eigvals(jac);
println("eigen values: $eigen_values");
t_end = 30_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_try3, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_try3, last_point[index_p1,index_p2,:], parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 50000);
println(LSE);

len_sol = length(sol.t);
ttr = t_truncate(len_sol);

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

min_dist = minimum(distance);
println("min distance: $min_dist");

labelsize = 85;
ticksize = 50;
t_plot_start = ttr;
t_plot_end = t_plot_start + 50_000;

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/"

CairoMakie.activate!();
t_start_plot_timeseries = 1; t_end_plot_timeseries = 40_000;
f = Figure(size = (1400, 450));
ax = Axis(f[1, 1], xlabel = L"t", ylabel = L"x_i",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], sol[1, t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :red)
lines!(ax, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], sol[3, t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :green)
#display(GLMakie.Screen(), f)
display(f)
save(path_to_save*"timeseries_R2.pdf", f)

f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 85, ylabeloffset = 85, zlabeloffset = 115,
    protrusions = (30, 30,120, 30),
    xticks = [-1.5, 0.0, 1.5], yticks = [-1.5, 0.0, 1.5], zticks = [-0.5, 0.0, 0.5]);

lines!(ax, sol[1, t_plot_start:t_plot_end], sol[3, t_plot_start:t_plot_end],
        sol[2, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[1], fixed_point[3], fixed_point[2], markersize = 15, color = :red)
text!(ax, fixed_point[1], fixed_point[3], fixed_point[2], text = L"O_1", fontsize = labelsize, align = (:center, :top))
display(GLMakie.Screen(), f);
display(f);
save(path_to_save*"atractor_R2.pdf", f)

index_x = 2; index_y = 4; index_z =5;
f = Figure(size = (1200,600));
ax = LScene(f[1, 1])
lines!(ax, sol[index_x, t_plot_start:t_plot_end], sol[index_y, t_plot_start:t_plot_end],
        sol[index_z, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], markersize = 10, color = :red)
display(GLMakie.Screen(), f)


f = Figure(size = (1200,600));
ax = Axis3(f[1, 1], xlabel = L"y_1", ylabel = L"y_2", zlabel = L"z",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 95, ylabeloffset = 75, zlabeloffset = 115,
    protrusions = (30, 30,120, 30),
    zticks = [-0.04, 0.0, 0.04], xticks = [-0.56, -0.60, -0.64], yticks = [-0.56, -0.60, -0.64]);
xlims!(-0.65, -0.55)
ylims!(-0.65, -0.55)
zlims!(-0.05, 0.05)
lines!(ax, sol[index_x, t_plot_start:t_plot_end], sol[index_y, t_plot_start:t_plot_end],
sol[index_z, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], markersize = 15, color = :red)
text!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], text = L"O_1", fontsize = labelsize, align = (:center, :top), offset = (0, -23))
display(GLMakie.Screen(), f)
display(f)
save(path_to_save*"atractor_zoom_R2.pdf", f)



