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

function d(p1, p2)
    p1x1, p1y1, p1x2, p1y2, p1z = p1
    p2x1, p2y1, p2x2, p2y2, p2z = p2
    dist = sqrt( (p1x1-p2x1)^2 + (p1y1-p2y1)^2 + (p1x2-p2x2)^2 + (p1y2-p2y2)^2 + (p1z-p2z)^2)
    return dist
end

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9();
adaptive = true;
abs_tol = 1e-13;
rel_tol = 1e-13;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[7] = 0.094589;
parameters[8] = 45;

u0_start = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071, 0.006134386010370992];
u0_start = SVector{5}(u0_start);
fixed_point = [-1.01, -0.6367552038434912, -1.01, -0.6367552038434885, -3.421437503947531e-13];

t_end = 50_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_try3, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

distance = zeros(length(sol[1, 1:end]))

for i in range(1, length(distance), step = 1)
    distance[i] = d(sol[i], fixed_point)
end

distance_to_fp = minimum(distance)
println("distance between fp and trajectory: $distance_to_fp");


t_plot_start =  1
t_plot_end = 20_500;

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/fig_7/";

index_x = 1; index_y = 3; index_z = 2;
#= f = Figure(size = (1200,600));
GLMakie.activate!();
ax = LScene(f[1, 1])
lines!(ax, sol[index_x, t_plot_start:t_plot_end], sol[index_y, t_plot_start:t_plot_end],
        sol[index_z, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], markersize = 10, color = :red)
scale!(ax.scene, 1, 1, 1)
display(GLMakie.Screen(), f) =#

#= labelsize = 105;
ticksize = 70; =#
labelsize = 85;
ticksize = 50;
CairoMakie.activate!();
f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 70, ylabeloffset = 60, zlabeloffset = 180,   
    protrusions = (30, 30,120, 30));
xlims!(ax, -1.1, -0.8);
ylims!(ax, -1.2, -0.6);
zlims!(ax, -0.64, -0.60)
lines!(ax, sol[index_x, t_plot_start:t_plot_end], sol[index_y, t_plot_start:t_plot_end],
        sol[index_z, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], markersize = 20, color = :red)
text!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], text = L"O_1", fontsize = labelsize, offset = (0, -120));
display(GLMakie.Screen(), f)

display(f)
save(path_to_save*"spiral_attractor_zoom.eps", f);

#--------------------------------------------------------------
labelsize = 85;
ticksize = 50;
index_x = 1; index_y = 2; index_z = 3;
t_plot_start =  1
t_plot_end = 100_000;
CairoMakie.activate!();
f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"y_1", zlabel = L"x_2",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 100,   
    protrusions = (30, 30,120, 30),
    xticks = [-1.7, 0.0, 1.7], zticks = [-1.7, 0.0, 1.7],
    yticks = [-0.6, 0.0, 0.6]);
lines!(ax, sol[index_x, t_plot_start:t_plot_end], sol[index_y, t_plot_start:t_plot_end],
        sol[index_z, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
scatter!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], markersize = 20, color = :red)
text!(ax, fixed_point[index_x], fixed_point[index_y], fixed_point[index_z], text = L"O_1", fontsize = labelsize,
offset = (-72, -95));

display(f)
save(path_to_save*"spiral_attractor.eps", f);


x_mean = (sol[1, :] + sol[3, :]) / 2;
labelsize = 85;
ticksize = 50;

t_start_plot_timeseries = 1; t_end_plot_timeseries = 200_000;
f = Figure(size = (1200, 700));
ax_x_mean = Axis(f[1, 1], ylabel = L"\overline{x}",
    xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, xticklabelsvisible = false,)
    #xticks = [0, 2000, 4000])

ax_x1_2 = Axis(f[2, 1], xlabel = L"t", ylabel = L"x_i",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
#xticks = [0, 2000, 4000])

lines!(ax_x_mean, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], x_mean[t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :blue)

lines!(ax_x1_2, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], sol[1, t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :red)
lines!(ax_x1_2, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], sol[3, t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :green)
display(f);

save(path_to_save*"x_mean_and_x1_x2_timeseries.pdf", f)