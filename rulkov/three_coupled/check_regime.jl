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
    include("/home/sergey/work/repo/dynamical-systems/system.jl");
end
using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie, JLD2

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function CALCPDF_debil(spikes, threshold, ϵ)
    ee_counter = [sum(i->s<=i<s+ϵ, spikes) for s in threshold]
    pdf = ee_counter ./ length(spikes)
    return pdf
end

path_to_save = "/home/sergey/work/repo/dynamical-systems/rulkov/three_coupled/diagram_LSE_change_g1_fix_g2=1/"

vector_LLE_501 = load(path_to_save * "LLE_g2=1_len_501.jld2")["vector_LLE"]
vector_u0s_501 = load(path_to_save * "u0s_g2=1_len_501.jld2")["vector_u0s"]

vector_LLE_1001 = load(path_to_save * "LLE_g2=1_len_1001.jld2")["vector_LLE"]
vector_u0s_1001 = load(path_to_save * "u0s_g2=1_len_1001.jld2")["vector_u0s"]


g1_range_501 = range(0.0, 10.0, 501);
g1_range_1001 = range(0.0, 10.0, 1001);

params = get_params_three_coupled_rulkov()
params[10] = 0.0;
params[11] = 7.0;
tspan = (0, 300_000);

#u0 = [-0.8653188666707976, -2.9297157423087055, -0.8883590191587176, -0.0, -0.0, -0.8462350817031328, -2.935527209875546, -0.8658876364860313, -0.0, -0.0, -0.9763100224806378, -2.9412769453476204, -0.9837730458388576, -0.0, -0.0]
#u0_first_iteration = three_coupled_rulkov_first_iteration(u0, params);


#[-1.0161473369754137, -2.950229273742396, -1.0154024921006022, -0.0, -0.0, -0.7949290114847254, -2.915437721040376, -0.8384543870469225, -0.0, -0.0, -1.0608144048515458, -2.953899489401151, -1.0589774970895833, -0.0, -0.0]
# SVector{length(u0_first_iteration)}(u0_first_iteration)


u0_first_iteration_atr1 = vector_u0s_1001[end, :]
prob_atr1 = DiscreteProblem(three_coupled_rulkov, u0_first_iteration_atr1, tspan, params);
sol_atr1 = solve(prob_atr1);
Ttr = 200_000;
point_from_attractor_atr1 = sol_atr1[:, end]

x_sum_atr1 = sol_atr1[1, Ttr:end] + sol_atr1[6, Ttr:end] + sol_atr1[11, Ttr:end]
sol_t_atr1 = sol_atr1.t[Ttr:end]

ds_atr1 = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor_atr1, params)
Λs_atr1 = lyapunov(ds_atr1, 50_000)
Λss_atr1 = lyapunovspectrum(ds_atr1, 50_000)
println("Λs atr1: $Λs_atr1");
println("Λss atr1: $Λss_atr1");


u0_first_iteration_atr2 = vector_u0s_501[end, :]
prob_atr2 = DiscreteProblem(three_coupled_rulkov, u0_first_iteration_atr2, tspan, params);
sol_atr2 = solve(prob_atr2);
Ttr = 200_000;
point_from_attractor_atr2 = sol_atr2[:, end]

x_sum_atr2 = sol_atr2[1, Ttr:end] + sol_atr2[6, Ttr:end] + sol_atr2[11, Ttr:end]
sol_t_atr2 = sol_atr2.t[Ttr:end]

ds_atr2 = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor_atr2, params)
Λs_atr2 = lyapunov(ds_atr2, 50_000)
Λss_atr2 = lyapunovspectrum(ds_atr2, 50_000)
println("Λs atr2: $Λs_atr2");
println("Λss atr2: $Λss_atr2");

labelsize = 20;
ticksize = 15;

#= # timeseries of xsum
t_plot_start = 1; t_plot_end = 50_000;
f = Figure(size = (1000, 400));
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_{sum}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(ax, sol_t_atr1[t_plot_start:t_plot_end], x_sum_atr1[t_plot_start:t_plot_end], linewidth = 1.0, color = :blue);
display(GLMakie.Screen(), f); =#


# timeseris x1, x2 ,x3
t_plot_start = 1; t_plot_end = 1_000;
f = Figure(size = (800, 500));
ax1 = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_{1}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
ax2 = Axis(f[2, 1], xlabel = L"time", ylabel = L"x_{2}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
ax3 = Axis(f[3, 1], xlabel = L"time", ylabel = L"x_{3}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
ylims!(ax1, -2, 2)
ylims!(ax2, -2, 2)
ylims!(ax3, -2, 2)
lines!(ax1, sol_atr1.t[Ttr:Ttr + t_plot_end], sol_atr1[1, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :red);
lines!(ax2, sol_atr1.t[Ttr:Ttr + t_plot_end], sol_atr1[6, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :red);
lines!(ax3, sol_atr1.t[Ttr:Ttr + t_plot_end], sol_atr1[11, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :red);

lines!(ax1, sol_atr2.t[Ttr:Ttr + t_plot_end], sol_atr2[1, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :blue);
lines!(ax2, sol_atr2.t[Ttr:Ttr + t_plot_end], sol_atr2[6, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :blue);
lines!(ax3, sol_atr2.t[Ttr:Ttr + t_plot_end], sol_atr2[11, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :blue);
display(GLMakie.Screen(), f);

# timeseris x1, x2, x3 in one axis
#= t_plot_start = 1; t_plot_end = 400;
f = Figure(size = (1400, 500));
ax1 = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_{i}",
xlabelsize = labelsize, ylabelsize = labelsize, 
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);

lines!(ax1, sol.t[Ttr:Ttr + t_plot_end], sol[1, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :red);
lines!(ax1, sol.t[Ttr:Ttr + t_plot_end], sol[6, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :blue);
lines!(ax1, sol.t[Ttr:Ttr + t_plot_end], sol[11, Ttr:Ttr + t_plot_end], linewidth = 1.0, color = :green);
display(GLMakie.Screen(), f); =#

f = Figure(size = (400, 400));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_{2}", zlabel = L"x_3");
scatter!(ax, sol_atr1[1, Ttr:end], sol_atr1[6, Ttr:end], sol_atr1[11, Ttr:end], markersize = 5.0, color = :red);
scatter!(ax, sol_atr2[1, Ttr:end], sol_atr2[6, Ttr:end], sol_atr2[11, Ttr:end], markersize = 5.0, color = :blue);
display(GLMakie.Screen(), f);