username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

include("/home/sergey/work/repo/dynamical-systems/system.jl");
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl");
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl");

tspan = (0, 2_500_000);
t_tr = 1_000_000;
t_window_plot = t_tr + 1_500_000;

params = get_params_rulkov_two_coupled_chem_mem()

params[1] = 3.9; # α
params[2] = 1.0; # σ
params[10] = 4.5; # g1
params[11] = 4.5; # g2

params[12] = 0.0; # k1
params[13] = 0.0; # k2

u0 = SVector(-1.953578330045283, -3.991607526888279, -1.9574210901468836,
            -1.97137793066347, -3.819163877171352, -1.9745518175123469,
            -0.11222999003131551);
    #= -1.953578330045283, -3.991607526888279, -1.9574210901468836,
 -1.46916966508604, -3.5865593516810397, -1.4738167754245777,
  -0.4836043147223059); =#

#= 1.0, 0.3, 0.01,
            1.5, 0.4, 0.5,
            1.0-1.5 =#

prob = DiscreteProblem(rulkov_two_coupled_chem_mem, SVector{7}(u0), tspan, params);
sol = solve(prob);

x_sum = sol[1, t_tr:tspan[2]] + sol[4, t_tr:tspan[2]]
t_range = sol.t[t_tr:tspan[2]]

ds = DeterministicIteratedMap(rulkov_two_coupled_chem_mem, sol[end], params)
Λs = lyapunovspectrum(ds, 200_000)
println("LSE: $Λs");

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

data = [x_sum, sol.t]
data_local_max = get_local_max(data)
data_local_min = get_local_min(data)

drop_artifacts(data_local_max, data_local_min)
Hs_xsum = Hs(data_local_max[1] ,6);

print("count EEs: $(count(data_local_max[1].>=Hs_xsum))");

index_x = 1;
index_y = 2;
index_z = 3;

CairoMakie.activate!();
f = Figure(size = (900, 600))
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_1", xlabelsize = 50, ylabelsize = 50, zlabelsize = 50)
scatter!(ax, sol[1, t_tr:t_window_plot], sol[4, t_tr:t_window_plot], sol[2, t_tr:t_window_plot],
markersize = 2.0, color = :black)
display(GLMakie.Screen(), f);

#= f = Figure(size = (400, 400))
ax = Axis(f[1, 1])
scatter!(ax, sol[index_x, t_tr:t_window_plot], sol[index_y, t_tr:t_window_plot], markersize = 5.0, color = :green)
scatter!(ax, sol[4, t_tr:t_window_plot], sol[5, t_tr:t_window_plot], markersize = 5.0, color = :blue)
display(GLMakie.Screen(), f); =#

f = Figure(size = (1000, 700))
ax1 = Axis(f[1, 1])
ax2 = Axis(f[2, 1])
axsum = Axis(f[3, 1])
t_start = t_tr;
t_end = t_start + 1_000;
lines!(ax1, sol.t[t_start:t_end], sol[1, t_start:t_end], linewidth = 1.0, color = :green)
lines!(ax2, sol.t[t_start:t_end], sol[4, t_start:t_end], linewidth = 1.0, color = :blue)
lines!(axsum, t_range[1:end], x_sum[1:end], linewidth = 1.0, color = :black)

display(GLMakie.Screen(), f);