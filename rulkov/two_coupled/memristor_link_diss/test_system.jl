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
t_window_plot = t_tr + 50_000;

params = get_params_rulkov_two_coupled_chem_mem()

params[1] = 3.9; # α
params[2] = 1.0; # σ
params[10] = 5.0; # g1
params[11] = 4.0; # g2

params[12] = 0.1; # k1
params[13] = 0.052; # k2

u0 = SVector(-1.953578330045283, -3.991607526888279, -1.9574210901468836, -1.46916966508604, -3.5865593516810397, -1.4738167754245777, -0.4836043147223059);

prob = DiscreteProblem(rulkov_two_coupled_chem_mem, SVector{7}(u0), tspan, params);
sol = solve(prob);

x_sum = sol[1, t_tr:end] + sol[4, t_tr:end]
t_range = sol.t[t_tr:end]

ds = DeterministicIteratedMap(rulkov_two_coupled_chem_mem, sol[end], params)
Λs = lyapunovspectrum(ds, 1_000_000)
println("LSE: $Λs");

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

data = [x_sum, sol.t]
data_local_max = get_local_max(data)
data_local_min = get_local_min(data)

drop_artifacts(data_local_max, data_local_min)
Hs_xsum = Hs(data_local_max[1], 8);

print("count EEs: $(count(data_local_max[1].>=Hs_xsum))");

CairoMakie.activate!();
f = Figure(size = (400, 400))
ax = Axis3(f[1, 1])
scatter!(ax, sol[1, t_tr:t_window_plot], sol[4, t_tr:t_window_plot], sol[2, t_tr:t_window_plot], markersize = 1.0, color = :black)
display(GLMakie.Screen(), f);

f = Figure(size = (1000, 400))
axx1 = Axis(f[1, 1])
axx2 = Axis(f[2, 1])
axxsum = Axis(f[3, 1])
lines!(axx1, sol.t[t_tr:t_window_plot], sol[1, t_tr:t_window_plot], linewidth = 1.0, color = :green)
lines!(axx2, sol.t[t_tr:t_window_plot], sol[4, t_tr:t_window_plot], linewidth = 1.0, color = :blue)
lines!(axxsum, t_range[1:t_window_plot], x_sum[1:t_window_plot], linewidth = 1.0, color = :black)
#scatter!(axxsum, data_local_max[2], data_local_max[1], markersize = 2.5, color = :red)
#xlims!(axxsum, t_range[1], t_range[t_window_plot])
hlines!(axxsum, Hs_xsum, linestyle = :dash, color = :red, linewidth = 3.0);
display(GLMakie.Screen(), f);