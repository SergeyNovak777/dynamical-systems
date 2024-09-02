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
end
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl")
include("/home/sergey/work/repo/dynamical-systems/system.jl")
using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie, JLD2

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

params = get_params_three_coupled_rulkov()

tspan = (0, 5_000_000); # 5000000

length_map = 400;
range_parameter_1 = range( 0.0, 10.0, length = length_map);
range_parameter_2 = range( 0.0, 10.0, length = length_map);

matrix_LSE = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/LSE_400x400_g_1_g_2.jld2")["λs"]
matrix_last_point = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/last_points_400x400_g_1_g_2.jld2")["last_points"]
matrix_first_point = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/u0s_400x400_g_1_g_2.jld2")["u0s"]
matrix_EEs = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/EEs/EEs_400x400_g_1_g_2.jld2")["matrix_EEs"]

index_p1 = 150;
index_p2 = 250;

println("g_1: $(range_parameter_1[index_p1])");
println("g_2: $(range_parameter_2[index_p2])");
println("first point: $(matrix_first_point[index_p1, index_p2, :])");
println("last point: $(matrix_last_point[index_p1, index_p2, :])");
println("LSE: $(matrix_LSE[index_p1, index_p2, :])");
println("EEs: $(matrix_EEs[index_p1, index_p2])");

length_range = 400;
g_1_range = range( 0.0, 10.0, length = length_range);
g_2_range = range( 0.0, 10.0, length = length_range);
absmax = maximum(abs.(matrix_EEs[:, :]))
mn, mx =  -100, 100
CairoMakie.activate!()  
f = Figure()    
ticksize = 35
tickpad = 10.0
textsize = 25
textsizecurve = 30
lw = 3.5
mksize = 12
ax = Axis(f[1, 1], xlabel = L"g_2",ylabel = L"g_1", xlabelsize = 50, ylabelsize = 50,
            xticklabelsize = ticksize, yticklabelsize = ticksize,
            xgridvisible  = false, ygridvisible = false,
            xticklabelpad = tickpad, yticklabelpad = tickpad)#,
hm = heatmap!(ax, g_1_range, g_2_range, matrix_EEs, colormap = :seismic,
            colorrange = (mn, mx))
scatter!(ax, range_parameter_1[index_p1], range_parameter_2[index_p2], markersize = 5, color = :black)
display(GLMakie.Screen(), f);



params[10] = range_parameter_1[index_p1];
params[11] = range_parameter_2[index_p2];

u0 =   matrix_last_point[index_p1, index_p2, :]
prob = DiscreteProblem(three_coupled_rulkov, u0, tspan, params);
sol = solve(prob);


Ttr = 1;
point_from_attractor = sol[:, Ttr]
x_sum = sol[1, :] + sol[6, :] + sol[11, :]
data = [x_sum[Ttr:end], sol.t[Ttr:end]]

data_local_max = get_local_max(data)
data_local_min = get_local_min(data)

drop_artifacts(data_local_max, data_local_min)

Hs_xsum = Hs(data_local_max[1] ,6);

ds = DeterministicIteratedMap(three_coupled_rulkov, u0, params)
Λs = lyapunovspectrum(ds, 10_000)
println("Λs: $Λs")

labelsize = 40;
ticksize = 30;

# timeseries of xsum
t_plot_start = 1; t_plot_end = 3_000_000;
f = Figure(size = (1000, 400));
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_{sum}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(ax, data[2][t_plot_start:t_plot_end], data[1][t_plot_start:t_plot_end], linewidth = 1.0, color = :blue);
hlines!(ax, Hs_xsum, linestyle = :dash, color = :red, linewidth = 3.0);
display(GLMakie.Screen(), f);

f = Figure(size = (400, 400));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_{2}", zlabel = L"x_3");
scatter!(ax, sol[1, Ttr:end], sol[6, Ttr:end], sol[11, Ttr:end], markersize = 1.0, color = :blue);
display(GLMakie.Screen(), f);