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


function get_params_three_coupled_rulkov()
    α = 3.9; σ =1.0; μ = 0.001;
    β_syn = 0.0001; σ_syn = 1.0;
    x_rp = -1.5; x_th = -0.8;
    γ_1 = 0.0; γ_2 = 0.0;
    g1 = 1.0; g2 = 7.0;
    return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2]
end

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function CALCPDF_debil(spikes, threshold, ϵ)
    ee_counter = [sum(i->s<=i<s+ϵ, spikes) for s in threshold]
    pdf = ee_counter ./ length(spikes)
    return pdf
end

params = get_params_three_coupled_rulkov()
params[8] = 0.5;
params[9] = 0.5;
params[10] = 3.0;
params[11] = 2.0;
tspan = (0, 5000000); # 5000000

u0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
u0_first_iteration = three_coupled_rulkov_first_iteration(u0, params);
prob = DiscreteProblem(three_coupled_rulkov, SVector{length(u0_first_iteration)}(u0_first_iteration), tspan, params);
sol = solve(prob);
Ttr = 1_000_000; # 1_000_000
point_from_attractor = sol[:, Ttr]
x_sum = sol[1, :] + sol[6, :] + sol[11, :]

ds = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor, params)
Λs = lyapunovspectrum(ds, 50000)
println("Λs: $Λs")

data = [x_sum[Ttr:end], sol.t[Ttr:end]]
sol = nothing;
x_sum = nothing;

data_local_max = get_local_max(data)
data_local_min = get_local_min(data)

drop_artifacts(data_local_max, data_local_min)

Hs_xsum = Hs(data_local_max[1] ,6);
thesholds = range(0.1, 3, 5_000_000);
PDF_old = CALCPDF_debil(data_local_max[1], thesholds, 0.05);

EE_mapcopy = PDF_old;
EE_mapcopy = [ iszero(x) ? NaN : x for x in EE_mapcopy ];

labelsize = 40;
ticksize = 30;
CairoMakie.activate!();

path_to_save_timeseries = "/home/sergey/MEGA/dynamical-systems/Rulkov/Images/timeseries/";
path_to_save_PDF = "/home/sergey/MEGA/dynamical-systems/Rulkov/Images/PDF/";
path_to_save_data = "/home/sergey/MEGA/dynamical-systems/Rulkov/data/PDF/";

# timeseries of xsum
t_plot_start = 1; t_plot_end = 200_000;
f = Figure(size = (1000, 400));
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_{sum}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(ax, data[2][t_plot_start:t_plot_end], data[1][t_plot_start:t_plot_end], linewidth = 1.0, color = :blue);
hlines!(ax, Hs_xsum, linestyle = :dash, color = :red, linewidth = 3.0);
display(GLMakie.Screen(), f);
#save(path_to_save_timeseries*"timeseries__xsum_g1=3.0_g2=2.0_with_gamma.eps", f)

# pdf old version
f = Figure();
ax = Axis(f[1, 1], xlabel = L"x_{sum}", ylabel = L"PDF",
xscale = log10, yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(thesholds, EE_mapcopy, linewidth = 1.0, color = :blue);
vlines!(ax, Hs_xsum, linewidth = 3.0, linestyle = :dash, color = :red);
display(GLMakie.Screen(), f);
#save(path_to_save_PDF*"PDF_g1=3.0_g2=2.0_with_gamma.eps", f)