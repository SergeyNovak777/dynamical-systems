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

function CALCPDF_debil(spikes, threshold, ϵ)
    ee_counter = [sum(i->s<=i<s+ϵ, spikes) for s in threshold]
    pdf = ee_counter ./ length(spikes)
    return pdf
end


#= path_to_load_data = "/home/sergey/MEGA/dynamical-systems/Rulkov/data/PDF/";
filename_sol = "sol__g1=5_g2=4_k1=0.1_k2=0.05.jld2"
filename_PDF = "PDF__g1=5_g2=4_k1=0.1_k2=0.05.jld2"

sol = load(path_to_load_data*filename_sol)["sol"];
EE_mapcopy = load(path_to_load_data*filename_PDF)["EE_mapcopy"]; =#

params = get_params_rulkov_two_coupled_chem_mem()

params[1] = 3.9; # α
params[2] = 1.0; # σ
params[10] = 5.0; # g1
params[11] = 4.0; # g2

params[12] = 0.1; # k1
params[13] = 0.065; # k2

tspan = (0, 4_000_000);
u0 = SVector(-1.953578330045283, -3.991607526888279, -1.9574210901468836,
            -1.97137793066347, -3.819163877171352, -1.9745518175123469,
            -0.11222999003131551);
prob = DiscreteProblem(rulkov_two_coupled_chem_mem, SVector{7}(u0), tspan, params);
sol = solve(prob);
Ttr = 1_000_000; # 1_000_000
point_from_attractor = sol[:, Ttr]
x_sum = (sol[1, :] + sol[4, :])/2

ds = DeterministicIteratedMap(rulkov_two_coupled_chem_mem, point_from_attractor, params)
Λs = lyapunovspectrum(ds, 250_000)
println("Λs: $Λs")

data = [x_sum[Ttr:tspan[2]-Ttr], sol.t[Ttr:tspan[2]-Ttr]]

data_local_max = get_local_max(data)
data_local_min = get_local_min(data)

drop_artifacts(data_local_max, data_local_min)

Hs_xsum = Hs(data_local_max[1] ,6);

println("count EEs: $(length(data_local_max[1][ data_local_max[1] .>= Hs_xsum ]))")

# length(data_local_max[1][ data_local_max[1] .>= Hs_xsum ])

thesholds = range(0.1, 2, 2_500_000);

PDF_old = CALCPDF_debil(data_local_max[1], thesholds, 0.05);

EE_mapcopy = PDF_old;
EE_mapcopy = [ iszero(x) ? NaN : x for x in EE_mapcopy ];

labelsize = 45;
ylabelsize = 50;
ticksize = 35;
CairoMakie.activate!();

path_to_save_timeseries = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/images/";
path_to_save_PDF = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/images/";
path_to_save_data = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/";
filename_sol = "sol__g1=5_g2=4_k1=0.1_k2=0.07.jld2"
filename_PDF = "PDF__g1=5_g2=4_k1=0.1_k2=0.07.jld2"


t_window_plot = 10_000;
xticks_x1_x2 = [-1, 0, 1];
xticks_xsum = [-2, -1, 0, 1, 2];
# timeseries of x1, x2, xsum
f = Figure(size = (1200, 900));

axx1 = Axis(f[1, 1], xlabel = L"time", ylabel = L"x",
xlabelsize = labelsize, ylabelsize = ylabelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
yticks = xticks_x1_x2)

axx2 = Axis(f[2, 1], xlabel = L"time", ylabel = L"x_{2}",
xlabelsize = labelsize, ylabelsize = ylabelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
yticks = xticks_x1_x2)

ax_sum = Axis(f[3, 1], xlabel = L"time", ylabel = L"\bar{x}",
xlabelsize = labelsize, ylabelsize = ylabelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
yticks = xticks_xsum);

lines!(axx1, sol.t[Ttr:Ttr + t_window_plot], sol[1, Ttr:Ttr + t_window_plot], linewidth = 1.0, color = :green)
lines!(axx2, sol.t[Ttr:Ttr + t_window_plot], sol[4, Ttr:Ttr + t_window_plot], linewidth = 1.0, color = :blue)

lines!(ax_sum, data[2][1:t_window_plot], data[1][1:t_window_plot], linewidth = 1.0, color = :black);
hlines!(ax_sum, Hs_xsum, linestyle = :dash, color = :red, linewidth = 5.0);

display(GLMakie.Screen(), f);
save(path_to_save_timeseries*"timeseries_xsum__g1=5_g2=4_k1=0.1_k2=0.065.eps", f)

# pdf old version
f = Figure();
ax = Axis(f[1, 1], xlabel = L"\bar{x}", ylabel = L"PDF", xscale = log10, yscale = log10,
xlabelsize = 60, ylabelsize = 50,
xticklabelsize = 35, yticklabelsize = 35,
xgridvisible = false, ygridvisible = false) #, xticks = [0.12, 0.5, 1.0,  1.5, 2.0, 2.5], xtickformat = "{:.1f}");
lines!(thesholds, EE_mapcopy, linewidth = 1.0, color = :blue);
vlines!(ax, Hs_xsum, linewidth = 5.0, linestyle = :dash, color = :red);
display(GLMakie.Screen(), f);
save(path_to_save_PDF*"PDF__g1=5_g2=4_k1=0.1_k2=0.065.eps", f)

jldsave(path_to_save_data*filename_sol; sol);
jldsave(path_to_save_data*filename_PDF; EE_mapcopy);