username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie
using JLD2;
include("/home/sergey/work/repo/dynamical-systems/system.jl");
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl");
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl");

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

t_sol = (0, 5_000_000);
t_tr = 2_000_00;
t_calc_LSE = 500_000;

u0 = SVector(-2.083638390440308, -3.9302148554862937, -2.0867818075436624,
            -1.97137793066347, -3.819163877171352, -1.9745518175123469,
            -0.11222999003131551)

params = get_params_rulkov_two_coupled_chem_mem()
params[10] = 5.0; # g1
params[11] = 2.0; # g2

params[12] = 0.0; # k1
params[13] = 0.000; # k2

length_range_k1 = 1_000;
index_control_parameter = 12;
range_k1 = range(0.0, 0.3, length = length_range_k1);

array_LSEs = zeros(length_range_k1, length(u0));
array_u0s = zeros(length_range_k1, length(u0)); 
array_EEs = zeros(length_range_k1);

params[index_control_parameter] = range_k1[1]; # k1

prob = DiscreteProblem(rulkov_two_coupled_chem_mem, SVector{7}(u0), t_sol, params);
ds = DeterministicIteratedMap(rulkov_two_coupled_chem_mem, SVector{7}(u0), params);

#= array_probs = [deepcopy(prob) for _ in 1:Threads.nthreads()-1];
array_ds = [deepcopy(ds) for _ in 1:Threads.nthreads()-1];
pushfirst!(array_probs, prob);
pushfirst!(array_ds, ds); =#

# first iteration
sol = solve(prob);
Λs = lyapunovspectrum(ds, t_calc_LSE, u0 = sol[end])
Λs = sort(Λs, rev = true);

array_u0s[1, :] = sol[end];
array_LSEs[1, :] = Λs;

if Λs[1] >= 0.001
    x_sum = sol[1, t_tr:end] + sol[4, t_tr:end];
    data = [x_sum, sol.t];

    data_local_max = get_local_max(data);
    data_local_min = get_local_min(data);

    drop_artifacts(data_local_max, data_local_min)
    Hs_xsum = Hs(data_local_max[1] ,6);

    count_EEs = count(data_local_max[1].>=Hs_xsum)
    array_EEs[1] = count_EEs;
end

println("control_param: $(params[index_control_parameter])");
println("u0: $u0");
println("last point: $(array_u0s[1, :])");
println("Λs: $(array_LSEs[1, :])");
println("count EEs: $(array_EEs[1])");
println("---------------------------------");
println("");

for index_cycle in range(2, length_range_k1, step = 1)

    newp = copy(prob.p)
    newp[index_control_parameter] = range_k1[index_cycle];
    u0_local = array_u0s[index_cycle - 1 , :];

    set_parameter!(ds, index_control_parameter, newp[index_control_parameter]);
    prob_local = remake(prob, u0 = u0_local, p = newp);
    sol_local = solve(prob_local);

    Λs_local = lyapunovspectrum(ds, t_calc_LSE, u0 = u0_local);
    Λs_local = sort(Λs_local, rev = true);

    if Λs_local[1] >= 0.001
        println("detect EEs");
        x_sum_local = sol_local[1, t_tr:end] + sol_local[4, t_tr:end];
        data_local = [x_sum_local, sol_local.t];

        data_local_max_lc = get_local_max(data_local);
        data_local_min_lc = get_local_min(data_local);

        drop_artifacts(data_local_max_lc, data_local_min_lc)
        Hs_xsum_lc = Hs(data_local_max_lc[1] ,6);

        count_EEs_lc = count(data_local_max_lc[1].>=Hs_xsum_lc)
        array_EEs[index_cycle] = count_EEs_lc;
    end
    
    array_u0s[index_cycle, :] = sol_local[end];
    array_LSEs[index_cycle, :] = Λs_local;

    println("index cycle: $(index_cycle)");
    println("control param prob: $(prob_local.p[index_control_parameter])");
    println("control param ds: $(ds.p[index_control_parameter])");
    println("u0: $u0_local");
    println("last point: $(array_u0s[index_cycle, :])")
    println("Λs: $(array_LSEs[index_cycle, :])");
    println("count EEs: $(array_EEs[index_cycle])");
    println("---------------------------------");
    println("");
end

path_to_save = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/";
file_name_LSEs = "g1=$(params[10])_g2=$(params[11])_change_k1_LSE_diagram_LSEs.jld2";
file_name_u0s = "g1=$(params[10])_g2=$(params[11])_change_k1_LSE_diagram_u0s.jld2";
file_name_EEs = "g1=$(params[10])_g2=$(params[11])_change_k1_LSE_diagram_EEs.jld2";


jldsave(path_to_save*file_name_LSEs; array_LSEs);
jldsave(path_to_save*file_name_u0s; array_u0s);
jldsave(path_to_save*file_name_EEs; array_EEs);



window_height, window_width = 400, 1000;
xlabel = L"k_1";
ylabel_LSEs, ylabelEEs = L"LSE", L"EE_{count}"
label_size = 35;
tickssize = 25;
linewidth = 2.0;
color_LLE1 = :red;
color_LLE2 = :green;
color_LLE3 = :blue;
color_LLE4 = :black;


CairoMakie.activate!();
figure = Figure(figsize = (window_height, window_width));
ax_LSE = Axis(figure[1,1],
        xlabel = xlabel, ylabel = ylabel_LSEs,
        xlabelsize = label_size, ylabelsize = label_size,
        xticklabelsize = tickssize, yticklabelsize = tickssize,
);
ax_EE = Axis(figure[2, 1],
        xlabel = xlabel, ylabel = ylabelEEs,
        xlabelsize = label_size, ylabelsize = label_size,
        xticklabelsize = tickssize, yticklabelsize = tickssize,
);

lines!(ax_LSE, range_k1, array_LSEs[:, 1], linewidth = linewidth, color = color_LLE1);
lines!(ax_LSE, range_k1, array_LSEs[:, 2], linewidth = linewidth, color = color_LLE2);
lines!(ax_LSE, range_k1, array_LSEs[:, 3], linewidth = linewidth, color = color_LLE3);

lines!(ax_EE, range_k1, array_EEs, linewidth = linewidth, color = color_LLE1);
display(GLMakie.Screen(), figure);