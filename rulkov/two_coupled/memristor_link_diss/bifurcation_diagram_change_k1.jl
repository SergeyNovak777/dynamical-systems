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
params[1] = 3.9; # α
params[2] = 1.0; # σ
params[10] = 5.0; # g1
params[11] = 1.0; # g2

params[12] = 0.0; # k1
params[13] = 0.000; # k2

index_save_variable = 1;
index_change_parameter = 12;
length_range_change_parameter = 5000;
range_change_parameter = range(0.0, 0.1, length = length_range_change_parameter);

params[index_change_parameter] = range_change_parameter[1]; # k1

ds = DeterministicIteratedMap(rulkov_two_coupled_chem_mem, SVector{7}(u0), params);

amount_saved_points = 2000; # amount saved points for each value parameter
Ttr = 10_000;

output = orbitdiagram(ds, index_save_variable, index_change_parameter, range_change_parameter,
        n = amount_saved_points, Ttr = Ttr);


L = length(range_change_parameter)
x = Vector{Float64}(undef, amount_saved_points*L)
y = copy(x)
for j in 1:L
    x[(1 + (j-1)*amount_saved_points):j*amount_saved_points] .= range_change_parameter[j]
    y[(1 + (j-1)*amount_saved_points):j*amount_saved_points] .= output[j]
end

fig, ax = scatter(x, y; axis = (xlabel = L"k_1", ylabel = L"x_1"),
    markersize = 0.8, color = ("black", 0.05),
)
display(GLMakie.Screen(), fig);
#= array_probs = [deepcopy(prob) for _ in 1:Threads.nthreads()-1];
array_ds = [deepcopy(ds) for _ in 1:Threads.nthreads()-1];
pushfirst!(array_probs, prob);
pushfirst!(array_ds, ds); =#



#= path_to_save = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/";
file_name_LSEs = "g1=$(params[10])_g2=$(params[11])_change_k1_LSE_diagram_LSEs.jld2";
file_name_u0s = "g1=$(params[10])_g2=$(params[11])_change_k1_LSE_diagram_u0s.jld2";
file_name_EEs = "g1=$(params[10])_g2=$(params[11])_change_k1_LSE_diagram_EEs.jld2";


jldsave(path_to_save*file_name_LSEs; array_LSEs);
jldsave(path_to_save*file_name_u0s; array_u0s);
jldsave(path_to_save*file_name_EEs; array_EEs); =#