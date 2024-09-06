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
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
end
using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie, JLD2, Statistics


Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

params = get_params_three_coupled_rulkov()
params[11] =  6.0;

path_to_save = "/home/sergey/work/repo/dynamical-systems/rulkov/three_coupled/diagram_LSE_change_g1_fix_g2=1/";
path_to_save_image = "/home/sergey/MEGA/dynamical-systems/Rulkov/Images/diagram/";
name_image = "diagram_g_2=6.0.eps"

tspan_Ttr = (0.0, 100_000.0);
tspan_EEs = (0.0, 1_000_000.0);
t_LLE = 20_000;

g_start = 0.0; g_end = 10.0; length_range_g = 1001;
g1_range = range(g_start, g_end, length_range_g);

vector_LLE = zeros(length_range_g, 15);
vector_u0s = zeros(length_range_g, 15);
vector_EEs = zeros(length_range_g);

u0_first= [-1.1587815277382845, -2.962698971867748, -1.1605714358991257, -0.0, -0.0, -0.8984043405608482, -2.9270460934995044, -0.9215955997670924, -0.0, -0.0, -0.9492331950541846, -2.9346719152031735, -0.9633487404593777, -0.0, -0.0]
u0_first = get_x_y_z(u0_first)
#[-1.0, 0.8, 1.0, 1.0, -0.8, -1.0, 2.0, 1.0, 0.3];

function first_iter_out_loop(u0_first)

    params[10] =  g1_range[1];
    u0_first_iteration_f = three_coupled_rulkov_first_iteration(u0_first, params);

    prob = DiscreteProblem(three_coupled_rulkov, SVector{length(u0_first_iteration_f)}(u0_first_iteration_f), tspan_Ttr, params);
    sol = solve(prob);
    point_from_attractor = sol[:, end];
    #-----------------------------------------------------

    ds = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor, params);
    Λs = lyapunovspectrum(ds, t_LLE);
    Λs = sort(Λs, rev = true);
    Λs = Λs
    #-----------------------------------------------------
    if Λs[1] > 0.0 
        prob = DiscreteProblem(three_coupled_rulkov, point_from_attractor, tspan_EEs, params);
        sol = solve(prob);
        xsum = sol[1, :] + sol[6, :] + sol[11, :];
        data = [xsum, sol.t]
        sol = nothing;
        data_local_max = get_local_max(data)
        data_local_min = get_local_min(data)
        drop_artifacts(data_local_max, data_local_min)
        Hs_xsum = Hs(data_local_max[1] ,6);
        count_EE = length(data_local_max[1][ data_local_max[1] .>= Hs_xsum ]);
    else
        count_EE = 0.0;
    end
    #-----------------------------------------------------

    vector_LLE[1, :] = Λs;
    vector_u0s[1, :] = point_from_attractor;
    vector_EEs[1] = count_EE;
    #-----------------------------------------------------

    println("index: $(1)"); flush(stdout);
    println("g: $(params[10])"); flush(stdout);
    println("λs: $(vector_LLE[1, 1])")
    println("u0: $(u0_first)"); flush(stdout);
    println("last point: $(point_from_attractor)"); flush(stdout);
    println("------------------"); flush(stdout);
    println(""); flush(stdout);
end

first_iter_out_loop(u0_first)

for index in 2:length_range_g

    params[10] = g1_range[index];
    u0 = vector_u0s[index-1, :];

    println("index: $(index)"); flush(stdout);
    println("g: $(params[10])"); flush(stdout);
    println("u0: $(u0)"); flush(stdout);

    u0 = get_x_y_z(u0);
    u0 = three_coupled_rulkov_first_iteration(u0, params);

    prob = DiscreteProblem(three_coupled_rulkov, u0, tspan_Ttr, params);
    sol = solve(prob);
    point_from_attractor = sol[:, end];

    ds = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor, params);
    Λs = lyapunovspectrum(ds, t_LLE);
    Λs = sort(Λs, rev = true);
    Λs = Λs
    if Λs[1] > 0.0
        prob = DiscreteProblem(three_coupled_rulkov, point_from_attractor, tspan_EEs, params);
        sol = solve(prob);
        xsum = sol[1, :] + sol[6, :] + sol[11, :];
        data = [xsum, sol.t]
        sol = nothing;
        data_local_max = get_local_max(data)
        data_local_min = get_local_min(data)
        drop_artifacts(data_local_max, data_local_min)
        Hs_xsum = Hs(data_local_max[1] ,6);
        count_EE = length(data_local_max[1][ data_local_max[1] .>= Hs_xsum ]);
    else
        count_EE = 0.0;
    end

    vector_LLE[index, :] = Λs;
    vector_u0s[index, :] = point_from_attractor;
    vector_EEs[index] = count_EE;
   
    println("λs: $(vector_LLE[index, 1])");
    println("EEs: $(vector_EEs[index])"); flush(stdout);
    println("last point: $(vector_u0s[index, :])"); flush(stdout);
    println("------------------"); flush(stdout);
    println(""); flush(stdout);

end

x_y_label_size = 50;
tick_size = 40;
x_ticks = [0, 2, 4, 6, 8];
window_width = 1000; window_height = 600;
CairoMakie.activate!();
f = Figure(size = (window_width, window_height))
ax_LSE = Axis(f[2, 1], xlabel = L"g_1", ylabel = L"λ_1",
        xlabelsize = x_y_label_size, ylabelsize = x_y_label_size,
        xticklabelsize = tick_size, yticklabelsize = tick_size,
        xticks = x_ticks, yticks = [-0.1, 0.0, 0.1],
        xgridvisible = false, ygridvisible = false)
ax_EE = Axis(f[1, 1], xlabel = L"g_1", ylabel = L"count \ EE",
xlabelsize = x_y_label_size, ylabelsize = x_y_label_size-5,
xticklabelsize = tick_size, yticklabelsize = tick_size, xticks = x_ticks, yticks = [0, 400, 800],
xgridvisible = false, ygridvisible = false)

lines!(ax_LSE, g1_range, vector_LLE[:, 1], linewidth = 1.5, color = :red)
lines!(ax_EE, g1_range, vector_EEs, linewidth = 1.5, color = :green)
hlines!(ax_LSE, 0.0, linewidth = 1.0, color = :black)
xlims!(ax_LSE, 0.0, 10.0)
ylims!(ax_LSE, -0.1, 0.1)
xlims!(ax_EE, 0.0, 10.0)
display(GLMakie.Screen(), f)    
save(path_to_save_image*name_image, f);

jldsave(path_to_save * "LLE_g2=6_len_1001.jld2"; vector_LLE)
jldsave(path_to_save * "u0s_g2=6_len_1001.jld2"; vector_u0s)
jldsave(path_to_save * "EEs_g2=6_len_1001.jld2"; vector_EEs)