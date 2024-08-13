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
using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie, JLD2

function get_x_y_z(u0)
    x1, y1, z1 = u0[1:3]
    x2, y2, z2 = u0[6:8]
    x3, y3, z3 = u0[11:13]

    return SVector(x1, y1, z1, x2, y2, z2, x3, y3, z3)
end

path_to_save = "/home/sergey/work/repo/dynamical-systems/rulkov/three_coupled/diagram_LSE_change_g1_fix_g2=1/"

params = get_params_three_coupled_rulkov()
#params[11] =  1.0;
tspan = (0, 20_000);
t_LLE = 20_000;

g_start = 10.0; g_end = 0.0; length_range_g = 1001;
g1_range = range(g_start, g_end, length_range_g);

vector_LLE = zeros(length_range_g, 15);
vector_u0s = zeros(length_range_g, 15);

u0_first= [-1.209552480291383, -2.974892687599158, -1.2093294078516483, -0.0, -2.9067059214835167, -1.209552480291383, -2.974892687599158, -1.2093294078516483, -0.0, -2.9067059214835167, 1.203416761796254, -2.6976364869798553, 1.0532487761090477, -0.0, -0.0]
u0_first = get_x_y_z(u0_first)
#[-1.0, 0.8, 1.0, 1.0, -0.8, -1.0, 2.0, 1.0, 0.3];



function first_iter_out_loop(u0_first)

    params[10] =  g1_range[1];
    params[11] =  g1_range[1];
    
    u0_first_iteration_f = three_coupled_rulkov_first_iteration(u0_first, params);
    # SVector{length(u0_first_iteration_f)}(u0_first_iteration_f)
    prob = DiscreteProblem(three_coupled_rulkov, u0_first_iteration_f, tspan, params);
    sol = solve(prob);
    point_from_attractor = sol[:, end];

    ds = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor, params);
    Λs = lyapunovspectrum(ds, t_LLE);
    vector_LLE[1, :] = Λs;
    vector_u0s[1, :] = point_from_attractor;

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
    params[11] = g1_range[index];

    u0 = vector_u0s[index-1, :];

    println("index: $(index)"); flush(stdout);
    println("g: $(params[10])"); flush(stdout);
    println("u0: $(u0)"); flush(stdout);

    u0 = get_x_y_z(u0);
    u0 = three_coupled_rulkov_first_iteration(u0, params);

    prob = DiscreteProblem(three_coupled_rulkov, u0, tspan, params);
    sol = solve(prob);
    point_from_attractor = sol[:, end];

    ds = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor, params);
    Λs = lyapunovspectrum(ds, t_LLE);
    vector_LLE[index, :] = Λs;
    vector_u0s[index, :] = point_from_attractor;

    println("λs: $(vector_LLE[index, 1])")
    println("last point: $(point_from_attractor)"); flush(stdout);
    println("------------------"); flush(stdout);
    println(""); flush(stdout);

end

window_width = 1000; window_height = 350;
f = Figure(size = (window_width, window_height))
ax = Axis(f[1, 1], xlabel = L"g_i", ylabel = L"λ_1",
xticks = [0.0, 2.5, 5.0, 7.5, 10.0], yticks = [-0.1, -0.05, 0.0, 0.05, 0.1])
scatter!(ax, g1_range, vector_LLE[:, 1], markersize = 2.5, color = :blue)
hlines!(ax, 0.0, linewidth = 1.0, color = :black)
xlims!(ax, 0.0, 10.0)
ylims!(ax, -0.1, 0.1)
display(GLMakie.Screen(), f)

#= jldsave(path_to_save * "LLE_g2=1_len_1001.jld2"; vector_LLE)
jldsave(path_to_save * "u0s_g2=1_len_1001.jld2"; vector_u0s) =#