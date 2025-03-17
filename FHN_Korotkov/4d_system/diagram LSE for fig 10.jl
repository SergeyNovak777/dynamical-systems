if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "irrito"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
    include("/home/irrito/work/repo/dynamical-systems/system.jl")
end

using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie, JLD2

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

u0 = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];
params = FHN2_try3_params()

params[7] = 0.09;

integ_set = (alg = Vern9(), adaptive = true, abstol=1e-13, reltol=1e-13, maxiters = 1e8)

t_span_Ttr = (0.0, 1000.0)
t_calc_LSE = 10_000;

k2_start = 0.0
k2_end = 100.0
len_range = 1000
rangek2 = range(k2_start, k2_end, length = len_range)

matrix_u0s = zeros(len_range, length(u0));
matrix_LSE = zeros(len_range, length(u0));

matrix_u0s[1, :] = u0;

for (index, value_parameter) in enumerate(rangek2)

    println("index: $(index)"); flush(stdout)
    println("value parameter : $(value_parameter)"); flush(stdout);

    u0_local = matrix_u0s[index, :]; # get u0
    params[8] = value_parameter;

    prob = ODEProblem(FHN2_4d, SVector{4}(u0_local), t_span_Ttr, params)
    sol = solve(prob, integ_set.alg, adaptive = integ_set.adaptive,
                abstol = integ_set.abstol, reltol = integ_set.reltol, 
                maxiters = integ_set.maxiters,
                save_everystep = false, save_start = false);
    
    last_point = sol[end];

    sol = nothing;
    u0_local = nothing;

    ds = CoupledODEs(FHN2_4d, last_point, params, diffeq = integ_set);
    LSE = lyapunovspectrum(ds, t_calc_LSE);

    println("LSE: $(LSE);"); flush(stdout);
    println("-----------------------"); flush(stdout);
    matrix_LSE[index, :] = LSE;

    if index < len_range
        matrix_u0s[index+1, :] = last_point;
    end
end

file_name_LSE = "matrix_LSE_fix_k1=0.09_change_k2_from_0_to_100.jld2";
file_name_u0s = "matrix_u0s_fix_k1=0.09_change_k2_from_0_to_100.jld2";

path_to_save = "/home/irrito/MEGA/dynamical_systems/FHN/data/";
jldsave(path_to_save*file_name_LSE; matrix_LSE);
jldsave(path_to_save*file_name_u0s; matrix_u0s);


test_load_LSE = load(path_to_save*file_name_LSE)["matrix_LSE"]

markersize = 1.5;
lbsize = 50;
ticksize = 35;
CairoMakie.activate!();
fig = Figure(size = (1200, 350))
axis = Axis(fig[1,1],
        xlabel = L"k_2",  ylabel = L"x_1",
        xlabelsize = lbsize, ylabelsize = lbsize,
        xticklabelsize = ticksize,yticklabelsize = ticksize,
        xgridvisible = false, ygridvisible = false);

lines!(axis, rangek2, test_load_LSE[:, 1], linewidth = 1.5, color = :red);
lines!(axis, rangek2, test_load_LSE[:, 2], linewidth = 1.5, color = :blue);
hlines!(axis, 0, rangek2, linestyle = :dash, linewidth = 1.5, color = :black)
display(GLMakie.Screen(), fig)

path_to_save = "/home/irrito/MEGA/dynamical_systems/FHN/images/"

file_name_image = "diagram_LSE_fig_10.eps";
full_path = path_to_save * file_name_image;

save(full_path, fig)