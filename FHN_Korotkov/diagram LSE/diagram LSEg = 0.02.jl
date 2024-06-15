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
#---------------------------------------------------------------------
# PACKAGES
using DifferentialEquations, DynamicalSystems, StaticArrays, JLD2

function goto_attractor(prob, integrator_setting)

    if integrator_setting.adaptive == true
        point_from_attractor = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, maxiters = integrator_setting.maxiters,
        save_everystep = false, save_start = false)
    else
        point_from_attractor = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        dt = integrator_setting.dt,
        save_everystep = false, save_start = false, maxiters = integrator_setting.maxiters)
    end

    return point_from_attractor[end]
end

function save_in_files(filename_LSEs_vector, filename_u0s_vector, filename_last_points_vector,
    u0, last_point, LSE)
    jldsave(filename_LSEs_vector; u0);
    jldsave(filename_u0s_vector; last_point);
    jldsave(filename_last_points_vector; LSE);
end

function save_in_vectors(index, LSE, u0, last_point)
    LSEs_vector[index, :] = LSE;
    u0s_vector[index, :] = u0;
    last_points_vector[index, :] = last_point;
end

function output(u0, last_point, LSE)
    println("u0: $u0"); flush(stdout)
    println("last_point: $last_point"); flush(stdout)
    println("LSE: $LSE"); flush(stdout)
end

sys = FHN2_try3;
params = FHN2_try3_params()
u0 = [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525]

length_range = 100;
range_parameter_1 = range( 0.001, 0.008, length = length_range);
index_parameter_1 = 7;

k2 = 0.75;
params[8] = k2;

time_attract = 2000.0
tspan = (0.0, time_attract)
time_calculate_LSE = 10000
abstl = 1e-10; reltl = 1e-10;
time_setting = Time_setting(time_attract, time_calculate_LSE);
integrator_setting = (alg = DP8(), adaptive = true, abstol = abstl, reltol = reltl, maxiters = 5e7);

LSEs_vector = zeros(5, length_range);
u0s_vector = zeros(5, length_range);
last_points_vector = zeros(5, length_range);

path_for_saved_files = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/diagram_LSE/g=0.02_fix_k2_change_k1/";
filename_LSEs_vector = "LSEs.jld2";
filename_u0s_vector = "u0s_vector.jld2";
filename_last_points_vector = "last_points.jld2"

for (index, value) in enumerate(range_parameter_1)

    if index == 1
        global u0_cycle = u0;
    end

    params[index_parameter_1] = value;
    prob = ODEProblem(sys, u0_cycle, tspan, params);
    last_point = goto_attractor(prob, integrator_setting);
    ds = CoupledODEs(sys, last_point, params,
        diffeq = integrator_setting);
    LSE = lyapunovspectrum(ds, time_calculate_LSE);

    save_in_vectors(index, LSE, u0_cycle, last_point);
    save_in_files(filename_LSEs_vector, filename_u0s_vector, filename_last_points_vector,
    u0_cycle, last_point, LSE);

    output(u0_cycle, last_point, LSE);

    u0_cycle = last_point;
end

