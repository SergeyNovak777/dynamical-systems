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
    include(pathtorepo * "/system.jl");
end

include("/home/sergey/work/repo/dynamical-systems/rate model/Map LSE/main.jl")

using DifferentialEquations, DynamicalSystems, StaticArrays, JLD2


struct Time_setting

    time_calculate_LSE::Int64
    time_attract::Int64
    tstep::Float64

end

cd("/home/sergey/work/repo/dynamical-systems/rate model/Map LSE/u0_zero/")

sys = rate_model;
params = rate_model_get_params();
values, index_params = rate_model_help(params);
u0 = zeros(5); #ones(5) * 0.25;

length_range = 100;
range_parameter_1 = range( 0.0, 10.0, length = length_range);
range_parameter_2 = range(0.0, 2.0, length = length_range);

index_parameter_1 = 21;
index_parameter_2 = 9;

name_parameter_1 = "γY";
name_parameter_2 = "IE";

time_attract = 1000
time_calculate_LSE = 1000
Δt = 1e-3
time_setting = Time_setting(time_attract, time_calculate_LSE, Δt);
integrator_setting = (alg = RK4(), adaptive = false, dt = Δt, maxiters = 5.0e7);

map_LSE(sys, params, u0,
range_parameter_1, range_parameter_2, index_parameter_1, index_parameter_2,
name_parameter_1, name_parameter_2,
time_setting, integrator_setting; printing = true)