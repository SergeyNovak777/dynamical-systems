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
#---------------------------------------------------------------------
#INCLUDE
include("/home/sergey/work/repo/dynamical-systems/Map LSE src v.2/main.jl")
include(pathtorepo * "/system.jl");

struct Time_setting
    time_calculate_LSE::Int64
    time_attract::Int64
end

#---------------------------------------------------------------------
function main()
        
    path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/g=0.02_change_k1_k_increase_time/"
    cd(path_to_save)
        
    sys = FHN2_try3;

    params = FHN2_try3_params()
    params[3] = 0.02
    u0 = [1.9075684503044907, -0.3904496392396389, 1.9413139567288633, -0.4925421548500994, 0.10209251561045568];
    
    length_range = 300;
    range_parameter_1 = range( 0.0, 0.15, length = length_range);
    range_parameter_2 = range(0.0, 0.3, length = length_range);
    
    index_parameter_1 = 7;
    index_parameter_2 = 8;
    
    name_parameter_1 = "k_1";
    name_parameter_2 = "k_2";
    
    time_attract = 5000
    time_calculate_LSE = 15000
    abstl = 1e-11; reltl = 1e-11;
    time_setting = Time_setting(time_attract, time_calculate_LSE);
    integrator_setting = (alg = DP8(), adaptive = true, abstol = abstl, reltol = reltl, maxiters = 10e6);
    
    type_inheritance = "move to side detect fp"
    
    map_LSE(sys, params, u0,
    range_parameter_1, range_parameter_2, index_parameter_1, index_parameter_2,
    name_parameter_1, name_parameter_2,
    time_setting, integrator_setting, type_inheritance; printing = true, ϵ = 1e-8)

end


