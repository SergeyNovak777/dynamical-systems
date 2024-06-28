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
        
    path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/zoom_chaos_k1_k2_g=0.1/"
    cd(path_to_save)
        
    sys = FHN2_try3;

    params = FHN2_try3_params()
    u0 = [-1.0421826508045808, -0.6270622406846756, -1.0152278251163065, -0.632604603509773, 0.005542362825099361];
    
    length_range = 350;
    range_parameter_1 = range( 0.010, 0.016, length = length_range);
    range_parameter_2 = range(0.9, 1.2, length = length_range);
    
    index_parameter_1 = 7;
    index_parameter_2 = 8;
    
    name_parameter_1 = "k_1";
    name_parameter_2 = "k_2";
    
    time_attract = 3000
    time_calculate_LSE = 10000
    abstl = 1e-11; reltl = 1e-11;
    time_setting = Time_setting(time_attract, time_calculate_LSE);
    integrator_setting = (alg = DP8(), adaptive = true, abstol = abstl, reltol = reltl, maxiters = 10e6);
    
    type_inheritance = "move to side detect fp"
    
    map_LSE(sys, params, u0,
    range_parameter_1, range_parameter_2, index_parameter_1, index_parameter_2,
    name_parameter_1, name_parameter_2,
    time_setting, integrator_setting, type_inheritance; printing = true, ϵ = 1e-8)

end


