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
include("/home/sergey/work/repo/dynamical-systems/Map LSE src/main.jl")
include(pathtorepo * "/system.jl");

struct Time_setting
    time_calculate_LSE::Int64
    time_attract::Int64
    tstep::Float64
end

#---------------------------------------------------------------------
function main()
        
    path_to_save = "/home/sergey/work/repo/dynamical-systems/Tsodyks Markram/Levanova/3 набор параметров/Map LSE/extended_3_map"
    cd(path_to_save)
        
    sys = TM;

    params = TM_model_get_params();
    u0 = [8.254587677389917, 0.7099299111721875, 0.549798431121074];
    
    length_range = 500;
    range_parameter_1 = range( -1.58, -1.78, length = length_range);
    range_parameter_2 = range(0.3, 0.26, length = length_range);
    
    index_parameter_1 = 11;
    index_parameter_2 = 8;
    
    name_parameter_1 = "I_0";
    name_parameter_2 = "U_0";
    
    time_attract = 500
    time_calculate_LSE = 500
    Δt = 1e-3
    maxiters = 1e6
    time_setting = Time_setting(time_attract, time_calculate_LSE, Δt);
    integrator_setting = (alg = RK4(), adaptive = false, dt = Δt, maxiters = maxiters);
    
    map_LSE(sys, params, u0,
    range_parameter_1, range_parameter_2, index_parameter_1, index_parameter_2,
    name_parameter_1, name_parameter_2,
    time_setting, integrator_setting; printing = true)

end