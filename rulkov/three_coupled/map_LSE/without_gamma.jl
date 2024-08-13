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
include("/home/sergey/work/repo/dynamical-systems/rulkov/Map LSE inheritance symmetrical/main.jl")
include(pathtorepo * "/system.jl");

#---------------------------------------------------------------------
function main()
        
    path_to_save = "/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma/"
    cd(path_to_save)
        
    sys = three_coupled_rulkov;

    params = get_params_three_coupled_rulkov()

    u0 = [-1.209552480291383, -2.974892687599158, -1.2093294078516483, -0.0, -2.9067059214835167,
        -1.209552480291383, -2.974892687599158, -1.2093294078516483, -0.0, -2.9067059214835167,
        1.203416761796254, -2.6976364869798553, 1.0532487761090477, -0.0, -0.0];
    
    length_range = 350;
    range_parameter_1 = range( 0.0, 10.0, length = length_range);
    range_parameter_2 = range(0.0, 10.0, length = length_range);
    
    index_parameter_1 = 10;1
    index_parameter_2 = 11;
    
    name_parameter_1 = "g_1";
    name_parameter_2 = "g_2";

    tspan = (0.0, 10_000);
    t_LSE = 10_000;
    time_setting = (tspan = tspan, t_LSE = t_LSE);

    map_LSE(sys, params, u0,
    range_parameter_1, range_parameter_2, index_parameter_1, index_parameter_2, name_parameter_1, name_parameter_2,
    time_setting; printing = true)

end