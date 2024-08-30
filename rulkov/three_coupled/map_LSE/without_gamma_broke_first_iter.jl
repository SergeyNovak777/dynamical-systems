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
        
    path_to_save = "/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/"
    cd(path_to_save)
        
    sys = three_coupled_rulkov;

    params = get_params_three_coupled_rulkov()

    u0 = [-0.5217736139644187, -2.4707736139644187, -1.0, -0.0, -0.0,
            -1.0, -2.4717736139644186, 1.4297481596495456, -0.0, -0.0, 1.4297481596495456,
            -2.470343865804769, 0.09202545431485776, -0.0, -0.0]
    
    length_map = 400;
    range_parameter_1 = range( 0.0, 10.0, length = length_map);
    range_parameter_2 = range( 0.0, 10.0, length = length_map);
    
    index_parameter_1 = 10;
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