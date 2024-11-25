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
    path_to_save = "/home/sergey/MEGA/dynamical-systems/rulkov_2_elements_with_mem_chem/data/map_LSE/"
    cd(path_to_save)
        
    sys = rulkov_two_coupled_chem_mem;

    params = get_params_rulkov_two_coupled_chem_mem()

    params[1] = 3.9; # α
    params[2] = 1.0; # σ
    params[10] = 4.8; # g1
    params[11] = 0.0; # g2
    
    params[12] = 0.1; # k1
    params[13] = 0.0; # k2

    u0 = [-1.953578330045283, -3.991607526888279, -1.9574210901468836,
            -1.97137793066347, -3.819163877171352, -1.9745518175123469,
            -0.11222999003131551];
    
    index_p1 = 11;
    index_p2 = 13;
    
    name_p1 = "g_2";
    name_p2 = "k_2";

    length_map = 200;
    range_p1 = range( 0.0, 10.0, length = length_map);
    range_p2 = range( 0.0, 0.1, length = length_map);
    
    tspan = (0.0, 500_000);
    t_LSE = 250_000;
    time_setting = (tspan = tspan, t_LSE = t_LSE);

    map_LSE(sys, params, u0,
    range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
    time_setting; printing = false, inheritance = "move to side")

end