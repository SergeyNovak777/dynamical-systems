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
#INCLUDE
include("/home/sergey/work/repo/dynamical-systems/rulkov/Map EE inheritance from matrix/main.jl")
include(pathtorepo * "/system.jl");

#---------------------------------------------------------------------
function main()
        
    path_to_save = "/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/EEs/"
    cd(path_to_save)
    

    matrix_LSE = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/LSE_400x400_g_1_g_2.jld2")["λs"]
    matrix_u0 = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/last_points_400x400_g_1_g_2.jld2")["last_points"]
    matrix_first_point = load("/home/sergey/MEGA/dynamical-systems/Rulkov/data/Maps/without_gamma_broke_first_iter/u0s_400x400_g_1_g_2.jld2")["u0s"]
    sys = three_coupled_rulkov;

    params = get_params_three_coupled_rulkov()

    length_map = 400;
    range_parameter_1 = range( 0.0, 10.0, length = length_map);
    range_parameter_2 = range( 0.0, 10.0, length = length_map);
    
    index_parameter_1 = 10;
    index_parameter_2 = 11;
    
    name_parameter_1 = "g_1";
    name_parameter_2 = "g_2";

    tspan = (0.0, 1_000_000);
    t_LSE = 10_000;
    time_setting = (tspan = tspan, t_LSE = t_LSE);

    map_LSE(sys, params,
    range_parameter_1, range_parameter_2, index_parameter_1, index_parameter_2, name_parameter_1, name_parameter_2,
    time_setting,
    matrix_u0, matrix_first_point, matrix_LSE; printing = true)

end