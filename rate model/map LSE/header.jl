username = "admin";
env = "integrate";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");
include(pathtorepo * "dynamical-systems\\system.jl");

function get_file_name(len_p1, len_p2, name_p1, name_p2)
    
    map_dim = " $(len_p1)x$(len_p2) "
    name = " $(name_p1) $(name_p2) "
    format = ".jld2"
    namefile_LSE = "LSE" * map_dim * name * format
    namefile_u0s = "u0s" * map_dim * name * format
    return namefile_LSE, namefile_u0s

end