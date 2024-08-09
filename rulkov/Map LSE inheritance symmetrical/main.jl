#---------------------------------------
#INCLUDE
include(joinpath(@__DIR__, "header.jl"))

function map_LSE(sys, params, u0,
    range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
    time_setting, integrator_setting; printing = false)

    #---------------------------------------
    len_p1 = length(range_p1)
    len_p2 = length(range_p2)
    dim = length(u0)

    #---------------------------------------
    global λs = zeros(len_p1, len_p2, dim)
    global u0s = zeros(len_p1, len_p2, dim)
    global last_points = zeros(len_p1, len_p2, dim)
    #---------------------------------------
    namefile_LSE, namefile_u0s = get_file_name(len_p1, len_p2,
     name_p1, name_p2)
    
    #---------------------------------------

    
end