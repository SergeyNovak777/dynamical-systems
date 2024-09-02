#---------------------------------------
#INCLUDE
include(joinpath(@__DIR__, "header.jl"))
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")

using DifferentialEquations, DynamicalSystems, StaticArrays, JLD2, Statistics

function map_LSE(sys, params,
    range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
    time_setting,
    matrix_u0, matrix_first_points, matrix_LSE; printing = false)
    #---------------------------------------
    len_p1 = length(range_p1)
    len_p2 = length(range_p2)
    dim = length(matrix_u0[1, 1, :])
    #---------------------------------------
    global matrix_EEs = zeros(len_p1, len_p2);
    #---------------------------------------
    namefile_EEs, namefile_print = get_file_name(len_p1, len_p2,
     name_p1, name_p2)
    #---------------------------------------
    if printing == true
        inheritance_from_matrix(sys, params, range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
        time_setting, matrix_u0, matrix_first_points, matrix_LSE,
        namefile_EEs, namefile_print);
    end
end