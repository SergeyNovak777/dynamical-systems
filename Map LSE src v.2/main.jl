#---------------------------------------
#INCLUDE
include(joinpath(@__DIR__, "header.jl"))

function map_LSE(sys, params, u0,
    range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
    time_setting, integrator_setting, type_inheritance; printing = false, ϵ = 1e-9)

    #---------------------------------------
    len_p1 = length(range_p1)
    len_p2 = length(range_p2)
    dim = length(u0)

    #---------------------------------------
    global λs = zeros(len_p1, len_p2, dim)
    global u0s = zeros(len_p1, len_p2, dim*2)

    #---------------------------------------
    namefile_LSE, namefile_u0s = get_file_name(len_p1, len_p2,
     name_p1, name_p2)
    
    #---------------------------------------
    index_fixed_parameter = index_p1
    value_fixed_parameter = range_p1[1]
    index_control_parameter = index_p2

    if type_inheritance == "move to side"

        pre_broaching_inh_side(sys, params, u0, time_setting, integrator_setting,
            index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
            name_p1, name_p2,
            namefile_LSE, namefile_u0s, printing)

        calculate_map_inh_side(sys, params, u0, time_setting, integrator_setting,
                index_p1, index_p2, range_p1, range_p2,
                name_p1, name_p2,
                namefile_LSE, namefile_u0s,
                printing)
                
    elseif type_inheritance == "no inheritance"
        calculate_map_LSE_without_inheritance(sys, params, u0, time_setting, integrator_setting,
        index_parameter_1, index_parameter_2, range_p1, range_p2,
        namefile_LSE, namefile_u0s, printing)

    elseif type_inheritance == "move to side detect fp"

        pre_broaching_inh_side(sys, params, u0, time_setting, integrator_setting,
            index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
            name_p1, name_p2,
            namefile_LSE, namefile_u0s, printing,
            ϵ)

        calculate_map_inh_side(sys, params, u0, time_setting, integrator_setting,
                index_p1, index_p2, range_p1, range_p2,
                name_p1, name_p2,
                namefile_LSE, namefile_u0s,
                printing,
                ϵ)
        
    end
end