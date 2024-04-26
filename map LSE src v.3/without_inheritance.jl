include(joinpath(@__DIR__, "general_function.jl"))

function map_without_inheritance(sys, parameters, u0,
    index_control_parameter_p1, index_control_parameter_p2, range_p1, range_p2,
    time_setting, integrator_setting; printing = false)

    #---------------------------------------
    len_p1 = length(range_p1)
    len_p2 = length(range_p2)
    dim = length(u0)
    #---------------------------------------
    global λs = zeros(len_p1, len_p2, dim)
    global u0s = zeros(len_p1, len_p2, dim)
    #---------------------------------------
    namefile_LSE, namefile_u0s = get_file_name(len_p1, len_p2,
    name_p1, name_p2)
    #---------------------------------------

    for (cycle_index_parameter_p1, value_parameter_p1) in enumerate(range_p1)
        for (cycle_index_parameter_p2, value_parameter_p2) in enumerate(range_p2)

            parameters = change_control_parameters(parameters, index_control_parameter_p1, value_parameter_p1,
                                                    index_control_parameter_p2, value_parameter_p2)

            point_from_attractor = get_point_from_attractor(sys, parameters, u0, time_setting.t_integrate, integrator_setting)
            LSE = get_LSE(sys, parameters, point_from_attractor, time_setting.t_LSE, integrator_setting)

            save_in_matrix(λs, u0s, cycle_index_parameter_p1, cycle_index_parameter_p2, dim,
            LSE, point_from_attractor)
    
            if mod(cycle_index_parameter_p2, 10) == 0
                save_in_files(namefile_LSE, namefile_u0s, λs, point_from_attractor, dim)
            end
        end
    end
end


