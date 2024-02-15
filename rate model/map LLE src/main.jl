include("/home/sergey/work/repo/dynamical-systems/rate model/map LLE/header.jl")

function map_LLE(sys, params, u0,
    range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
    time_setting, integrator_setting; printing = false)

    #---------------------------------------
    len_p1 = length(range_p1)
    len_p2 = length(range_p2)
    dim = length(u0)

    #---------------------------------------
    global λs = zeros(len_p1, len_p2, 1)
    global u0s = zeros(len_p1, len_p2, dim*2)

    #---------------------------------------
    namefile_LLE, namefile_u0s = get_file_name(len_p1, len_p2,
     name_p1, name_p2)
    
    #---------------------------------------
    index_fixed_parameter = index_p1
    value_fixed_parameter = range_p1[1]
    index_control_parameter = index_p2

    if printing == false
        pre_broaching_without_print(sys, params, u0, time_setting, integrator_setting,
          index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
          namefile_LLE, namefile_u0s)
    else
        pre_broaching_with_print(sys, params, u0, time_setting, integrator_setting,
        index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
        name_p1, name_p2,
        namefile_LLE, namefile_u0s)
    end

    #---------------------------------------
    if printing == false
        calculate_map_without_print(sys, params, u0, time_setting, integrator_setting,
            index_parameter_1, index_parameter_2, range_p1, range_p2,
            namefile_LLE, namefile_u0s)
    else
        calculate_map_with_print(sys, params, u0, time_setting, integrator_setting,
            index_parameter_1, index_parameter_2, range_p1, range_p2,
            name_p1, name_p2,
            namefile_LLE, namefile_u0s)
    end

end

function calculate_map_without_print(sys, params, u0, time_setting, integrator_setting,
    index_parameter_1, index_parameter_2, range_p1, range_p2,
    namefile_LLE, namefile_u0s, dim = length(u0))

    for (index_cycle_parameter_2, value_parameter_2) in enumerate(range_p2)
        for (index_cycle_parameter_1, value_parameter_1) in enumerate(range_p1)
            
            if index_parameter_1 == 1
                continue
            end

            u0_local_map = u0s[index_parameter_1 - 1, index_parameter_2, dim+1:end]

            ds = init_ODE_prob(sys, params, u0_local_map, time_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)
            
            point_from_attractor = goto_attractor(ds, integrator_setting)

            ds = init_Coupled_ODE(sys, params, point_from_attractor, integrator_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)

            LLE = calculate_LLE(ds, time_setting.time_calculate_LLE)

            save_in_matrix(index_cycle_parameter_1, index_cycle_parameter_2, dim,
                LLE, u0_local_pre_broach, point_from_attractor)
        
                u0_local_map = point_from_attractor
    
            if mod(local_index_p2, 10) == 0
                save_in_files(namefile_LLE, namefile_u0s, dim)
            end

        end
    end

end

function calculate_map_with_print(sys, params, u0, time_setting, integrator_setting,
    index_parameter_1, index_parameter_2, range_p1, range_p2,
    name_p1, name_p2,
    namefile_LLE, namefile_u0s, dim = length(u0))

    for (index_cycle_parameter_2, value_parameter_2) in enumerate(range_p2)
        for (index_cycle_parameter_1, value_parameter_1) in enumerate(range_p1)
            
            if index_cycle_parameter_1 == 1
                continue
            end

            u0_local_map = u0s[index_cycle_parameter_1 - 1, index_cycle_parameter_2, dim+1:end]

            ds = init_ODE_prob(sys, params, u0_local_map, time_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)
            
            point_from_attractor = goto_attractor(ds, integrator_setting)

            ds = init_Coupled_ODE(sys, params, point_from_attractor, integrator_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)

            LLE = calculate_LLE(ds, time_setting.time_calculate_LLE)

            save_in_matrix(index_cycle_parameter_1, index_cycle_parameter_2, dim,
                LLE, u0_local_map, point_from_attractor)
                
            println("index_cycle_$(name_p2): $(index_cycle_parameter_2); value_$(name_p2): $(value_parameter_2)"); flush(stdout)
            println("index_cycle_$(name_p1): $(index_cycle_parameter_1); value_$(name_p1): $(value_parameter_1)");flush(stdout)

            println("init_point: $(u0_local_map)");flush(stdout)
            println("last_point: $(point_from_attractor)");flush(stdout)

            println("LLE: $(LLE)");flush(stdout)

            println("-------------------------------------------");flush(stdout)
            println("");flush(stdout)

            u0_local_map = point_from_attractor
        end
        save_in_files(namefile_LLE, namefile_u0s, dim)
    end

end

function pre_broaching_without_print(sys, params, u0, time_setting, integrator_setting,
    index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
    namefile_LLE, namefile_u0s, dim = length(u0))

    for (local_index_p2, value_p2) in enumerate(range_p2)

        if local_index_p2 == 1
            global u0_local_pre_broach = u0
        end
    
        ds = init_ODE_prob(sys, params, u0_local_pre_broach, time_setting,
        index_control_parameter, value_p2,
        index_fixed_parameter, value_fixed_parameter)
    
        point_from_attractor = goto_attractor(ds, integrator_setting)

        ds = init_Coupled_ODE(sys, params, point_from_attractor, integrator_setting,
        index_control_parameter, value_p2,
        index_fixed_parameter, value_fixed_parameter)
    
        LLE = calculate_LLE(ds, time_setting.time_calculate_LLE)
    
        save_in_matrix(1, local_index_p2, dim,
        LLE, u0_local_pre_broach, point_from_attractor)
        
        u0_local_pre_broach = point_from_attractor
    
        if mod(local_index_p2, 10) == 0
            save_in_files(namefile_LLE, namefile_u0s, dim)
        end
    
    end

end

function pre_broaching_with_print(sys, params, u0, time_setting, integrator_setting,
    index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
    name_p1, name_p2,
    namefile_LLE, namefile_u0s, dim = length(u0))

    for (local_index_p2, value_p2) in enumerate(range_p2)

        if local_index_p2 == 1
            global u0_local_pre_broach = u0
        end
    
        ds = init_ODE_prob(sys, params, u0_local_pre_broach, time_setting,
        index_control_parameter, value_p2,
        index_fixed_parameter, value_fixed_parameter)
    
        point_from_attractor = goto_attractor(ds, integrator_setting)

        ds = init_Coupled_ODE(sys, params, point_from_attractor, integrator_setting,
        index_control_parameter, value_p2,
        index_fixed_parameter, value_fixed_parameter)
    
        LLE = calculate_LLE(ds, time_setting.time_calculate_LLE)
    
        save_in_matrix(1,local_index_p2, dim,
        LLE, u0_local_pre_broach, point_from_attractor)
        
        println(params);flush(stdout)

        println("$(name_p1): $(params[index_fixed_parameter])");flush(stdout)
        println("index_cycle_$(name_p2): $(local_index_p2), value_$(name_p2): $(value_p2)");flush(stdout)

        println("init_point: $(u0_local_pre_broach)");flush(stdout)
        println("last_point: $(point_from_attractor)");flush(stdout)

        println("LLE: $(LLE)");flush(stdout)

        println("-------------------------------------------");flush(stdout)
        println("");flush(stdout)
        u0_local_pre_broach = point_from_attractor
    
        if mod(local_index_p2, 10) == 0
            save_in_files(namefile_LLE, namefile_u0s, dim)
        end
    
    end

end