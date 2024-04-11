function get_file_name(len_p1, len_p2, name_p1, name_p2)
    
    map_dim = "_$(len_p1)x$(len_p2)_"
    name = "$(name_p1)_$(name_p2)"
    format = ".jld2"
    namefile_LSE = "LSE" * map_dim * name * format
    namefile_u0s = "u0s" * map_dim * name * format
    return namefile_LSE, namefile_u0s
end

function get_percent(number, percent)
    return floor(Int64, (number / 100) * percent )
end

function get_point_from_attractor(prob, integrator_setting)

    if integrator_setting.adaptive == true
        last_point_from_attractor = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, maxiters = integrator_setting.maxiters,
        save_everystep = false, save_start = false)
    else
        last_point_from_attractor = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        dt = integrator_setting.dt,
        save_everystep = false, save_start = false, maxiters = integrator_setting.maxiters)
    end
    return last_point_from_attractor[end]
end

function init_ODE_prob(sys, parameters, u0, time_setting,
    index_parameter_p2, value_parameter_p2,index_parameter_p1, value_parameter_p1)

    parameters[index_parameter_p2] = value_parameter_p2
    parameters[index_parameter_p1] = value_parameter_p1
    tspan = (0.0, time_setting.time_attract)
    prob = ODEProblem(sys, SVector{length(u0)}(u0), tspan, parameters)
    return prob
end

function init_Coupled_ODE(sys, parameters, u0, integ_set,
     index_parameter_p2, value_parameter_p2,index_parameter_p1, value_parameter_p1)

    parameters[index_parameter_p2] = value_parameter_p2;
    parameters[index_parameter_p1] = value_parameter_p1;
    ds = CoupledODEs(sys, u0, parameters, diffeq = integ_set);
    return ds
end

function solver(prob, integrator_setting)
    if integrator_setting.adaptive == true
        solution = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, maxiters = integrator_setting.maxiters)
    else
        solution = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        dt = integrator_setting.dt,
        maxiters = integrator_setting.maxiters)
    end
    return solution
end

function distance_between_point(X1, X2)
    return sqrt( (X1[1] - X2[1])^2 + ( X1[2] - X2[2] )^2 + ( X1[3] - X2[3] )^2 )
end

function difference_between_points(solution, len_matrix)
    matrix_difference = zeros(len_matrix)
    for index in range(1, len_matrix - 1, step = 2)
        point1 = solution[index]
        point2 = solution[index+1]
        matrix_difference[index] = distance_between_point(point1, point2)
    end
    return matrix_difference
end

function calculate_LSE(ds, time_calculate_LSE)
    LSE = lyapunovspectrum(ds, time_calculate_LSE)
    return LSE
end

function save_in_matrix(index_cycle_p1, index_cycle_p2, dim,
    LSE, u0, last_point_from_attractor)

    λs[index_cycle_p1, index_cycle_p2, :] = LSE
    u0s[index_cycle_p1, index_cycle_p2, 1:dim] = u0
    u0s[index_cycle_p1, index_cycle_p2, dim+1:end] = last_point_from_attractor
end

function save_in_matrix(index_cycle_p1, index_cycle_p2, dim,
    LSE, last_point_from_attractor)

    λs[index_cycle_p1, index_cycle_p2, :] = LSE
    u0s[index_cycle_p1, index_cycle_p2, dim+1:end] = last_point_from_attractor
end

function save_in_files(namefile_LSE, namefile_u0s, dim)

    jldsave(namefile_LSE; λs)
    jldsave(namefile_u0s, init_points = u0s[:, :, 1:dim], last_points = u0s[:, :, dim+1:end])
end

function print_output_inh(name_p1, name_p2, index_cycle_parameter_1, index_cycle_parameter_2, value_parameter_1, value_parameter_2,
    u0, point_from_attractor, LSE)

    println("index_cycle_$(name_p1): $(index_cycle_parameter_1); value_$(name_p1): $(value_parameter_1)"); flush(stdout)
    println("index_cycle_$(name_p2): $(index_cycle_parameter_2); value_$(name_p2): $(value_parameter_2)"); flush(stdout)
    println("init_point: $(u0)"); flush(stdout)
    println("last_point: $(point_from_attractor)"); flush(stdout)
    println("LSE: $(LSE)"); flush(stdout)
    println("-------------------------------------------"); flush(stdout)
    println(""); flush(stdout)
end

function print_output_inh(name_p1, name_p2,
    index_parameter_in_cycle, value_parameter_in_cycle, value_fixed_parameter,
    u0, u0end, LSE
    )

    println("$name_p1: $value_fixed_parameter"); flush(stdout)
    println("index parameter $name_p2 in cycle: $(index_parameter_in_cycle); value $name_p2: $(value_parameter_in_cycle)"); flush(stdout)
    println("init_point: $u0"); flush(stdout)
    println("last_point: $u0end"); flush(stdout)
    println("LSE: $LSE"); flush(stdout)
    println("-------------------------------------------"); flush(stdout)
    println(""); flush(stdout)
end
#--------------------------------------------------------------------------------------------------
#Inheritance to side

function pre_broaching_inh_side(sys, params, u0, time_setting, integrator_setting,
    index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
    name_p1, name_p2,
    namefile_LLE, namefile_u0s, 
    flag_output, dim = length(u0))

    for (index_parameter_p2_in_cycle, value_parameter_p2_in_cycle) in enumerate(range_p2)

        if index_parameter_p2_in_cycle == 1
            global u0_local_pre_broach = u0
        end
    
        ds = init_ODE_prob(sys, params, u0_local_pre_broach, time_setting,
        index_control_parameter, value_parameter_p2_in_cycle,
        index_fixed_parameter, value_fixed_parameter)
    
        point_from_attractor = get_point_from_attractor(ds, integrator_setting)

        ds = init_Coupled_ODE(sys, params, point_from_attractor, integrator_setting,
        index_control_parameter, value_parameter_p2_in_cycle,
        index_fixed_parameter, value_fixed_parameter)
    
        LSE = calculate_LSE(ds, time_setting.time_calculate_LSE)
    
        save_in_matrix(1,index_parameter_p2_in_cycle, dim,
        LSE, u0_local_pre_broach, point_from_attractor)

        if flag_output == true
            print_output_inh(name_p1, name_p2,
            index_parameter_p2_in_cycle, value_parameter_p2_in_cycle, params[index_fixed_parameter],
            u0_local_pre_broach, point_from_attractor, LSE)
        end
        
        u0_local_pre_broach = point_from_attractor
    
        if mod(index_parameter_p2_in_cycle, 10) == 0
            save_in_files(namefile_LLE, namefile_u0s, dim)
        end
    end
end

function calculate_map_inh_side(sys, params, u0, time_setting, integrator_setting,
    index_parameter_1, index_parameter_2, range_p1, range_p2,
    name_p1, name_p2,
    namefile_LLE, namefile_u0s,
    flag_output, dim = length(u0))

    for (index_cycle_parameter_2, value_parameter_2) in enumerate(range_p2)
        for (index_cycle_parameter_1, value_parameter_1) in enumerate(range_p1)
            
            if index_cycle_parameter_1 == 1
                continue
            end

            u0_local_map = u0s[index_cycle_parameter_1 - 1, index_cycle_parameter_2, dim+1:end]

            ds = init_ODE_prob(sys, params, u0_local_map, time_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)
            
            point_from_attractor = get_point_from_attractor(ds, integrator_setting)

            ds = init_Coupled_ODE(sys, params, point_from_attractor, integrator_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)

            LSE = calculate_LSE(ds, time_setting.time_calculate_LSE)

            save_in_matrix(index_cycle_parameter_1, index_cycle_parameter_2, dim,
                LSE, u0_local_map, point_from_attractor)

            if flag_output == true
                print_output_inh(name_p1, name_p2, index_cycle_parameter_1, index_cycle_parameter_2, value_parameter_1, value_parameter_2,
                u0, point_from_attractor, LSE)
            end

            u0_local_map = point_from_attractor
        end
        save_in_files(namefile_LLE, namefile_u0s, dim)
    end
end

#--------------------------------------------------------------------------------------------------
#Inheritance to side with detect FP

function pre_broaching_inh_side(sys, params, u0, time_setting, integrator_setting,
    index_fixed_parameter, value_fixed_parameter, index_control_parameter, range_p2,
    name_p1, name_p2,
    namefile_LLE, namefile_u0s, 
    flag_output, ϵ; dim = length(u0))

    for (index_parameter_p2_in_cycle, value_parameter_p2_in_cycle) in enumerate(range_p2)

        if index_parameter_p2_in_cycle == 1
            global u0_local_pre_broach = u0
        end
    
        ds = init_ODE_prob(sys, params, u0_local_pre_broach, time_setting,
        index_control_parameter, value_parameter_p2_in_cycle,
        index_fixed_parameter, value_fixed_parameter)
    
        solution = solver(ds, integrator_setting)

        len_sol = length(solution)
        percent = get_percent(len_sol, 100-50)
        len_matrix = length(solution[:, percent:end])
        percent_sol = solution[:, percent:end]

        if mod(len_matrix, 2) != 0
            len_matrix = len_matrix - 1
        end
    
        matrix_difference = difference_between_points(percent_sol, len_matrix)

        if all(matrix_difference .<= ϵ) == true
            LSE = ones(dim) * -1
        else
            ds = init_Coupled_ODE(sys, params, solution[end], integrator_setting,
            index_control_parameter, value_parameter_p2_in_cycle,
            index_fixed_parameter, value_fixed_parameter)
        
            LSE = calculate_LSE(ds, time_setting.time_calculate_LSE)
        end
       
    
        save_in_matrix(1,index_parameter_p2_in_cycle, dim,
        LSE, u0_local_pre_broach, solution[end])

        if flag_output == true
            print_output_inh(name_p1, name_p2,
            index_parameter_p2_in_cycle, value_parameter_p2_in_cycle, params[index_fixed_parameter],
            u0_local_pre_broach, solution[end], LSE)
        end
        
        u0_local_pre_broach = solution[end]
    
        if mod(index_parameter_p2_in_cycle, 10) == 0
            save_in_files(namefile_LLE, namefile_u0s, dim)
        end
    end
end

function calculate_map_inh_side(sys, params, u0, time_setting, integrator_setting,
    index_parameter_1, index_parameter_2, range_p1, range_p2,
    name_p1, name_p2,
    namefile_LLE, namefile_u0s,
    flag_output, ϵ; dim = length(u0))

    for (index_cycle_parameter_2, value_parameter_2) in enumerate(range_p2)
        for (index_cycle_parameter_1, value_parameter_1) in enumerate(range_p1)
            
            if index_cycle_parameter_1 == 1
                continue
            end

            u0_local_map = u0s[index_cycle_parameter_1 - 1, index_cycle_parameter_2, dim+1:end]

            ds = init_ODE_prob(sys, params, u0_local_map, time_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)
            
            solution = solver(ds, integrator_setting)

            len_sol = length(solution)
            percent = get_percent(len_sol, 100-50)
            len_matrix = length(solution[:, percent:end])
            percent_sol = solution[:, percent:end]

            if mod(len_matrix, 2) != 0
                len_matrix = len_matrix - 1
            end
    
            matrix_difference = difference_between_points(percent_sol, len_matrix)

            if all(matrix_difference .<= ϵ) == true
                LSE = ones(dim) * -1
            else
                ds = init_Coupled_ODE(sys, params, solution[end], integrator_setting,
                index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)

                LSE = calculate_LSE(ds, time_setting.time_calculate_LSE)
            end

            save_in_matrix(index_cycle_parameter_1, index_cycle_parameter_2, dim,
                LSE, u0_local_map, solution[end])

            if flag_output == true
                print_output_inh(name_p1, name_p2, index_cycle_parameter_1, index_cycle_parameter_2, value_parameter_1, value_parameter_2,
                u0_local_map, solution[end], LSE)
            end

            u0_local_map = solution[end]
        end
        save_in_files(namefile_LLE, namefile_u0s, dim)
    end
end
#--------------------------------------------------------------------------------------------------
#Without Inheritance

function calculate_map_LSE_without_inheritance(sys, params, u0, time_setting, integrator_setting,
    index_parameter_1, index_parameter_2, range_p1, range_p2,
    namefile_LSE, namefile_u0s,
    flag_output,  dim = length(u0))

    for (index_cycle_parameter_2, value_parameter_2) in enumerate(range_p2)
        for (index_cycle_parameter_1, value_parameter_1) in enumerate(range_p1)
            
            ds = init_ODE_prob(sys, params, u0, time_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)
            
            point_from_attractor = get_point_from_attractor(ds, integrator_setting)

            ds = init_Coupled_ODE(sys, params, point_from_attractor, integrator_setting,
            index_parameter_2, value_parameter_2, index_parameter_1, value_parameter_1)

            LSE = calculate_LSE(ds, time_setting.time_calculate_LSE)

            save_in_matrix(index_cycle_parameter_1, index_cycle_parameter_2, dim,
                LSE, point_from_attractor)
            
            if flag_output == true
                print_output_w_inh(name_p1, name_p2, index_cycle_parameter_1, index_cycle_parameter_2, value_parameter_1, value_parameter_2,
                point_from_attractor, LSE)
            end

            if mod(local_index_p2, 10) == 0
                save_in_files(namefile_LSE, namefile_u0s, dim)
            end
        end
    end
end

function print_output_w_inh(name_p1, name_p2, index_cycle_parameter_1, index_cycle_parameter_2, value_parameter_1, value_parameter_2,
    point_from_attractor, LSE)

    println("index_cycle_$(name_p1): $(index_cycle_parameter_1); value_$(name_p1): $(value_parameter_1)"); flush(stdout)
    println("index_cycle_$(name_p2): $(index_cycle_parameter_2); value_$(name_p2): $(value_parameter_2)"); flush(stdout)
    println("last_point: $(point_from_attractor)"); flush(stdout)
    println("LSE: $(LSE)"); flush(stdout)
    println("-------------------------------------------"); flush(stdout)
    println(""); flush(stdout)
end