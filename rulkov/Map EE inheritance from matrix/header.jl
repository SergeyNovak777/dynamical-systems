Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function get_file_name(len_p1, len_p2, name_p1, name_p2)
    
    map_dim = "_$(len_p1)x$(len_p2)_"
    name = "$(name_p1)_$(name_p2)"
    format = ".jld2"
    namefile_EEs = "EEs" * map_dim * name * format
    namefile_print = "printing" * map_dim * name * ".txt"
    return namefile_EEs,namefile_print
end
    
function get_sol(sys, params, u0, tspan)
    prob = DiscreteProblem(sys, u0, tspan, params);
    sol = solve(prob);
    xsum = sol[1, :] + sol[6, :] + sol[11, :];
    data = [xsum, sol.t]
    sol = nothing;
    return data;
    
end

function get_count_EE(data)

    data_local_max = get_local_max(data)
    data_local_min = get_local_min(data)

    drop_artifacts(data_local_max, data_local_min)

    Hs_xsum = Hs(data_local_max[1] ,6);
    count_EE = length(data_local_max[1][ data_local_max[1] .>= Hs_xsum ]);

    return count_EE
end

function save_in_matrix(index_cycle_p1, index_cycle_p2, count_EE)
    matrix_EEs[index_cycle_p1, index_cycle_p2] = count_EE;
end

function save_in_files(namefile_EEs)
    jldsave(namefile_EEs; matrix_EEs);
end

function print_in_file(index_cycle_p1, index_cycle_p2, value_p1, value_p2, name_p1, name_p2,
    u0, first_point, count_EE,
    namefile_print)

    open(namefile_print, "a") do io
        println(io, "index cycle $(name_p1): $index_cycle_p1); value p1 : $(value_p1)");
        println(io, "index cycle $(name_p2): $index_cycle_p2); value p2 : $(value_p2)");
        println(io, "first point $first_point");
        println(io, "u0: $(u0);");
        println(io, "count EE: $(count_EE)");
        println(io, "--------------------------------------------------");
        println(io, "");
    end
end

function inheritance_from_matrix(sys, params, range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
    time_setting, matrix_u0, matrix_first_points, matrix_LSE,
    namefile_EEs, namefile_print)

        length_range_p1 = length(range_p1);
        length_range_p2 = length(range_p2);

        for index_cycle_p1 in 1:length_range_p1
            for index_cycle_p2 in 1:length_range_p2

                    
                    #= u0_before = u0_local;
                    u0_local = get_x_y_z(u0_local);
                    u0_local = three_coupled_rulkov_first_iteration(u0_local, params); =#

                    params[index_p1] = range_p1[index_cycle_p1];
                    params[index_p2] = range_p2[index_cycle_p2];
                    u0 = matrix_u0[index_cycle_p1, index_cycle_p2, :];
                    first_point = matrix_first_points[index_cycle_p1, index_cycle_p2, :];

                    if matrix_LSE[index_cycle_p1, index_cycle_p2, 1] >= 0 
                        
                        data = get_sol(sys, params, u0, time_setting.tspan);
                        count_EE = get_count_EE(data);
                        save_in_matrix(index_cycle_p1, index_cycle_p2, count_EE);

                        print_in_file(index_cycle_p1, index_cycle_p2, params[index_p1], params[index_p2], name_p1, name_p2,
                        u0, first_point, count_EE,
                        namefile_print)
                    else
                        print_in_file(index_cycle_p1, index_cycle_p2, params[index_p1], params[index_p2], name_p1, name_p2,
                        u0, first_point, 0.0,
                        namefile_print)
                    end

                    if mod(index_cycle_p2, 10) == 0
                        save_in_files(namefile_EEs)
                    end
            end
        end         
end
