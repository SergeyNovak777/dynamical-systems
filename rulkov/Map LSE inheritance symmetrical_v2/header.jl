function get_file_name(len_p1, len_p2, name_p1, name_p2)
    
    map_dim = "_$(len_p1)x$(len_p2)_"
    name = "$(name_p1)_$(name_p2)"
    format = ".jld2"
    namefile_LSE = "LSE" * map_dim * name * format;
    namefile_u0s = "u0s" * map_dim * name * format;
    namefile_last_points = "last_points" * map_dim * name * format
    namefile_print = "printing" * map_dim * name * ".txt"
    return namefile_LSE, namefile_u0s, namefile_last_points,namefile_print
end

function get_last_point_from_attractor(sys, params, u0, tspan)
    prob = DiscreteProblem(sys, u0, tspan, params, save_everystep = false, save_start = false);
    sol = solve(prob);
    return sol[:, end];
end

function get_LSE(sys, params, u0, t_LSE)
    ds = DeterministicIteratedMap(sys, u0, params);
    return lyapunovspectrum(ds, t_LSE);
end

function save_in_matrix(index_cycle_p1, index_cycle_p2, u0, last_point_from_attractor, LSE)
    λs[index_cycle_p1, index_cycle_p2, :] = LSE;
    u0s[index_cycle_p1, index_cycle_p2, :] = u0;
    last_points[index_cycle_p1, index_cycle_p2, :] = last_point_from_attractor;
end

function save_in_files(namefile_LSE, namefile_u0s, namefile_last_points)
    jldsave(namefile_LSE; λs);
    jldsave(namefile_u0s; u0s);
    jldsave(namefile_last_points; last_points);
end

function print_in_file(index_cycle_p1, index_cycle_p2, value_p1, value_p2, name_p1, name_p2,
    u0, last_point_from_attractor, LSE,
    namefile_print)

    open(namefile_print, "a") do io
        println(io, "index cycle $(name_p1): $index_cycle_p1); value p1 : $(value_p1)");
        println(io, "index cycle $(name_p2): $index_cycle_p2); value p2 : $(value_p2)");
        println(io, "u0: $(u0);");
        println(io, "last point: $(last_point_from_attractor)");
        println(io, "LSE: $(LSE)");
        println(io, "--------------------------------------------------");
        println(io, "");
    end
end

function inheritance_diagonal_print(sys, params, u0,
    range_p1, range_p2, index_p1, index_p2, name_p1, name_p2, length_range_p1, length_range_p2,
    time_setting,
    namefile_LSE, namefile_u0s, namefile_last_points, namefile_print)

        u0_local = u0;

        for index_cycle_p1 in 1:length_range_p1
            for index_cycle_p2 in 1:length_range_p2

                if index_cycle_p1 == index_cycle_p2
                    
                    u0_before = u0_local;
                    u0_local = get_x_y_z(u0_local);
                    u0_local = three_coupled_rulkov_first_iteration(u0_local, params);

                    params[index_p1] = range_p1[index_cycle_p1];
                    params[index_p2] = range_p2[index_cycle_p2];

                    last_point_from_attractor = get_last_point_from_attractor(sys, params, u0_local, time_setting.tspan);
                    LSE = get_LSE(sys, params, last_point_from_attractor, time_setting.t_LSE);

                    save_in_matrix(index_cycle_p1, index_cycle_p2, u0_local, last_point_from_attractor, LSE);

                    if mod(index_cycle_p1, 10) == 0
                        save_in_files(namefile_LSE, namefile_u0s, namefile_last_points)
                    end

                    print_in_file(index_cycle_p1, index_cycle_p2, params[index_p1], params[index_p2], name_p1, name_p2,
                    u0_before, last_point_from_attractor, LSE, namefile_print)

                    u0_local = last_point_from_attractor;

                    if index_cycle_p1 != length_range_p1
                        
                        inheritance_side1_print(sys, params, range_p1, range_p2,
                        index_p1, index_p2, index_cycle_p1 + 1, index_cycle_p2, name_p1, name_p2,
                        time_setting,
                        namefile_LSE, namefile_u0s, namefile_last_points, namefile_print);

                        inheritance_side2_print(sys, params, range_p1, range_p2,
                        index_p1, index_p2, index_cycle_p1, index_cycle_p2 + 1, name_p1, name_p2,
                        time_setting,
                        namefile_LSE, namefile_u0s, namefile_last_points, namefile_print)
                    end
                end
            end
        end         
end

function inheritance_side1_print(sys, params, range_p1, range_p2,
    index_p1, index_p2, index_cycle_p1, index_cycle_p2, name_p1, name_p2,
    time_setting,
    namefile_LSE, namefile_u0s, namefile_last_points, namefile_print)

    params[index_p2] = range_p2[index_cycle_p2];

    for index_cycle_second_p1 in index_cycle_p1:length(range_p1)

        params[index_p1] = range_p1[index_cycle_second_p1];

        if index_cycle_second_p1 == index_cycle_p1
            global u0_local = last_points[index_cycle_second_p1 - 1, index_cycle_p2, :];
        end

        u0_before = u0_local;

        u0_local = get_x_y_z(u0_local);
        u0_local = three_coupled_rulkov_first_iteration(u0_local, params);

        last_point_from_attractor = get_last_point_from_attractor(sys, params, u0_local, time_setting.tspan);
        LSE = get_LSE(sys, params, last_point_from_attractor, time_setting.t_LSE);

        save_in_matrix(index_cycle_second_p1, index_cycle_p2, u0_local, last_point_from_attractor, LSE);

        if mod(index_cycle_second_p1, 10) == 0
            save_in_files(namefile_LSE, namefile_u0s, namefile_last_points)
        end

        print_in_file(index_cycle_second_p1, index_cycle_p2, params[index_p1], params[index_p2], name_p1, name_p2,
                    u0_before, last_point_from_attractor, LSE, namefile_print)

        u0_local = last_point_from_attractor;
    end
end

function inheritance_side2_print(sys, params, range_p1, range_p2,
    index_p1, index_p2, index_cycle_p1, index_cycle_p2, name_p1, name_p2,
    time_setting,
    namefile_LSE, namefile_u0s, namefile_last_points, namefile_print)

    params[index_p1] = range_p1[index_cycle_p1];

    for index_cycle_second_p2 in index_cycle_p2:length(range_p2)
        
        params[index_p2] = range_p2[index_cycle_second_p2];

        if index_cycle_second_p2 == index_cycle_p2
            global u0_local = last_points[index_cycle_p1, index_cycle_second_p2 - 1, :];
        end

        u0_before = u0_local;

        u0_local = get_x_y_z(u0_local);
        u0_local = three_coupled_rulkov_first_iteration(u0_local, params);

        last_point_from_attractor = get_last_point_from_attractor(sys, params, u0_local, time_setting.tspan);
        LSE = get_LSE(sys, params, last_point_from_attractor, time_setting.t_LSE);

        save_in_matrix(index_cycle_p1, index_cycle_second_p2, u0_local, last_point_from_attractor, LSE);

        if mod(index_cycle_second_p2, 10) == 0
            save_in_files(namefile_LSE, namefile_u0s, namefile_last_points)
        end

        print_in_file(index_cycle_p1, index_cycle_second_p2, params[index_p1], params[index_p2], name_p1, name_p2,
                    u0_before, last_point_from_attractor, LSE, namefile_print)

        u0_local = last_point_from_attractor;
    end
end