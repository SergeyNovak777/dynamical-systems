function solver(prob, integrator_setting)

    if integrator_setting.adaptive == true
        last_point = solve(prob, alg = integrator_setting.alg, adaptive = true,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, maxiters = integrator_setting.maxiters,
        save_everystep = false, save_start = false);
    else
        last_point = solve(prob, alg = integrator_setting.alg, adaptive = false,
        dt = integrator_setting.dt, maxiters = integrator_setting.maxiters,
        save_everystep = false, save_start = false);
    end

    return last_point[end];
end

function get_file_name(len_p1, len_p2, name_p1, name_p2)
    
    map_dim = "_$(len_p1)x$(len_p2)_"
    name = "$(name_p1)_$(name_p2)"
    format = ".jld2"
    namefile_LSE = "LSE" * map_dim * name * format
    namefile_u0s = "u0s" * map_dim * name * format
    
    return namefile_LSE, namefile_u0s

end

function save_in_matrix(cycle_index_p1, cycle_index_p2, u0, last_point, LSE)
    u0s[cycle_index_p1, cycle_index_p2, 1] = u0;
    u0s[cycle_index_p1, cycle_index_p2, 1] = last_point;
    λs[cycle_index_p1, cycle_index_p2, :] = LSE;
end

function save_in_files()
    jldsave(namefile_u0s; u0s);
    jldsave(namefile_LSE; λs);
end