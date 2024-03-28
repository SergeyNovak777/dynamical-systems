function get_file_name(len_p1, len_p2, name_p1, name_p2)
    
    map_dim = "_$(len_p1)x$(len_p2)_"
    name = "$(name_p1)_$(name_p2)"
    format = ".jld2"
    namefile_LSE = "LSE" * map_dim * name * format
    namefile_u0s = "u0s" * map_dim * name * format
    return namefile_LSE, namefile_u0s

end

function init_ODE_prob(sys, params, u0_lc, time_setting,
    index_parameter_p2, value_parameter_p2,index_parameter_p1, value_parameter_p1)

    params[index_parameter_p2] = value_parameter_p2
    params[index_parameter_p1] = value_parameter_p1
    tspan = (0.0, time_setting.time_attract)
    prob = ODEProblem(sys, SVector{length(u0_lc)}(u0_lc), tspan, params)

    return prob
end

function init_Coupled_ODE(sys, params, u0_lc, integ_set,
     index_parameter_p2, value_parameter_p2,index_parameter_p1, value_parameter_p1)

    params[index_parameter_p2] = value_parameter_p2;
    params[index_parameter_p1] = value_parameter_p1;
    
    ds = CoupledODEs(sys, u0_lc, params, diffeq = integ_set);

    return ds
end

function goto_attractor(prob, integrator_setting)

    if integrator_setting.adaptive == true
        point_from_attractor = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, maxiters = integrator_setting.maxiters,
        save_everystep = false, save_start = false)
    else
        point_from_attractor = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        dt = integrator_setting.dt,
        save_everystep = false, save_start = false, maxiters = integrator_setting.maxiters)
    end

    return point_from_attractor[end]
end

function calculate_LSE(ds, time_calculate_LSE)
    LSE = lyapunovspectrum(ds, time_calculate_LSE)
    return LSE
end

function save_in_matrix(index_cycle_p1, index_cycle_p2, dim,
    LSE, u0_local_pre_broach, point_from_attractor)

    λs[index_cycle_p1, index_cycle_p2, :] = LSE
    u0s[index_cycle_p1, index_cycle_p2, 1:dim] = u0_local_pre_broach
    u0s[index_cycle_p1, index_cycle_p2, dim+1:end] = point_from_attractor

end

function save_in_files(namefile_LSE, namefile_u0s, dim)
    jldsave(namefile_LSE; λs)
    jldsave(namefile_u0s, init_points = u0s[:, :, 1:dim], last_points = u0s[:, :, dim+1:end])
end