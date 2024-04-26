function change_control_parameters(parameters, cycle_index_parameter_p1, value_parameter_p1, cycle_index_parameter_p2, value_parameter_p2)
    parametersp[cycle_index_parameter_p1] = value_parameter_p1
    parameters[cycle_index_parameter_p2] = value_parameter_p2
    return parameters
end

function get_point_from_attractor(sys, parameters, u0, t_integrate, integrator_setting)

    prob = ODEProblem(sys, SVector{length(u0_lc)}(u0), (0.0, t_integrate), parameters)

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

function get_LSE(sys, parameters, u0, t_LSE, integrator_setting)
    ds = CoupledODEs(sys, u0, parameters, diffeq = integrator_setting);
    LSE = lyapunovspectrum(ds, t_LSE)
    return LSE
end

function get_file_name(len_p1, len_p2, name_p1, name_p2)
    map_dim = "_$(len_p1)x$(len_p2)_"
    name = "$(name_p1)_$(name_p2)"
    format = ".jld2"
    namefile_LSE = "LSE" * map_dim * name * format
    namefile_u0s = "u0s" * map_dim * name * format
    return namefile_LSE, namefile_u0s
end

function save_in_matrix(λs, u0s, index_cycle_p1, index_cycle_p2, dim,
    LSE, u0_local_pre_broach, point_from_attractor)
    λs[index_cycle_p1, index_cycle_p2, :] = LSE
    u0s[index_cycle_p1, index_cycle_p2, 1:dim] = u0_local_pre_broach
    u0s[index_cycle_p1, index_cycle_p2, dim+1:end] = point_from_attractor
end

function save_in_files(namefile_LSE, namefile_u0s, λs, init_points, last_points, dim)
    jldsave(namefile_LSE; λs)
    jldsave(namefile_u0s, init_points = u0s[:, :, 1:dim], last_points = u0s[:, :, dim+1:end])
end

function save_in_matrix(λs, u0s, index_cycle_p1, index_cycle_p2, dim,
    LSE, point_from_attractor)
    λs[index_cycle_p1, index_cycle_p2, :] = LSE
    u0s[index_cycle_p1, index_cycle_p2, dim+1:end] = point_from_attractor
end

function save_in_files(namefile_LSE, namefile_u0s, λs, last_points, dim)
    jldsave(namefile_LSE; λs)
    jldsave(namefile_u0s, init_points = u0s[:, :, 1:dim], last_points = u0s[:, :, dim+1:end])
end