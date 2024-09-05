# typeof(prob).name.name == :DiscreteProblem
# typeof(prob).name.name == :ODEProblem

function get_last_point_from_attractor(prob, integrate_setting)
    if typeof(prob).name.name == :ODEProblem 
        if integrate_setting.adaptive == true
            sol = solve(prob, integrate_setting.alg, adaptive = integrate_setting.adaptive,
                    abstol =  integrate_setting.abstol, reltol = integrate_setting.reltol, maxiters = integrate_setting.maxiters,
                    save_everystep = false, save_start = false);
        else
            sol = solve(prob, integrate_setting.alg, adaptive = integrate_setting.adaptive,
                    dt = integrate_setting.dt, maxiters = integrate_setting.maxiters,
                    save_everystep = false, save_start = false);
        end
    elseif typeof(prob).name.name == :DiscreteProblem
        sol = solve(prob, save_everystep = false, save_start = false);
    end

    return sol[:, end];
end

function get_LSE(ds, t_LSE)
    LSE = lyapunovpsectrum(ds, t_LSE);
    LSE = sort(LSE, rev = true);
    return LSE;
end

function save_in_matrix_dia(LSE, u0s, index_p)
    matrix_LSE[index_p] = LSE;
    matrix_u0s[index_p] = u0s;
end

function save_in_files(filename_matrix_LSE, filename_matrix_u0s, matrix_LSE, matrix_u0s)
    jldsave(filename_matrix_LSE; matrix_LSE);
    jldsave(filename_matrix_u0s; matrix_u0s);
end

function diagram_LSE(prob, ds, range_parameter, index_parameter, name_parameter)

    for index in range(1, len_range, step = 1 )
        
        local_parameter = range_parameter[index];
        prob.p[index_parameter] = local_parameter;
        set_parameter!(system, index_parameter, local_parameter);

        
    end
end