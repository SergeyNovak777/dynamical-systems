function get_solve(prob, integrator_setting; last_point = false)

    if last_point == false
        sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
            abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
            maxiters = integrator_setting.maxiters);
    else
        sol = sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
            abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
            maxiters = integrator_setting.maxiters,
            ave_everystep = false, save_start = false);
    end

    return sol;
end

