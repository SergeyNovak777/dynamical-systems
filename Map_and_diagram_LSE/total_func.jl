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

