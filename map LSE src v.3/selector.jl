include(joinpath(@__DIR__, "without_inheritance.jl"))

function selector(type_inheritance)
    if type_inheritance == "no inheritance"
        func = map_without_inheritance
    end

    return func
end

function calculate_map_LSE(sys, parameters, u0,
    index_control_parameter_p1, index_control_parameter_p2, range_p1, range_p2,
    time_setting, integrator_settingprinting)

    func = selector(type_inheritance)

    func(sys, parameters, u0,
    index_control_parameter_p1, index_control_parameter_p2, range_p1, range_p2,
    time_setting, integrator_setting; printing = false)
end