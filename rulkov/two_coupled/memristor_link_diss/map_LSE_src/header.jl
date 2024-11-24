#= for discrete system =#
function get_solve_sys_d(type_get, sys)

    if type_get == "last point"
        sol = get_last_point_d(sys);
    elseif type_get == "full solve"
        sol = get_full_solve(sys);
    end

    return sol;
end

function get_last_point_d(sys)
    return solve(sys, save_everystep = false, save_start = false);
end

function get_full_solve_d(sys)
    return solve(sys);
end

#= ------------------------------------------------------------- =#
#= total function =#
function get_LSE(sys, u0, t_LSE)
    LSE = lyapunovspectrum(sys, t_LSE, u0 = u0);
end
#= ------------------------------------------------------------- =#
#= for cont system =#
#= ------------------------------------------------------------- =#


#= inheritance =#
function left_right(prob, ds, u0, params,
    time,
    range_p1, range_p2, name_p1, name_p2,
    index_p1, index_p2;
    flag_print = false)

    

    for index_p1 in range_p1
        for index_p2 in range_p2
            
        end
    end
end
