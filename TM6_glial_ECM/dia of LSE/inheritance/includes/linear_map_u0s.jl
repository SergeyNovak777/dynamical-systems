#=
name_p1 = "I0";
name_p2 = "αE";
index_p1 = 11;
index_p2 = 6;
return [name_p1, name_p2, index_p1, index_p2];

integ_set = (alg, adaptive, abstol, reltol)
=#

function solver(prob, integ_set)

    alg = integ_set.alg;
    adapt = integ_set.adaptive;
    abtol = integ_set.abstol;
    retol = integ_set.reltol;

    lastpoint = solve(prob,
    alg, adaptive = adapt, abstol = abtol, reltol = retol,
    save_everystep = false, save_start = false);
    
    return lastpoint;
end

#=
FUNCTION linear_2pmap
ARGUMENTS

*   sys -- system 
*   params -- vector of parameters
*   u0_initial -- initial condition
*   tspan -- Tuple which containt initial and end time
*   integral_setting -- Tuple, which contain setting for integrator

*   control_params -- vector, which contain names of control parameters and their index;
    example:
    [name_p1, name_p2, index_p1, index_p2] 

*   p1_range -- range first control parameters
*   p2_range -- range second control parameters
=#


function linear_2pmap(sys, params, u0_initial, tspan, integral_setting,
    control_params, p1_range, p2_range)

    D = length(u0_initial);
    u0_loop = zeros(D);
    #--------------------------------------------

    name_p1, name_p2, index_p1, index_p2 = control_params;
    #--------------------------------------------
    len_range = length(p1_range);
    matrix_u0_start = zeros(len_range, len_range, D);
    matrix_u0_end = zeros(len_range, len_range, D);
    #---------------------------------------------

    filename_matrix_u0_start = "linear map u0start $(name_p1) $(name_p2) $(len_range)x$(len_range)"; 
    filename_matrix_u0_end = "linear map u0end $(name_p1) $(name_p2) $(len_range)x$(len_range)";
    #---------------------------------------------

    exitcode = 0; # for detect MaxIters and etc
    #---------------------------------------------


    for p2_index_range in eachindex(p2_range)

        if p2_index_range == 1
            global u0_start = u0_initial;
        end;

        loc_params = copy(params);
        loc_params[index_p2] = p2_range[p2_index_range];

        prob = ODEProblem(sys, u0_start, tspan, loc_params);
        u0_end = solver(prob, integral_setting);

        if u0_end.retcode == ReturnCode.MaxIters
            exitcode = 1;
        end
        u0_end = u0_end[end]
        matrix_u0_start[1, p2_index_range, :] = u0_start;
        matrix_u0_end[1, p2_index_range, :] = u0_end; 

        u0_start = u0_end;

        if exitcode == 1
            exit();
        end;
    end

    return matrix_u0_start, matrix_u0_end;
end