#=
    time_calculate_LSE
    time_attract
    tstep

    algo
    adapt
    dt
=#
function map_LSE(sys, params, u0,
    range_p1, range_p2, index_p1, index_p2, name_p1, name_p2,
    time_setting, integrator_setting, printing = false)

    #---------------------------------------
    len_p1 = length(range_p1)
    len_p2 = length(range_p2)
    dim = length(u0)

    #---------------------------------------
    global λs = zeros(len_p1, len_p2, dim)
    global u0s = zeros(len_p1, len_p2, dim)

    #---------------------------------------
    namefile_LSE, namefile_u0s = get_file_name(len_p1, len_p2,
     name_p1, name_p2)
    
    #---------------------------------------
    index_fix = index_p1
    var_fix = range_p1[1]
    index_control = index_p2

    if printing == false
        pre_broaching_without_print(sys, params, u0, time_setting, integrator_setting,
          index_fix, var_fix, index_control, range_p2)
    else
        pre_broaching_with_print(sys, params, u0, time_setting, integrator_setting,
        index_fix, var_fix, index_control, range_p2)
    end

    #---------------------------------------



end

function pre_broaching_without_print(sys, params, u0, time_setting, integrator_setting,
    index_fix, var_fix, index_control, range_p2)

    for (p2_loc_index, p2_loc) in enumerate(range_p2)
    
        if p2_loc_index == 1
            global u0_lc = u0
        end
        
        ds = init_ds_(sys, params, index_control, p2_loc,
        index_fix, var_fix, u0_lc, integ_set)

        u0_lc = goto_attractor(ds, time_attract, integ_set)

        ds = init_ds_(rate_model, p, index_control, p2_loc,
        index_fix, var_fix, u0_lc, integ_set)
        
        ΛΛ = spectrum(ds, time_LSE)
        

        
        save_output(p2_loc_index, ΛΛ, u0_lc)
        save_tofile(namefile_LSE, namefile_u0s)

    end

end

function pre_broaching_with_print(sys, params, u0, time_setting, integrator_setting,
    index_fix, var_fix, index_control, range_p2)

end